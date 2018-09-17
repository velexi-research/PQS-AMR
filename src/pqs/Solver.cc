/*! \file Solver.cc
 *
 * \brief
 * Implementation file for Solver class.
 */

/*
 * ---------------------------------------------------------------------
 * COPYRIGHT/LICENSE.  This file is part of the XYZ package.  It is
 * subject to the license terms in the LICENSE file found in the
 * top-level directory of this distribution.  No part of the XYZ
 * package, including this file, may be copied, modified, propagated,
 * or distributed except according to the terms contained in the
 * LICENSE file.
 * ---------------------------------------------------------------------
 */

// --- Headers, namespaces, and type declarations

// Standard library
#include <cstddef>
#include <limits>
#include <math.h>
#include <memory>
#include <sstream>
#include <string>

// MPI
#include "mpi.h"

// SAMRAI
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/math/LSMAlgorithms.h"
#include "PQS/math/LSMToolbox.h"
#include "PQS/math/Toolbox.h"
#include "PQS/math/TimeIntegration.h"
#include "PQS/pqs/Algorithms.h"
#include "PQS/pqs/Solver.h"
#include "PQS/pqs/TagInitAndDataTransferModule.h"
#include "PQS/utilities/error.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class Patch; } }
namespace SAMRAI { namespace hier { class VariableContext; } }
namespace PQS { namespace pqs { class InterfaceInitStrategy; } }
namespace PQS { namespace pqs { class PoreInitStrategy; } }

// --- Class implementation

namespace PQS {
namespace pqs {

// --- Public methods

Solver::Solver(
        const shared_ptr<tbox::Database>& config_db,
        const shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy,
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy)
{
    // Check arguments
    if (config_db == NULL) {
        PQS_ERROR(this, "Solver", "'config_db' must not be NULL");
    }
    verifyConfigurationDatabase(config_db);

    if (pore_init_strategy == NULL) {
        PQS_ERROR(this, "Solver", "'pore_init_strategy' must not be NULL");
    }
    if (interface_init_strategy == NULL) {
        PQS_ERROR(this, "Solver", "'interface_init_strategy' must not be NULL");
    }

    // Set data members
    if (patch_hierarchy != NULL) {
        d_patch_hierarchy = patch_hierarchy;
    } else {
        createPatchHierarchy(config_db);
    }

    d_curvature = 0;  // TODO fix to be correct when restarting simulation
    d_step_count = 0;  // TODO fix to be correct when restarting simulation

    // Load configuration from config_db
    loadConfiguration(config_db);

    // Set up simulation variables
    setupSimulationVariables();

    // Construct pqs::Algorithms object
    shared_ptr<tbox::Database> pqs_config_db =
        config_db->getDatabase("PQS");
    d_pqs_algorithms = shared_ptr<pqs::Algorithms>(
            new pqs::Algorithms(pqs_config_db->getDatabase("Algorithms"),
                                d_lse_rhs_id, d_psi_id, d_grad_psi_id));

    // Construct math::LSM::Algorithms object
    d_lsm_algorithms = shared_ptr<PQS::math::LSM::Algorithms>(
            new PQS::math::LSM::Algorithms(d_patch_hierarchy,
                                           d_max_stencil_width));

    // Set up grid management objects
    setupGridManagement(config_db, pore_init_strategy, interface_init_strategy);

    // Initialize simulation
    initializeSimulation();

    // Emit contents of variable database to log file.
    tbox::plog << "\nVariable database..." << endl;
    hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

} // Solver::Solver()

Solver::~Solver()
{
    // Free memory allocated for simulation data
    for (int level_num = 0;
            level_num < d_patch_hierarchy->getNumberOfLevels(); level_num++) {
        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);
        patch_level->deallocatePatchData(d_permanent_variables);
        patch_level->deallocatePatchData(d_intermediate_variables);
    }
} // Solver::~Solver()

void Solver::resetHierarchyConfiguration(
        const int coarsest_level_num,
        const int finest_level_num)
{
    // --- Check arguments

    if (coarsest_level_num < 0) {
        PQS_ERROR(this, "resetHierarchyConfiguration",
                  "'coarsest_level_num' must be non-negative");
    }
    if (coarsest_level_num > finest_level_num) {
        PQS_ERROR(this, "resetHierarchyConfiguration",
                  string("'coarsest_level_num' must be less ") +
                  string("than or equal to 'finest_level_num'"));
    }
    for (int level_num = 0;
            level_num <= finest_level_num;
            level_num++) {

        if (d_patch_hierarchy->getPatchLevel(level_num) == NULL) {
            PQS_ERROR(this, "resetHierarchyConfiguration",
                      string("PatchLevel ") +
                      to_string(level_num) +
                      string(" is NULL in PatchHierarchy"));
        }
    }

    // --- Call resetHierarchyConfiguration() for LSM::Algorithms object

    d_lsm_algorithms->resetHierarchyConfiguration(coarsest_level_num,
                                                  finest_level_num);

} // Solver::resetHierarchyConfiguration()

void Solver::equilibrateInterface(
        const double curvature,
        const PQS_ALGORITHM_TYPE algorithm_type)
{
    // --- Check arguments

    if ( (algorithm_type != PRESCRIBED_CURVATURE_MODEL) &&
         (algorithm_type != SLIGHTLY_COMPRESSIBLE_MODEL) ) {
        PQS_ERROR(this, "equilibrateInterface",
                  string("Invalid PQS algorithm (") +
                  to_string(algorithm_type) +
                  string(")"));
    }

    // Get dimensionality of problem
    const int dim = d_patch_hierarchy->getDim().getValue();
    if ((dim != 2) && (dim != 3)) {
        PQS_ERROR(this, "equilibrateInterface",
                  string("Invalid number of spatial dimensions (=") +
                  to_string(dim) +
                  string("). Valid values: 2, 3."));
    }

    // --- Preparations

    // Create SAMRAI::math::HierarchyCellDataOpsReal object
    shared_ptr< SAMRAI::math::HierarchyCellDataOpsReal<PQS_REAL> > math_ops =
            shared_ptr< SAMRAI::math::HierarchyCellDataOpsReal<PQS_REAL> >(
                    new SAMRAI::math::HierarchyCellDataOpsReal<PQS_REAL>(
                             d_patch_hierarchy));

    // Compute volume of pore space
    const double pore_space_volume =
            math::LSM::computeVolume(
                    d_patch_hierarchy, d_psi_id,
                    -1, // compute volume for psi < 0
                    d_control_volume_id);

    // Initialize loop variables
    double t = 0.0;
    int step = 0;
    double previous_saturation = 0.0;

    double delta_phi = d_lsm_min_delta_phi + 1;
    double delta_saturation = d_lsm_min_delta_saturation + 1;

    // Allocate PatchData
    for (int level_num = 0;
            level_num < d_patch_hierarchy->getNumberOfLevels();
            level_num++) {

        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);

        patch_level->allocatePatchData(d_intermediate_variables);
    }

    // Copy phi data from LSM to PQS context
    math_ops->copyData(d_phi_lsm_current_id, d_phi_pqs_id);

    // --- Perform level set method computation

    while ( (step < d_lsm_max_iterations) &&
            (t < d_lsm_t_max) &&
            (delta_phi > d_lsm_min_delta_phi) &&
            (delta_saturation > d_lsm_min_delta_saturation) ) {

        // --- Preparations

        double dt;

        // Compute volume of non-wettting phase
        double non_wetting_phase_volume =
                math::LSM::computeVolume(
                        d_patch_hierarchy, d_phi_lsm_current_id,
                        -1, // compute volume for phi < 0
                        d_control_volume_id);

        // Use TVD Runge-Kutta integration in time to compute phi(t+dt)
        for (int rk_stage = 1; rk_stage <= d_time_integration_order; rk_stage++)
        {
            // --- Preparations

            // Time step variables
            double max_stable_dt_on_proc = std::numeric_limits<double>::max();

            // Configuration for computation of RHS of evolution equation
            int phi_id;
            if (d_time_integration_order == 1) {
                phi_id = d_phi_lsm_current_id;
            } else if (d_time_integration_order == 2) {
                if (rk_stage == 1) {
                    phi_id = d_phi_lsm_current_id;
                } else if (rk_stage == 2) {
                    phi_id = d_phi_lsm_next_id;
                }
            } else if (d_time_integration_order == 3) {
                if (rk_stage == 1) {
                    phi_id = d_phi_lsm_current_id;
                } else if (rk_stage == 2) {
                    phi_id = d_phi_lsm_next_id;
                } else if (rk_stage == 3) {
                    phi_id = d_phi_lsm_next_id;
                }
            }

            int ghost_cell_fill_context;
            if (rk_stage == 1) {
                ghost_cell_fill_context = LSM_CURRENT;
            } else {
                ghost_cell_fill_context = LSM_NEXT;
            }

            // --- Fill ghost cells

            d_tag_init_and_data_xfer_module->fillGhostCells(
                    ghost_cell_fill_context);

            // --- Compute RHS of level set evolution equation

            for (int level_num = 0;
                    level_num < d_patch_hierarchy->getNumberOfLevels();
                    level_num++) {

                shared_ptr<hier::PatchLevel> patch_level =
                    d_patch_hierarchy->getPatchLevel(level_num);

                // Compute RHS of level set evolution equation
                for (hier::PatchLevel::Iterator pi(patch_level->begin());
                        pi!=patch_level->end(); pi++) {

                    shared_ptr<hier::Patch> patch = *pi;

                    // Compute RHS of level set evolution equation on Patch
                    double stable_dt_on_patch;
                    if (algorithm_type == PRESCRIBED_CURVATURE_MODEL) {
                        stable_dt_on_patch = d_pqs_algorithms->
                                computePrescribedCurvatureModelRHS(
                                        patch, phi_id);
                    } else if (algorithm_type == SLIGHTLY_COMPRESSIBLE_MODEL) {
                        stable_dt_on_patch = d_pqs_algorithms->
                                computeSlightlyCompressibleModelRHS(
                                        patch, phi_id,
                                        non_wetting_phase_volume);
                    }

                    // Update maximum stable time step
                    if (stable_dt_on_patch < max_stable_dt_on_proc) {
                        max_stable_dt_on_proc = stable_dt_on_patch;
                    }
                }
            }

            // --- Compute stable time step (across processors)

            dt = max_stable_dt_on_proc;
            int err =
                tbox::SAMRAI_MPI::getSAMRAIWorld().AllReduce(&dt, 1, MPI_MIN);
            if (err != 0) {
                PQS_ERROR(this, "equilibrateInterface",
                          "MPI.AllReduce() to compute 'dt' failed");
            }

            // --- Advance phi

            if (d_time_integration_order == 1) {
                math::TimeIntegration::RK1Step(
                        d_patch_hierarchy,
                        d_phi_lsm_next_id,
                        d_phi_lsm_current_id,
                        d_lse_rhs_id,
                        dt);

            } else if (d_time_integration_order == 2) {
                if (rk_stage == 1) {
                    math::TimeIntegration::TVDRK2Stage1(
                            d_patch_hierarchy,
                            d_phi_lsm_next_id,
                            d_phi_lsm_current_id,
                            d_lse_rhs_id,
                            dt);

                } else if (rk_stage == 2) {
                    math::TimeIntegration::TVDRK2Stage2(
                            d_patch_hierarchy,
                            d_phi_lsm_next_id,
                            d_phi_lsm_next_id,
                            d_phi_lsm_current_id,
                            d_lse_rhs_id,
                            dt);
                }
            } else if (d_time_integration_order == 3) {
                if (rk_stage == 1) {
                    math::TimeIntegration::TVDRK3Stage1(
                            d_patch_hierarchy,
                            d_phi_lsm_next_id,
                            d_phi_lsm_current_id,
                            d_lse_rhs_id,
                            dt);

                } else if (rk_stage == 2) {
                    math::TimeIntegration::TVDRK3Stage2(
                            d_patch_hierarchy,
                            d_phi_lsm_next_id,
                            d_phi_lsm_next_id,
                            d_phi_lsm_current_id,
                            d_lse_rhs_id,
                            dt);

                } else if (rk_stage == 3) {
                    math::TimeIntegration::TVDRK3Stage3(
                            d_patch_hierarchy,
                            d_phi_lsm_next_id,
                            d_phi_lsm_next_id,
                            d_phi_lsm_current_id,
                            d_lse_rhs_id,
                            dt);
                }
            }
        }

        // --- Update metrics used in stopping criteria

        // Compute max norm of change in phi
        delta_phi = math::computeMaxNormDiff(d_patch_hierarchy,
                                             d_phi_lsm_next_id,
                                             d_phi_lsm_current_id,
                                             d_control_volume_id);

        // Compute change in saturation
        if (d_lsm_min_delta_saturation) {
            // Compute current saturation
            double saturation = non_wetting_phase_volume / pore_space_volume;

            // Compute change in saturation
            delta_saturation = fabs(saturation - previous_saturation);

            // Save previous saturation
            previous_saturation = saturation;
        }

        cout << step << ":" << dt << ":" << delta_phi << ", "
             << d_lsm_min_delta_phi << ", "
             << non_wetting_phase_volume << ", "
             << pore_space_volume << endl;

        // --- Prepare for next iteration

        // Update time and step
        t += dt;
        step++;

        // Swap phi data from LSM current and LSM next contexts
        math_ops->swapData(d_phi_lsm_next_id, d_phi_lsm_current_id);
    }

    // --- Clean up

    // Copy phi data from LSM to PQS context
    math_ops->copyData(d_phi_pqs_id, d_phi_lsm_current_id);

    // Deallocate PatchData
    for (int level_num = 0;
            level_num < d_patch_hierarchy->getNumberOfLevels();
            level_num++) {

        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);

        patch_level->deallocatePatchData(d_intermediate_variables);
    }
} // Solver::equilibrateInterface()

void Solver::reinitializeInterface(
        const math::LSM::REINIT_ALGORITHM_TYPE algorithm_type,
        const int max_time_steps,
        const double steady_state_condition,
        const double stop_distance)
{
    // --- Check arguments

    if ( (algorithm_type != math::LSM::REINIT_EQN_SGN_PHI0) &&
         (algorithm_type != math::LSM::REINIT_EQN_SGN_PHI) ) {
        PQS_ERROR(this, "reinitializeInterface",
                  string("Invalid reinitialization algorithm (") +
                  to_string(algorithm_type) +
                  string(")"));
    }

    // --- Use LSM::Algorithm object to reinitialize interface

    d_lsm_algorithms->reinitializeLevelSetFunction(
            d_phi_pqs_id,
            d_control_volume_id,
            d_time_integration_order,
            algorithm_type,
            max_time_steps,
            steady_state_condition,
            stop_distance);

} // Solver::reinitializeInterface()

double Solver::getCurvature() const
{
    return d_curvature;
} // Solver::getCurvature()

double Solver::getInitialCurvature() const
{
    return d_initial_curvature;
} // Solver::getInitialCurvature()

double Solver::getFinalCurvature() const
{
    return d_final_curvature;
} // Solver::getFinalCCurvature()

double Solver::getCurvatureStep() const
{
    return d_curvature_step;
} // Solver::getCurvatureStep()

int Solver::getStepCount() const
{
    return d_step_count;
} // Solver::getStepCount()

shared_ptr<hier::PatchHierarchy> Solver::getPatchHierarchy() const
{
    return d_patch_hierarchy;
} // Solver::getPatchHierarchy()

int Solver::getPoreSpacePatchDataId() const
{
    return d_psi_id;
} // Solver::getPoreSpacePatchDataId()

int Solver::getInterfacePatchDataId() const
{
    return d_phi_pqs_id;
} // Solver::getInterfacePatchDataId()

int Solver::getControlVolumePatchDataId() const
{
    return d_control_volume_id;
} // Solver::getControlVolumePatchDataId()

int Solver::getLSERHSPatchDataId() const
{
    return d_lse_rhs_id;
} // Solver::getLSERHSPatchDataId()

void Solver::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::Solver::printClassData..." << endl;
    os << "(Solver*) this = " << (Solver*) this << endl;
    os << "d_patch_hierarchy = " << d_patch_hierarchy.get() << endl;
    os << "d_gridding_algorithm = " << d_gridding_algorithm.get() << endl;
    os << "d_tag_init_and_data_xfer_module = "
       << d_tag_init_and_data_xfer_module.get() << endl;
    os << "d_pqs_algorithms = " << d_pqs_algorithms.get() << endl;

    os << endl;
    d_gridding_algorithm->printClassData(os);

    os << endl;
    d_tag_init_and_data_xfer_module->printClassData(os);

    os << endl;
    d_pqs_algorithms->printClassData(os);

} // Solver::printClassData()


// --- Private methods

void Solver::verifyConfigurationDatabase(
        const shared_ptr<tbox::Database>& config_db,
        const bool verify_patch_hierarchy) const
{
    // --- Check arguments

    if (config_db == NULL) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'config_db' must not be NULL");
    }

    // --- Verify PQS database

    if (!config_db->isDatabase("PQS")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'PQS' database missing from 'config_db'");
    }
    shared_ptr<tbox::Database> pqs_config_db =
        config_db->getDatabase("PQS");

    // Algorithms database
    if (!pqs_config_db->isDatabase("Algorithms")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'Algorithms' database missing from 'config_db'");
    }

    // Physical parameters
    if (!pqs_config_db->isDouble("initial_curvature")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'initial_curvature' missing from 'PQS' database");
    }
    if (!pqs_config_db->isDouble("final_curvature")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'final_curvature' missing from 'PQS' database");
    }
    if (!pqs_config_db->isDouble("curvature_step")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'curvature_step' missing from 'PQS' database");
    }

    // Level set method parameters
    if (!pqs_config_db->isDouble("lsm_t_max")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'lsm_t_max' missing from 'PQS' database");
    }
    if (!pqs_config_db->isInteger("lsm_max_iterations")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'lsm_max_iterations' missing from 'PQS' database");
    }
    if (!pqs_config_db->isDouble("lsm_min_delta_phi")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'lsm_min_delta_phi' missing from 'PQS' database");
    }
    if (!pqs_config_db->isDouble("lsm_min_delta_saturation")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'lsm_min_delta_saturation' missing from 'PQS' database");
    }

    // Numerical method parameters
    if (!pqs_config_db->isString("lsm_spatial_derivative_type")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'lsm_spatial_derivative_type' missing from 'PQS' database");
    }
    string lsm_spatial_derivative_type =
        pqs_config_db->getString("lsm_spatial_derivative_type");
    if ( !( (lsm_spatial_derivative_type == "ENO1") ||
            (lsm_spatial_derivative_type == "ENO2") ||
            (lsm_spatial_derivative_type == "ENO3") ||
            (lsm_spatial_derivative_type == "WENO5") ) ) {

        PQS_ERROR(this, "verifyConfigurationDatabase",
                  string("Invalid 'lsm_spatial_derivative_order'. ") +
                  string("Valid values: \"ENO1\", \"ENO2\", ") +
                  string("\"ENO3\", \"WENO5\"."));
    }

    if (!pqs_config_db->isInteger("time_integration_order")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'time_integration_order' missing from 'PQS' database");
    }
    int time_integration_order =
        pqs_config_db->getInteger("time_integration_order");
    if ( !( (time_integration_order == 1) ||
            (time_integration_order == 2) ||
            (time_integration_order == 3) ) ) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  string("Inavlid 'time_integration_order'. ") +
                  string("Valid values: 1, 2, 3."));
    }

    // Verify SAMRAI database
    if (!config_db->isDatabase("SAMRAI")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'SAMRAI' database missing from 'config_db'");
    }

    // Verify Geometry and PatchHierarchy databases
    if (verify_patch_hierarchy) {
        shared_ptr<tbox::Database> samrai_config_db =
            config_db->getDatabase("SAMRAI");
        if (!samrai_config_db->isDatabase("Geometry")) {
            PQS_ERROR(this, "verifyConfigurationDatabase",
                      "'Geometry' database missing from 'SAMRAI' database");
        }
        if (!samrai_config_db->isDatabase("PatchHierarchy")) {
            PQS_ERROR(this, "verifyConfigurationDatabase",
                      "'PatchHierarchy' database missing from 'SAMRAI' "
                      "database");
        }

        // Verify SAMRAI::Geometry database
        shared_ptr<tbox::Database> geometry_config_db =
            samrai_config_db->getDatabase("Geometry");
        if (!geometry_config_db->isInteger("dim")) {
            PQS_ERROR(this, "verifyConfigurationDatabase",
                      "'dim' database missing from 'SAMRAI::Geometry' "
                      "database");
        }
    }
}

void Solver::loadConfiguration(
        const shared_ptr<tbox::Database>& config_db)
{
    // --- Check arguments

    if (config_db == NULL) {
        PQS_ERROR(this, "loadConfiguration",
                  "'config_db' must not be NULL");
    }

    // --- Load configuration parameters

    shared_ptr<tbox::Database> pqs_config_db =
        config_db->getDatabase("PQS");

    // Physical parameters
    d_initial_curvature = pqs_config_db->getDouble("initial_curvature");
    d_final_curvature = pqs_config_db->getDouble("final_curvature");
    d_curvature_step = pqs_config_db->getDouble("curvature_step");

    // Level set method parameters
    d_lsm_t_max = pqs_config_db->getDouble("lsm_t_max");
    if (d_lsm_t_max <= 0) {
        PQS_ERROR(this, "loadConfiguration",
                  "'lsm_t_max' must be positive.");
    }

    d_lsm_max_iterations = pqs_config_db->getInteger("lsm_max_iterations");
    if (d_lsm_max_iterations <= 0) {
        PQS_ERROR(this, "loadConfiguration",
                  "'lsm_max_iterations' must be positive.");
    }

    d_lsm_min_delta_phi = pqs_config_db->getDouble("lsm_min_delta_phi");
    if (d_lsm_min_delta_phi < 0) {
        PQS_ERROR(this, "loadConfiguration",
                  "'lsm_min_delta_phi' must be non-negative.");
    }

    d_lsm_min_delta_saturation =
        pqs_config_db->getDouble("lsm_min_delta_saturation");
    if (d_lsm_min_delta_saturation < 0) {
        PQS_ERROR(this, "loadConfiguration",
                  "'lsm_min_delta_saturation' must be non-negative.");
    }

    // Numerical set method parameters
    string lsm_spatial_derivative_type =
        pqs_config_db->getString("lsm_spatial_derivative_type");
    if (lsm_spatial_derivative_type == "ENO1") {
        d_lsm_spatial_derivative_type = ENO1;
    } else if (lsm_spatial_derivative_type == "ENO2") {
        d_lsm_spatial_derivative_type = ENO2;
    } else if (lsm_spatial_derivative_type == "ENO3") {
        d_lsm_spatial_derivative_type = ENO3;
    } else if (lsm_spatial_derivative_type == "WENO5") {
        d_lsm_spatial_derivative_type = WENO5;
    }

    d_time_integration_order =
        pqs_config_db->getInteger("time_integration_order");

    // --- Set simulation parameters computed from configuration parameters

    // Set maximum stencil width
    int stencil_width;
    switch (d_lsm_spatial_derivative_type) {
        case ENO1: {
            stencil_width = 1;
            break;
        }
        case ENO2: {
            stencil_width = 2;
            break;
        }
        case ENO3: {
            stencil_width = 3;
            break;
        }
        case WENO5: {
            stencil_width = 3;
            break;
        }
        default: {
            PQS_ERROR(this, "setupGridManagement",
                      string("Invalid 'd_lsm_spatial_derivative_type' ") +
                      string("value: ") +
                      to_string(d_lsm_spatial_derivative_type));
        }
    }

    d_max_stencil_width = shared_ptr<hier::IntVector>(
        new hier::IntVector(d_patch_hierarchy->getDim(), stencil_width));

} // Solver::loadConfiguration()

void Solver::createPatchHierarchy(
        const shared_ptr<tbox::Database>& config_db)
{
    // Preparations
    shared_ptr<tbox::Database> samrai_config_db =
        config_db->getDatabase("SAMRAI");
    shared_ptr<tbox::Database> geometry_config_db =
        samrai_config_db->getDatabase("Geometry");

    // Create CartesianGridGeometry
    const tbox::Dimension dim(geometry_config_db->getInteger("dim"));
    shared_ptr<geom::CartesianGridGeometry> grid_geometry =
        shared_ptr<geom::CartesianGridGeometry>(
            new geom::CartesianGridGeometry(
                dim, "CartesianGeometry",
                samrai_config_db->getDatabase("Geometry")));

    // Create PatchHierarchy
    d_patch_hierarchy = shared_ptr<hier::PatchHierarchy>(
        new hier::PatchHierarchy("PatchHierarchy", grid_geometry,
            samrai_config_db->getDatabase("PatchHierarchy")));

} // Solver::createPatchHierarchy()

void Solver::setupSimulationVariables()
{
    // --- Preparations

    // Get dimensionality of problem
    tbox::Dimension dim = d_patch_hierarchy->getDim();

    // Create IntVector for zero ghost cell widths
    hier::IntVector zero_ghost_cell_width(dim, 0);

    // Initialize PatchData component selectors
    d_permanent_variables.clrAllFlags();
    d_intermediate_variables.clrAllFlags();

    // --- Create PatchData for simulation variables

    // Get pointer to VariableDatabase
    hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

    // Get variable contexts
    shared_ptr<hier::VariableContext> pqs_context =
        var_db->getContext("pqs");
    shared_ptr<hier::VariableContext> lsm_current_context =
        var_db->getContext("lsm_current");
    shared_ptr<hier::VariableContext> lsm_next_context =
        var_db->getContext("lsm_next");
    shared_ptr<hier::VariableContext> computed_context =
        var_db->getContext("computed");  // computed from phi or psi
    shared_ptr<hier::VariableContext> samr_context =
        var_db->getContext("SAMR");  // variable that supports SAMR

    // phi (fluid-fluid interface)
    shared_ptr< pdat::CellVariable<PQS_REAL> > phi_variable;
    if (var_db->checkVariableExists("phi")) {
        phi_variable =
            SAMRAI_SHARED_PTR_CAST< pdat::CellVariable<PQS_REAL> >(
                var_db->getVariable("phi"));
    } else {
        const int depth = 1;
        phi_variable = shared_ptr< pdat::CellVariable<PQS_REAL> >(
            new pdat::CellVariable<PQS_REAL>(dim, "phi", depth));
    }

    d_phi_pqs_id =
        var_db->registerVariableAndContext(phi_variable,
                                           pqs_context,
                                           zero_ghost_cell_width);
    d_permanent_variables.setFlag(d_phi_pqs_id);

    d_phi_lsm_current_id =
        var_db->registerVariableAndContext(phi_variable,
                                           lsm_current_context,
                                           *d_max_stencil_width);
    d_phi_lsm_next_id =
        var_db->registerVariableAndContext(phi_variable,
                                           lsm_next_context,
                                           *d_max_stencil_width);
    d_intermediate_variables.setFlag(d_phi_lsm_current_id);
    d_intermediate_variables.setFlag(d_phi_lsm_next_id);

    // psi (solid-pore interface)
    shared_ptr< pdat::CellVariable<PQS_REAL> > psi_variable;
    if (var_db->checkVariableExists("psi")) {
        psi_variable =
            SAMRAI_SHARED_PTR_CAST< pdat::CellVariable<PQS_REAL> >(
                var_db->getVariable("psi"));
    } else {
        const int depth = 1;
        psi_variable = shared_ptr< pdat::CellVariable<PQS_REAL> >(
            new pdat::CellVariable<PQS_REAL>(dim, "psi", depth));
    }
    d_psi_id =
        var_db->registerVariableAndContext(psi_variable,
                                           pqs_context,
                                           *d_max_stencil_width);
    d_permanent_variables.setFlag(d_psi_id);

    // grad psi (solid-pore interface)
    shared_ptr< pdat::CellVariable<PQS_REAL> > grad_psi_variable;
    if (var_db->checkVariableExists("grad psi")) {
        grad_psi_variable =
            SAMRAI_SHARED_PTR_CAST< pdat::CellVariable<PQS_REAL> >(
                var_db->getVariable("grad psi"));
    } else {
        const int depth = dim.getValue();
        grad_psi_variable =
            shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "grad psi", depth));
    }
    d_grad_psi_id =
        var_db->registerVariableAndContext(grad_psi_variable,
                                           computed_context,
                                           zero_ghost_cell_width);
    d_intermediate_variables.setFlag(d_grad_psi_id);

    /* TODO: review to see if these are needed

    // normal velocity
    shared_ptr< pdat::CellVariable<PQS_REAL> > normal_velocity_variable;
    if (var_db->checkVariableExists("normal velocity")) {
        normal_velocity_variable =
            SAMRAI_SHARED_PTR_CAST< pdat::CellVariable<PQS_REAL> >(
                var_db->getVariable("normal velocity"));
    } else {
        const int depth = 1;
        normal_velocity_variable =
            shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "normal velocity",
                                                    depth));
    }
    d_normal_velocity_id =
        var_db->registerVariableAndContext(normal_velocity_variable,
                                           computed_context,
                                           zero_ghost_cell_width);
    d_intermediate_variables.setFlag(d_normal_velocity_id);

    // vector velocity
    shared_ptr< pdat::CellVariable<PQS_REAL> > vector_velocity_variable;
    if (var_db->checkVariableExists("vector velocity")) {
        vector_velocity_variable =
            SAMRAI_SHARED_PTR_CAST< pdat::CellVariable<PQS_REAL> >(
                var_db->getVariable("vector velocity"));
    } else {
        const int depth = dim.getValue();
        vector_velocity_variable =
            shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "vector velocity",
                                                 depth));
    }
    d_vector_velocity_id =
        var_db->registerVariableAndContext(vector_velocity_variable,
                                           computed_context,
                                           zero_ghost_cell_width);
    d_intermediate_variables.setFlag(d_vector_velocity_id);

    */

    // RHS of level set evolution equation
    shared_ptr< pdat::CellVariable<PQS_REAL> > rhs_variable;
    if (var_db->checkVariableExists("rhs")) {
        rhs_variable =
            SAMRAI_SHARED_PTR_CAST< pdat::CellVariable<PQS_REAL> >(
                var_db->getVariable("rhs"));
    } else {
        const int depth = 1;
        rhs_variable = shared_ptr< pdat::CellVariable<PQS_REAL> >(
            new pdat::CellVariable<PQS_REAL>(dim, "rhs", depth));
    }
    d_lse_rhs_id =
        var_db->registerVariableAndContext(rhs_variable,
                                           lsm_current_context,
                                           zero_ghost_cell_width);
    d_intermediate_variables.setFlag(d_lse_rhs_id);

    // control volume
    shared_ptr< pdat::CellVariable<PQS_REAL> > control_volume_variable;
    if (var_db->checkVariableExists("control volume")) {
        control_volume_variable =
            SAMRAI_SHARED_PTR_CAST< pdat::CellVariable<PQS_REAL> >(
                var_db->getVariable("control volume"));
    } else {
        const int depth = 1;
        control_volume_variable =
            shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "control volume", depth));
    }
    d_control_volume_id =
        var_db->registerVariableAndContext(control_volume_variable,
                                           samr_context,
                                           zero_ghost_cell_width);

} // Solver::setupSimulationVariables()

void Solver::setupGridManagement(
        const shared_ptr<tbox::Database>& config_db,
        const shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy)
{
    // --- Check arguments

    // TODO: check database structure
    // PQS{ }
    // SAMRAI{ BoxGenerator, LoadBalancer, GriddingAlgorithm }
    //
    // throw runtime_error(
    // "'TagInitAndDataTransferModule' section not found in configuration database");

    // --- Preparations

    // Get SAMRAI configuration database
    shared_ptr<tbox::Database> samrai_config_db =
        config_db->getDatabase("SAMRAI");

    // --- Construct SAMRAI objects

    // Construct box generator and load balancer objects
    shared_ptr<mesh::BergerRigoutsos> box_generator =
        shared_ptr<mesh::BergerRigoutsos>(
            new mesh::BergerRigoutsos(
                d_patch_hierarchy->getDim(),
                samrai_config_db->getDatabase("BoxGenerator")));

    shared_ptr<mesh::ChopAndPackLoadBalancer> load_balancer =
        shared_ptr<mesh::ChopAndPackLoadBalancer> (
            new mesh::ChopAndPackLoadBalancer(
                d_patch_hierarchy->getDim(),
                "LoadBalancer",
                samrai_config_db->getDatabase("LoadBalancer")));

    // Construct PQS::pqs::TagInitAndDataTransferModule
    d_tag_init_and_data_xfer_module =
        shared_ptr<pqs::TagInitAndDataTransferModule>(
            new pqs::TagInitAndDataTransferModule(
                config_db->getDatabase("PQS"),
                d_patch_hierarchy,
                this,
                pore_init_strategy,
                interface_init_strategy,
                d_phi_pqs_id,
                d_phi_lsm_current_id,
                d_phi_lsm_next_id,
                d_psi_id,
                d_control_volume_id,
                *d_max_stencil_width));

    // Construct SAMRAI::mesh::GriddingAlgorithm object
    d_gridding_algorithm = shared_ptr<mesh::GriddingAlgorithm> (
        new mesh::GriddingAlgorithm(
            d_patch_hierarchy,
            "GriddingAlgorithm",
            samrai_config_db->getDatabase("GriddingAlgorithm"),
            d_tag_init_and_data_xfer_module, box_generator, load_balancer));

} // Solver::setupGridManagement()

void Solver::initializeSimulation()
{
    // Construct and initialize the levels of PatchHierarchy.
    if (tbox::RestartManager::getManager()->isFromRestart()) {
        // TODO
    } else {
        double time = 0.0;  // 0.0 is an arbitrary simulation time

        d_gridding_algorithm->makeCoarsestLevel(time);

        for (int level_num = 0;
                d_patch_hierarchy->levelCanBeRefined(level_num);
                level_num++) {

            d_gridding_algorithm->makeFinerLevel(
                    0, // TODO: do we need tag buffer?
                    true, // initial_cycle=true
                    d_step_count,
                    time);
        }

        // TODO: synchronize coarser levels with finer levels that didn't
        // exist when the finer coarser level data was initialized.
    }
} // Solver::initializeSimulation()

} // PQS::pqs namespace
} // PQS namespace
