/*! \file LSMAlgorithms.cc
 *
 * \brief
 * Implementation file for LSM::Algorithms class.
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
#include <memory>
#include <sstream>

// SAMRAI
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/PIO.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/math/LSMAlgorithms.h"
#include "PQS/utilities/error.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class IntVector; } }

// --- Class implementation

namespace PQS {
namespace math {
namespace LSM {

// --- Public methods

Algorithms::Algorithms(
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const shared_ptr<hier::IntVector>& max_stencil_width)
{
    // Check arguments
    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "Algorithms", "'patch_hierarchy' must not be NULL");
    }

    // Set data members
    d_patch_hierarchy = patch_hierarchy;

    // Set up simulation variables
    setupSimulationVariables();

    // Emit contents of variable database to log file.
    tbox::plog << "\nVariable database..." << endl;
    hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

} // Algorithms::Algorithms()

Algorithms::~Algorithms()
{
    // Free memory allocated for simulation data
    for (int level_num = 0;
            level_num < d_patch_hierarchy->getNumberOfLevels(); level_num++) {
        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);
        patch_level->deallocatePatchData(d_lsm_algorithm_variables);
    }
} // Algorithms::~Algorithms()

void Algorithms::reinitializeLevelSetFunction(
        const shared_ptr<hier::PatchHierarchy> patch_hierarchy,
        const int phi_id,
        const int max_time_steps,
        const double steady_state_condition,
        const double stop_distance)
{
/*
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
            math::LSM::computeVolume(d_patch_hierarchy,
                                     d_psi_id,
                                     -1, // compute volume for psi < 0
                                     d_control_volume_id);

    // Initialize loop variables
    double t = 0.0;
    int step = 0;
    double previous_saturation = 0.0;

    double delta_phi = 2 * d_lsm_min_delta_phi;
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
                math::LSM::computeVolume(d_patch_hierarchy,
                                        d_phi_lsm_current_id,
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
*/
} // Algorithms::reinitializeLevelSetFunctions()

void Algorithms::computeExtensionField(
        const shared_ptr<hier::PatchHierarchy> patch_hierarchy,
        const int phi_id,
        const int S_id,
        const int max_time_steps,
        const double steady_state_condition,
        const double stop_distance)
{
} // Algorithms::computeExtensionField()

shared_ptr<hier::PatchHierarchy> Algorithms::getPatchHierarchy() const
{
    return d_patch_hierarchy;
} // Algorithms::getPatchHierarchy()

void Algorithms::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::math::Algorithms::printClassData..." << endl;
    os << "(Algorithms*) this = " << (Algorithms*) this << endl;
    os << "d_patch_hierarchy = " << d_patch_hierarchy.get() << endl;

} // Algorithms::printClassData()


// --- Private methods

void Algorithms::setupSimulationVariables()
{
/*
    // --- Preparations

    // Get dimensionality of problem
    tbox::Dimension dim = d_patch_hierarchy->getDim();

    // Create IntVector for ghost cell widths
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

    // RHS of level set evolution equation
    shared_ptr< pdat::CellVariable<PQS_REAL> > lse_rhs_variable;
    if (var_db->checkVariableExists("LSE RHS")) {
        lse_rhs_variable =
            SAMRAI_SHARED_PTR_CAST< pdat::CellVariable<PQS_REAL> >(
                var_db->getVariable("LSE RHS"));
    } else {
        const int depth = 1;
        lse_rhs_variable = shared_ptr< pdat::CellVariable<PQS_REAL> >(
            new pdat::CellVariable<PQS_REAL>(dim, "LSE RHS", depth));
    }
    d_lse_rhs_id =
        var_db->registerVariableAndContext(lse_rhs_variable,
                                           computed_context,
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
*/
} // Algorithms::setupSimulationVariables()

} // PQS::math::LSM namespace
} // PQS::math namespace
} // PQS namespace
