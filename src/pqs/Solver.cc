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
#include <sstream>
#include <string>

// Boost
#include <boost/smart_ptr/make_shared_object.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/utilities.h"
#include "PQS/pqs/Algorithms.h"
#include "PQS/pqs/Solver.h"
#include "PQS/pqs/TagInitAndDataTransferModule.h"

// Class/type declarations
namespace SAMRAI { namespace mesh { class TagAndInitializeStrategy; } }
namespace SAMRAI { namespace hier { class Patch; } }
namespace SAMRAI { namespace hier { class Variable; } }
namespace SAMRAI { namespace hier { class VariableContext; } }
namespace PQS { namespace pqs { class InterfaceInitStrategy; } }
namespace PQS { namespace pqs { class PoreInitStrategy; } }

// --- Class implementation

namespace PQS {
namespace pqs {

// --- Implementation of public methods

Solver::Solver(
        const boost::shared_ptr<tbox::Database>& config_db,
        const boost::shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const boost::shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy,
        const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy)
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
    d_cycle = 0;

    // Load configuration from config_db
    loadConfiguration(config_db);

    // Set up simulation variables
    setupSimulationVariables();

    // Set up grid management objects
    setupGridManagement(config_db, pore_init_strategy, interface_init_strategy);

    // Initialize simulation
    initializeSimulation();

} // Solver::Solver()

Solver::~Solver()
{
    // Free memory allocated for simulation data
    for (int level_num=0;
            level_num < d_patch_hierarchy->getNumberOfLevels(); level_num++) {
        boost::shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);
        patch_level->deallocatePatchData(d_permanent_variables);
        patch_level->deallocatePatchData(d_intermediate_variables);
    }
} // Solver::~Solver()

void Solver::equilibrateInterface(const double curvature)
{
    // --- Preparations

    // Iinitialize loop variables
    double t = 0.0;
    double dt;

    int step = 0;

    double delta_phi = 2 * d_lsm_min_delta_phi;
    double delta_saturation = 2 * d_lsm_min_delta_saturation;

    // --- Perform level set method computation

    while ( (step < d_lsm_max_iterations) &&
            (t < d_lsm_t_max) &&
            (delta_phi < d_lsm_min_delta_phi) &&
            (delta_saturation < d_lsm_min_delta_saturation) ) {

        // --- Fill ghost cells

        d_tag_init_and_data_xfer_module->fillGhostCells();

        // --- Compute RHS of level set evolution equation

        // Loop over PatchLevel in PatchHierarchy
        for (int level_num=0;
                level_num < d_patch_hierarchy->getNumberOfLevels();
                level_num++) {

            boost::shared_ptr<hier::PatchLevel> patch_level =
                d_patch_hierarchy->getPatchLevel(level_num);

            // Loop over Patches on PatchLevel
            for (hier::PatchLevel::Iterator pi(patch_level->begin());
                    pi!=patch_level->end(); pi++) {

                boost::shared_ptr<hier::Patch> patch = *pi;

                // Compute RHS of level set evolution equation on Patch
                Algorithms::computeSlightlyCompressibleModelRHS(patch);
            }
        }

        // --- Compute stable time step

        // TODO: dt =

        // --- Advance phi

        // TODO

        // --- Computing stopping criteria

        // TODO
    }
} // Solver::equilibrateInterface()

void Solver::advanceInterface(const double delta_curvature)
{
    // TODO
    bool done = true;
    while (!done) {
        // Loop over PatchLevel in PatchHierarchy
        for (int level_num=0;
                level_num < d_patch_hierarchy->getNumberOfLevels();
                level_num++) {

            boost::shared_ptr<hier::PatchLevel> patch_level =
                d_patch_hierarchy->getPatchLevel(level_num);

            // Loop over Patches on PatchLevel
            for (hier::PatchLevel::Iterator pi(patch_level->begin());
                    pi!=patch_level->end(); pi++) {

                boost::shared_ptr<hier::Patch> patch = *pi;
                Algorithms::computePrescribedCurvatureModelRHS(patch);
            }
        }

        // * compute RHS
        // * advance phi
        // * computing stopping criteria
    }
} // Solver::advanceInterface()

double Solver::getCurvature() const
{
    return d_curvature;
} // Solver::getCurvature()

int Solver::getCycle() const
{
    return d_cycle;
} // Solver::getCycle()

boost::shared_ptr<hier::PatchHierarchy> Solver::getPatchHierarchy() const
{
    return d_patch_hierarchy;
} // Solver::getPatchHierarchy()

int Solver::getPoreSpacePatchDataId() const
{
    return d_psi_id;
} // Solver::getPoreSpacePatchDataId()

int Solver::getInterfacePatchDataId() const
{
    return d_phi_pqs_current_id;
} // Solver::getInterfacePatchDataId()


void Solver::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::Solver::printClassData..." << endl;
    os << "(Solver*) this = " << (Solver*) this << endl;
    os << "d_patch_hierarchy = " << d_patch_hierarchy.get() << endl;
    os << "d_gridding_alg = " << d_gridding_alg.get() << endl;
    os << "d_tag_init_and_data_xfer_module = "
       << d_tag_init_and_data_xfer_module.get() << endl;

    os << endl;
    d_gridding_alg->printClassData(os);

    os << endl;
    d_tag_init_and_data_xfer_module->printClassData(os);

} // Solver::printClassData()

// --- Implementation of private methods

void Solver::verifyConfigurationDatabase(
        const boost::shared_ptr<tbox::Database>& config_db,
        const bool verify_patch_hierarchy) const
{
    if (config_db == NULL) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'config_db' must not be NULL");
    }

    // --- Verify PQS database

    if (!config_db->isDatabase("PQS")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'PQS' database missing from 'config_db'");
    }
    boost::shared_ptr<tbox::Database> pqs_config_db =
        config_db->getDatabase("PQS");

    // Physical parameters
    if (!pqs_config_db->isDouble("initial_curvature")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'initial_curvature' missing from 'PQS' database");
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
    std::string lsm_spatial_derivative_type =
        pqs_config_db->getString("lsm_spatial_derivative_type");
    if ( !( (lsm_spatial_derivative_type == "ENO1") ||
            (lsm_spatial_derivative_type == "ENO2") ||
            (lsm_spatial_derivative_type == "ENO3") ||
            (lsm_spatial_derivative_type == "WENO5") ) ) {

        PQS_ERROR(this, "verifyConfigurationDatabase",
                  std::string("Invalid 'lsm_spatial_derivative_order'. ") +
                  std::string("Valid values: \"ENO1\", \"ENO2\", ") +
                  std::string("\"ENO3\", \"WENO5\"."));
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
                  std::string("Inavlid 'time_integration_order'. ") +
                  std::string("Valid values: 1, 2, 3."));
    }

    // Verify SAMRAI database
    if (!config_db->isDatabase("SAMRAI")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'SAMRAI' database missing from 'config_db'");
    }

    // Verify Geometry and PatchHierarchy databases
    if (verify_patch_hierarchy) {
        boost::shared_ptr<tbox::Database> samrai_config_db =
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
        boost::shared_ptr<tbox::Database> geometry_config_db =
            samrai_config_db->getDatabase("Geometry");
        if (!geometry_config_db->isInteger("dim")) {
            PQS_ERROR(this, "verifyConfigurationDatabase",
                      "'dim' database missing from 'SAMRAI::Geometry' "
                      "database");
        }
    }
}

void Solver::loadConfiguration(
        const boost::shared_ptr<tbox::Database>& config_db)
{
    // --- Load configuration parameters

    boost::shared_ptr<tbox::Database> pqs_config_db =
        config_db->getDatabase("PQS");

    // Physical parameters
    d_curvature = pqs_config_db->getDouble("initial_curvature");

    // Level set method parameters
    d_lsm_t_max = pqs_config_db->getDouble("lsm_t_max");
    d_lsm_max_iterations = pqs_config_db->getInteger("lsm_max_iterations");
    d_lsm_min_delta_phi = pqs_config_db->getDouble("lsm_min_delta_phi");
    d_lsm_min_delta_saturation =
        pqs_config_db->getDouble("lsm_min_delta_saturation");

    // Numerical set method parameters
    std::string lsm_spatial_derivative_type =
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

    // Set maximum ghost cell width
    int ghost_cell_width;
    switch (d_lsm_spatial_derivative_type) {
        case ENO1: {
            ghost_cell_width = 1;
            break;
        }
        case ENO2: {
            ghost_cell_width = 2;
            break;
        }
        case ENO3: {
            ghost_cell_width = 3;
            break;
        }
        case WENO5: {
            ghost_cell_width = 3;
            break;
        }
        default: {
            PQS_ERROR(this, "setupGridManagement",
                      std::string("Invalid 'd_lsm_spatial_derivative_type' ") +
                      std::string("value: ") +
                      std::to_string(d_lsm_spatial_derivative_type));
        }
    }

    d_max_ghost_cell_width = boost::shared_ptr<hier::IntVector>(
        new hier::IntVector(d_patch_hierarchy->getDim(), ghost_cell_width));

} // Solver::loadConfiguration()

void Solver::createPatchHierarchy(
        const boost::shared_ptr<tbox::Database>& config_db)
{
    // Preparations
    boost::shared_ptr<tbox::Database> samrai_config_db =
        config_db->getDatabase("SAMRAI");
    boost::shared_ptr<tbox::Database> geometry_config_db =
        samrai_config_db->getDatabase("Geometry");

    // Create CartesianGridGeometry
    const tbox::Dimension dim(geometry_config_db->getInteger("dim"));
    boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry =
        boost::shared_ptr<geom::CartesianGridGeometry>(
            new geom::CartesianGridGeometry(
                dim, "CartesianGeometry",
                samrai_config_db->getDatabase("Geometry")));

    // Create PatchHierarchy
    d_patch_hierarchy = boost::shared_ptr<hier::PatchHierarchy>(
        new hier::PatchHierarchy("PatchHierarchy", grid_geometry,
            samrai_config_db->getDatabase("PatchHierarchy")));

} // Solver::createPatchHierarchy()

void Solver::setupSimulationVariables()
{
    // --- Preparations

    // Get dimensionality of problem
    tbox::Dimension dim = d_patch_hierarchy->getDim();

    // Create IntVector for ghost cell widths
    hier::IntVector max_ghost_cell_width(*d_max_ghost_cell_width);
    hier::IntVector zero_ghost_cell_width(dim, 0);

    // Initialize PatchData component selectors
    d_permanent_variables.clrAllFlags();
    d_intermediate_variables.clrAllFlags();

    // --- Create PatchData for simulation variables

    // Get pointer to VariableDatabase
    hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

    // Get default variable context
    boost::shared_ptr<hier::VariableContext> default_context =
        var_db->getContext("default");

    // phi (fluid-fluid interface)
    boost::shared_ptr< pdat::CellVariable<PQS_REAL> > phi_variable;
    if (var_db->checkVariableExists("phi")) {
        phi_variable =
            BOOST_CAST<pdat::CellVariable<PQS_REAL>, hier::Variable>(
                var_db->getVariable("phi"));
    } else {
        const int depth = 1;
        phi_variable = boost::shared_ptr< pdat::CellVariable<PQS_REAL> >(
            new pdat::CellVariable<PQS_REAL>(dim, "phi", depth));
    }
    boost::shared_ptr<hier::VariableContext> pqs_current_context =
        var_db->getContext("pqs_current");
    boost::shared_ptr<hier::VariableContext> pqs_next_context =
        var_db->getContext("pqs_next");
    boost::shared_ptr<hier::VariableContext> lsm_current_context =
        var_db->getContext("lsm_current");
    boost::shared_ptr<hier::VariableContext> lsm_next_context =
        var_db->getContext("lsm_next");

    d_phi_pqs_current_id =
        var_db->registerVariableAndContext(phi_variable,
                                           pqs_current_context,
                                           zero_ghost_cell_width);
    d_phi_pqs_next_id =
        var_db->registerVariableAndContext(phi_variable,
                                           pqs_next_context,
                                           zero_ghost_cell_width);
    d_permanent_variables.setFlag(d_phi_pqs_current_id);
    d_intermediate_variables.setFlag(d_phi_pqs_next_id);

    d_phi_lsm_current_id =
        var_db->registerVariableAndContext(phi_variable,
                                           lsm_current_context,
                                           max_ghost_cell_width);
    d_phi_lsm_next_id =
        var_db->registerVariableAndContext(phi_variable,
                                           lsm_next_context,
                                           max_ghost_cell_width);
    d_intermediate_variables.setFlag(d_phi_lsm_current_id);
    d_intermediate_variables.setFlag(d_phi_lsm_next_id);

    // psi (solid-pore interface)
    boost::shared_ptr< pdat::CellVariable<PQS_REAL> > psi_variable;
    if (var_db->checkVariableExists("psi")) {
        psi_variable =
            BOOST_CAST<pdat::CellVariable<PQS_REAL>, hier::Variable>(
                var_db->getVariable("psi"));
    } else {
        const int depth = 1;
        psi_variable = boost::shared_ptr< pdat::CellVariable<PQS_REAL> >(
            new pdat::CellVariable<PQS_REAL>(dim, "psi", depth));
    }
    d_psi_id =
        var_db->registerVariableAndContext(psi_variable,
                                           default_context,
                                           max_ghost_cell_width);
    d_permanent_variables.setFlag(d_psi_id);

    // grad psi (solid-pore interface)
    boost::shared_ptr< pdat::CellVariable<PQS_REAL> > grad_psi_variable;
    if (var_db->checkVariableExists("grad psi")) {
        grad_psi_variable =
            BOOST_CAST<pdat::CellVariable<PQS_REAL>, hier::Variable>(
                var_db->getVariable("grad psi"));
    } else {
        const int depth = dim.getValue();
        grad_psi_variable =
            boost::shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "grad psi", depth));
    }
    d_grad_psi_id =
        var_db->registerVariableAndContext(grad_psi_variable,
                                           default_context,
                                           zero_ghost_cell_width);
    d_intermediate_variables.setFlag(d_grad_psi_id);

    /* TODO: review to see if these are needed

    // normal velocity
    boost::shared_ptr< pdat::CellVariable<PQS_REAL> > normal_velocity_variable;
    if (var_db->checkVariableExists("normal velocity")) {
        normal_velocity_variable =
            BOOST_CAST<pdat::CellVariable<PQS_REAL>, hier::Variable>(
                var_db->getVariable("normal velocity"));
    } else {
        const int depth = 1;
        normal_velocity_variable =
            boost::shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "normal velocity",
                                                    depth));
    }
    d_normal_velocity_id =
        var_db->registerVariableAndContext(normal_velocity_variable,
                                           default_context,
                                           zero_ghost_cell_width);
    d_intermediate_variables.setFlag(d_normal_velocity_id);

    // vector velocity
    boost::shared_ptr< pdat::CellVariable<PQS_REAL> > vector_velocity_variable;
    if (var_db->checkVariableExists("vector velocity")) {
        vector_velocity_variable =
            BOOST_CAST<pdat::CellVariable<PQS_REAL>, hier::Variable>(
                var_db->getVariable("vector velocity"));
    } else {
        const int depth = dim.getValue();
        vector_velocity_variable =
            boost::shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "vector velocity",
                                                 depth));
    }
    d_vector_velocity_id =
        var_db->registerVariableAndContext(vector_velocity_variable,
                                           default_context,
                                           zero_ghost_cell_width);
    d_intermediate_variables.setFlag(d_vector_velocity_id);

    */

    // RHS of level set evolution equation
    boost::shared_ptr< pdat::CellVariable<PQS_REAL> > lse_rhs_variable;
    if (var_db->checkVariableExists("LSE RHS")) {
        lse_rhs_variable =
            BOOST_CAST<pdat::CellVariable<PQS_REAL>, hier::Variable>(
                var_db->getVariable("LSE RHS"));
    } else {
        const int depth = 1;
        lse_rhs_variable = boost::shared_ptr< pdat::CellVariable<PQS_REAL> >(
            new pdat::CellVariable<PQS_REAL>(dim, "LSE RHS", depth));
    }
    d_lse_rhs_id =
        var_db->registerVariableAndContext(lse_rhs_variable,
                                           default_context,
                                           zero_ghost_cell_width);
    d_intermediate_variables.setFlag(d_psi_id);

    // control volume
    boost::shared_ptr< pdat::CellVariable<PQS_REAL> > control_volume_variable;
    if (var_db->checkVariableExists("control volume")) {
        control_volume_variable =
            BOOST_CAST<pdat::CellVariable<PQS_REAL>, hier::Variable>(
                var_db->getVariable("control volume"));
    } else {
        const int depth = 1;
        control_volume_variable =
            boost::shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "control volume", depth));
    }
    d_control_volume_id =
        var_db->registerVariableAndContext(control_volume_variable,
                                           default_context,
                                           zero_ghost_cell_width);
    d_intermediate_variables.setFlag(d_control_volume_id);

} // Solver::setupSimulationVariables()

void Solver::setupGridManagement(
        const boost::shared_ptr<tbox::Database>& config_db,
        const boost::shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const boost::shared_ptr<pqs::InterfaceInitStrategy>&
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
    boost::shared_ptr<tbox::Database> samrai_config_db =
        config_db->getDatabase("SAMRAI");

    // --- Construct SAMRAI objects

    // Construct box generator and load balancer objects
    boost::shared_ptr<mesh::BergerRigoutsos> box_generator =
        boost::shared_ptr<mesh::BergerRigoutsos>(
            new mesh::BergerRigoutsos(
                d_patch_hierarchy->getDim(),
                samrai_config_db->getDatabase("BoxGenerator")));

    boost::shared_ptr<mesh::ChopAndPackLoadBalancer> load_balancer =
        boost::shared_ptr<mesh::ChopAndPackLoadBalancer> (
            new mesh::ChopAndPackLoadBalancer(
                d_patch_hierarchy->getDim(),
                "LoadBalancer",
                samrai_config_db->getDatabase("LoadBalancer")));

    // Construct PQS::pqs::TagInitAndDataTransferModule
    d_tag_init_and_data_xfer_module =
        boost::shared_ptr<pqs::TagInitAndDataTransferModule>(
            new pqs::TagInitAndDataTransferModule(
                config_db->getDatabase("PQS"),
                d_patch_hierarchy,
                pore_init_strategy,
                interface_init_strategy,
                d_phi_pqs_current_id, d_psi_id));

    // Construct SAMRAI::mesh::GriddingAlgorithm object
    d_gridding_alg = boost::shared_ptr<mesh::GriddingAlgorithm> (
        new mesh::GriddingAlgorithm(
            d_patch_hierarchy,
            "GriddingAlgorithm",
            samrai_config_db->getDatabase("GriddingAlgorithm"),
            boost::shared_ptr<mesh::TagAndInitializeStrategy>(
                d_tag_init_and_data_xfer_module),
            box_generator,
            load_balancer));

} // Solver::setupGridManagement()

void Solver::initializeSimulation()
{
    // Construct and initialize the levels of PatchHierarchy.
    if (tbox::RestartManager::getManager()->isFromRestart()) {
        // TODO
    } else {
        double time = 0.0;  // 0.0 is an arbitrary simulation time

        d_gridding_alg->makeCoarsestLevel(time);

        for (int level_num = 0;
                d_patch_hierarchy->levelCanBeRefined(level_num);
                level_num++) {

            d_gridding_alg->makeFinerLevel(0, // TODO: do we need tag buffer?
                                           true, // initial_cycle=true
                                           d_cycle,
                                           time);
        }

        // TODO: synchronize coarser levels with finer levels that didn't
        // exist when the finer coarser level data was initialized.
    }
} // Solver::initializeSimulation()

} // PQS::pqs namespace
} // PQS namespace
