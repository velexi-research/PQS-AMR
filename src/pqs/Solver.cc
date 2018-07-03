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
#include <stdexcept>

// Boost
#include <boost/smart_ptr/make_shared_object.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/utilities.h"
#include "PQS/pqs/Solver.h"
#include "PQS/pqs/TagAndInitModule.h"

// Class/type declarations
namespace SAMRAI { namespace mesh { class TagAndInitializeStrategy; } }
namespace SAMRAI { namespace hier { class Variable; } }
namespace SAMRAI { namespace hier { class VariableContext; } }
namespace PQS { namespace pqs { class DataInitStrategy; } }


// --- Class implementation

namespace PQS {
namespace pqs {

// --- Implementation of public methods

Solver::Solver(
        const boost::shared_ptr<tbox::Database>& config_db,
        const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const boost::shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const boost::shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy)
{
    // Check parameters
    if (config_db == NULL) {
        PQS_ERROR(this, "Solver", "'config_db' must not be NULL");
    }
    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "Solver", "'patch_hierarchy' must not be NULL");
    }
    if (pore_init_strategy == NULL) {
        PQS_ERROR(this, "Solver", "'pore_init_strategy' must not be NULL");
    }
    if (interface_init_strategy == NULL) {
        PQS_ERROR(this, "Solver", "'interface_init_strategy' must not be NULL");
    }

    // Set data members
    d_patch_hierarchy = patch_hierarchy;
    d_num_steps = 0;

    // Load configuration from config_db
    loadConfiguration(config_db);

    // Set up simulation variables
    setupSimulationVariables();

    // Set up grid management objects
    setupGridManagement(config_db, pore_init_strategy, interface_init_strategy);

    // Initialize simulation
    initializeSimulation();

} // Solver::Solver()

void Solver::equilibrateInterface(const double curvature)
{
    // TODO
} // Solver::equilibrateInterface()

void Solver::advanceInterface(const double delta_curvature)
{
    // TODO
} // Solver::advanceInterface()

double Solver::getCurvature() const
{
    return d_curvature;
} // Solver::getCurvature()

int Solver::getStep() const
{
    return d_num_steps;
} // Solver::getStep()

void Solver::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::Solver::printClassData..." << endl;
    os << "(Solver*) this = " << (Solver*) this << endl;
    os << "d_patch_hierarchy = " << d_patch_hierarchy.get() << endl;
    os << "d_gridding_alg = " << d_gridding_alg.get() << endl;
    os << "d_tag_and_init_module = " << d_tag_and_init_module.get() << endl;

    os << endl;
    d_gridding_alg->printClassData(os);

    os << endl;
    d_tag_and_init_module->printClassData(os);

} // Solver::printClassData()

// --- Implementation of private methods

void Solver::loadConfiguration(
        const boost::shared_ptr<tbox::Database>& config_db)
{
    // Check parameters
    if (config_db == NULL) {
        PQS_ERROR(this, "loadConfiguration", "'config_db' must not be NULL");
    }
    if (!config_db->isDatabase("PQS")) {
        PQS_ERROR(this, "loadConfiguration",
                  "'PQS' database missing from 'config_db'");
    }

    // Load configuration parameters
    boost::shared_ptr<tbox::Database> pqs_config_db =
        config_db->getDatabase("PQS");
    d_curvature = pqs_config_db->getDouble("initial_curvature");

} // Solver::loadConfiguration()

void Solver::setupSimulationVariables()
{
    // Initialize PatchData component selectors
    d_permanent_variables.clrAllFlags();
    d_intermediate_variables.clrAllFlags();

    // --- Create PatchData for simulation variables

    // Get dimensionality of problem
    tbox::Dimension dim = d_patch_hierarchy->getDim();

    // Get pointer to VariableDatabase
    hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

    // Create zero ghostcell width IntVector
    hier::IntVector zero_ghostcell_width(dim, 0);

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
                                           zero_ghostcell_width);
    d_phi_pqs_next_id =
        var_db->registerVariableAndContext(phi_variable,
                                           pqs_next_context,
                                           zero_ghostcell_width);
    d_permanent_variables.setFlag(d_phi_pqs_current_id);
    d_intermediate_variables.setFlag(d_phi_pqs_next_id);

    d_phi_lsm_current_id =
        var_db->registerVariableAndContext(phi_variable,
                                           lsm_current_context,
                                           zero_ghostcell_width);
    d_phi_lsm_next_id =
        var_db->registerVariableAndContext(phi_variable,
                                           lsm_next_context,
                                           zero_ghostcell_width);
    d_permanent_variables.setFlag(d_phi_lsm_current_id);
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
                                           zero_ghostcell_width);
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
                                           zero_ghostcell_width);
    d_intermediate_variables.setFlag(d_grad_psi_id);

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
                                           zero_ghostcell_width);
    d_intermediate_variables.setFlag(d_normal_velocity_id);

    // external velocity
    boost::shared_ptr< pdat::CellVariable<PQS_REAL> >
        external_velocity_variable;
    if (var_db->checkVariableExists("external velocity")) {
        external_velocity_variable =
            BOOST_CAST<pdat::CellVariable<PQS_REAL>, hier::Variable>(
                var_db->getVariable("external velocity"));
    } else {
        const int depth = dim.getValue();
        external_velocity_variable =
            boost::shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "external velocity",
                                                    depth));
    }
    d_external_velocity_id =
        var_db->registerVariableAndContext(external_velocity_variable,
                                           default_context,
                                           zero_ghostcell_width);
    d_intermediate_variables.setFlag(d_external_velocity_id);

    // curvature
    boost::shared_ptr< pdat::CellVariable<PQS_REAL> > curvature_variable;
    if (var_db->checkVariableExists("curvature")) {
        curvature_variable =
            BOOST_CAST<pdat::CellVariable<PQS_REAL>, hier::Variable>(
                var_db->getVariable("curvature"));
    } else {
        const int depth = 1;
        curvature_variable =
            boost::shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "curvature", depth));
    }
    d_curvature_id =
        var_db->registerVariableAndContext(curvature_variable,
                                           default_context,
                                           zero_ghostcell_width);
    d_intermediate_variables.setFlag(d_curvature_id);

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
                                           zero_ghostcell_width);
    d_intermediate_variables.setFlag(d_control_volume_id);

} // Solver::setupSimulationVariables()

void Solver::setupGridManagement(
        const boost::shared_ptr<tbox::Database>& config_db,
        const boost::shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const boost::shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy)
{
    // --- Check parameters
    // TODO: check database structure
    // PQS{ }
    // SAMRAI{ BoxGenerator, LoadBalancer, GriddingAlgorithm }
    //
    // throw runtime_error(
    // "'TagAndInitModule' section not found in configuration database");

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
            new mesh::ChopAndPackLoadBalancer(d_patch_hierarchy->getDim(),
                "LoadBalancer",
                samrai_config_db->getDatabase("LoadBalancer")));

    // Construct PQS::pqs::TagAndInitModule
    boost::shared_ptr<pqs::TagAndInitModule> d_tag_and_init_module =
        boost::shared_ptr<pqs::TagAndInitModule>(
            new pqs::TagAndInitModule(config_db->getDatabase("PQS"),
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
                d_tag_and_init_module),
            box_generator,
            load_balancer));

} // Solver::setupGridManagement()

void Solver::initializeSimulation()
{
    // TODO
} // Solver::initializeSimulation()

} // PQS::pqs namespace
} // PQS namespace
