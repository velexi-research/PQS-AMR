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

// Standard library headers
#include <cstddef>
#include <sstream>
#include <stdexcept>

// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI headers
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/tbox/Database.h"

// PQS headers
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/GridManager.h"
#include "PQS/pqs/Solver.h"

// Class/type declarations
namespace SAMRAI { namespace mesh { class TagAndInitializeStrategy; } }
namespace PQS { namespace pqs { class DataInitStrategy; } }



// --- Implementation for PQS::pqs::Solver methods

namespace PQS {
namespace pqs {

Solver::Solver(
        boost::shared_ptr<tbox::Database> config_db,
        boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
        boost::shared_ptr<pqs::DataInitStrategy> data_init_strategy) {

    // Check parameters
    if (config_db == NULL) {
        // TODO
    }
    if (patch_hierarchy == NULL) {
        // TODO
    }
    if (data_init_strategy == NULL) {
        // TODO
    }

    // Set data members
    d_patch_hierarchy = patch_hierarchy;

    // Load configuration
    loadConfiguration(config_db);

    /*
     * Construct SAMRAI::mesh::GriddingAlgorithm
     */

    // Construct box generator and load balancer objects
    boost::shared_ptr<tbox::Database> box_generator_config_db(0);
    if (config_db->isDatabase("BergerRigoutsos")) {
        box_generator_config_db = config_db->getDatabase("BergerRigoutsos");
    }
    boost::shared_ptr<mesh::BergerRigoutsos> box_generator =
        boost::shared_ptr<mesh::BergerRigoutsos>(
            new mesh::BergerRigoutsos(d_patch_hierarchy->getDim(),
                box_generator_config_db));

    boost::shared_ptr<tbox::Database> load_balancer_config_db(0);
    if (config_db->isDatabase("ChopAndPackLoadBalancer")) {
        load_balancer_config_db =
            config_db->getDatabase("ChopAndPackLoadBalancer");
    }
    boost::shared_ptr<mesh::ChopAndPackLoadBalancer> load_balancer =
        boost::shared_ptr<mesh::ChopAndPackLoadBalancer> (
            new mesh::ChopAndPackLoadBalancer(d_patch_hierarchy->getDim(),
                "LoadBalancer", load_balancer_config_db));

    // Construct PQS::pqs::GridManager
    boost::shared_ptr<tbox::Database> grid_manager_config_db;
    if (config_db->isDatabase("GridManager")) {
        grid_manager_config_db = config_db->getDatabase("GridManager");
    } else {
        throw runtime_error(
            "'GridManager' section not found in configuration database");
    }
    boost::shared_ptr<pqs::GridManager> grid_manager =
        boost::shared_ptr<pqs::GridManager>(
            new pqs::GridManager(grid_manager_config_db, d_patch_hierarchy));


    // Construct SAMRAI::mesh::GriddingAlgorithm object
    d_gridding_alg = boost::shared_ptr<mesh::GriddingAlgorithm> (
        new mesh::GriddingAlgorithm(
            d_patch_hierarchy,
            "GriddingAlgorithm",
            grid_manager_config_db,
            boost::shared_ptr<mesh::TagAndInitializeStrategy>(grid_manager),
            box_generator,
            load_balancer));
}

void Solver::loadConfiguration(
        const boost::shared_ptr<tbox::Database>& config_db)
{
    // Check parameters
}

void Solver::printClassData(ostream& os) const
{
    os << endl
       << "===================================" << endl;
    os << "PQS::Solver" << endl;

    os << "Object Pointers" << endl;
    os << "---------------" << endl;
    os << "(Solver*) this = " << (Solver*) this << endl;
    os << "d_patch_hierarchy = " << d_patch_hierarchy.get() << endl;

    os << "===================================" << endl << endl;
}

} // PQS::pqs namespace
} // PQS namespace
