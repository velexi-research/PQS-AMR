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
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/tbox/Database.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/Solver.h"
#include "PQS/pqs/TagAndInitModule.h"

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

    // Construct PQS::pqs::TagAndInitModule
    boost::shared_ptr<tbox::Database> tag_and_init_module_config_db;
    if (config_db->isDatabase("TagAndInitModule")) {
        tag_and_init_module_config_db =
            config_db->getDatabase("TagAndInitModule");
    } else {
        throw runtime_error(
            "'TagAndInitModule' section not found in configuration database");
    }
    boost::shared_ptr<pqs::TagAndInitModule> d_tag_and_init_module =
        boost::shared_ptr<pqs::TagAndInitModule>(
            new pqs::TagAndInitModule(tag_and_init_module_config_db,
                                      d_patch_hierarchy));


    // Construct SAMRAI::mesh::GriddingAlgorithm object
    d_gridding_alg = boost::shared_ptr<mesh::GriddingAlgorithm> (
        new mesh::GriddingAlgorithm(
            d_patch_hierarchy,
            "GriddingAlgorithm",
            tag_and_init_module_config_db, // TODO
            boost::shared_ptr<mesh::TagAndInitializeStrategy>(
                d_tag_and_init_module),
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
}

} // PQS::pqs namespace
} // PQS namespace
