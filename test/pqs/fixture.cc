/*! \file fixture.cc
 *
 * \brief
 * Implementation of fixture for unit tests for PQS::pqs classes.
 */

/*
 * ---------------------------------------------------------------------
 * COPYRIGHT/LICENSE. This file is part of the XYZ package. It is
 * subject to the license terms in the LICENSE file found in the
 * top-level directory of this distribution. No part of the XYZ
 * package, including this file, may be copied, modified, propagated,
 * or distributed except according to the terms contained in the
 * LICENSE file.
 * ---------------------------------------------------------------------
 */

// --- Headers, namespaces, and type declarations

// Standard library
#include <iosfwd>
#include <string>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep

// PQS test
#include "fixture.h"
#include "TestInterfaceInitModule.h"
#include "TestPoreInitModule.h"


// --- Fixtures

namespace pqsTests {

// Constructor (set up)
pqsTests::pqsTests() {

    // --- Initialize SAMRAI

    if (!tbox::SAMRAIManager::isInitialized()) {
        int argc = 0;
        char **argv = 0;
        tbox::SAMRAI_MPI::init(&argc, &argv);

        tbox::SAMRAIManager::initialize();
        tbox::SAMRAIManager::startup();
        const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    }

    // --- Construct configuration database

    config_db = boost::shared_ptr<tbox::MemoryDatabase>(
        new tbox::MemoryDatabase("Configuration Parameters"));

    // ------ PQS configuration

    // PQS database
    boost::shared_ptr<tbox::Database> pqs_config_db =
        config_db->putDatabase("PQS");

    pqs_config_db->putDouble("initial_curvature", 0.5);

    // ------ SAMRAI configuration

    // SAMRAI database
    boost::shared_ptr<tbox::Database> samrai_config_db =
        config_db->putDatabase("SAMRAI");

    // BoxGenerator database
    boost::shared_ptr<tbox::Database> box_generator_config_db =
        samrai_config_db->putDatabase("BoxGenerator");

    // LoadBalancer database
    boost::shared_ptr<tbox::Database> load_balancer_config_db =
        samrai_config_db->putDatabase("LoadBalancer");

    // GriddingAlgorithm database
    boost::shared_ptr<tbox::Database> gridding_algorithm_config_db =
        samrai_config_db->putDatabase("GriddingAlgorithm");

    // --------- Geometry configuration

    // Geometry database
    boost::shared_ptr<tbox::Database> geometry_config_db =
        samrai_config_db->putDatabase("Geometry");

    // Problem dimension
    const tbox::Dimension dim(3);

    // Physical bounds
    double x_lo[3] = {-1.0, -1.0, -1.0};
    double x_up[3] = {1.0, 1.0, 1.0};
    geometry_config_db->putDoubleArray("x_lo", x_lo, 3);
    geometry_config_db->putDoubleArray("x_up", x_up, 3);

    // Box
    int box_lower[3] = {0, 0, 0};
    int box_upper[3] = {49, 49, 49};
    tbox::DatabaseBox domain_boxes(dim, box_lower, box_upper);
    geometry_config_db->putDatabaseBoxArray("domain_boxes", &domain_boxes, 1);

    // --------- PatchHierarchy configuration

    // PatchHierarchy database
    boost::shared_ptr<tbox::Database> patch_hierarchy_config_db =
        samrai_config_db->putDatabase("PatchHierarchy");

    int max_levels = 3;
    patch_hierarchy_config_db->putInteger("max_levels", max_levels);
    patch_hierarchy_config_db->putDatabase("ratio_to_coarser");
    for (int ln=1; ln <= max_levels; ln++) {
        int ratio_to_coarser[3] = {2, 2, 2};
        std::string level_name = std::string("level_") + std::to_string(ln);
        patch_hierarchy_config_db->putIntegerArray(level_name,
                                                   ratio_to_coarser, 3);
    }

    // --- Initialize Geometry and PatchHierarchy

    // Geometry
    grid_geometry = boost::shared_ptr<geom::CartesianGridGeometry>(
        new geom::CartesianGridGeometry(
            dim, "CartesianGeometry", geometry_config_db));

    // PatchHierarchy
    patch_hierarchy = boost::shared_ptr<hier::PatchHierarchy>(
        new hier::PatchHierarchy("PatchHierarchy", grid_geometry,
                                 patch_hierarchy_config_db));

    // TestPoreInitModule (implements PoreInitStrategy)
    pore_init_strategy = boost::shared_ptr<pqs::PoreInitStrategy>(
        new TestPoreInitModule());

    // TestInterfaceInitModule (implements InterfaceInitStrategy)
    interface_init_strategy = boost::shared_ptr<pqs::InterfaceInitStrategy>(
        new TestInterfaceInitModule());
}

pqsTests::~pqsTests() {
    // Shutdown SAMRAI
    if (s_num_tests_remaining > 0) {
        s_num_tests_remaining--;
    } else {
        tbox::SAMRAIManager::shutdown();
        tbox::SAMRAIManager::finalize();
        tbox::SAMRAI_MPI::finalize();
    }
}

} // pqsTests namespace
