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
#include <stddef.h>
#include <string>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

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
#include "PQS/utilities.h"
#include "PQS/pqs/Solver.h"

// PQS test
#include "fixture.h"
#include "TestInterfaceInitModule.h"
#include "TestPoreInitModule.h"


// --- Fixtures

namespace pqsTests {

// Static data members
int pqsTest::s_num_tests = 7;
int pqsTest::s_num_tests_remaining = pqsTest::s_num_tests;

// Methods
pqsTest::pqsTest() {

    // --- Initialize SAMRAI

    if (!tbox::SAMRAIManager::isInitialized()) {
        int argc = 0;
        char **argv = 0;
        tbox::SAMRAI_MPI::init(&argc, &argv);

        tbox::SAMRAIManager::initialize();
        const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    }
    tbox::SAMRAIManager::startup();

    // --- Construct configuration database

    config_db = boost::shared_ptr<tbox::Database>(
        new tbox::MemoryDatabase("Configuration Parameters"));

    // ------ PQS configuration

    // PQS database
    boost::shared_ptr<tbox::Database> pqs_config_db =
        config_db->putDatabase("PQS");

    // Physical parameters
    pqs_config_db->putDouble("initial_curvature", 0.5);

    // Level set method parameters
    pqs_config_db->putInteger("lsm_t_max", 1.0);
    pqs_config_db->putInteger("lsm_max_iterations", 5);
    pqs_config_db->putDouble("lsm_min_delta_phi", 0.1);
    pqs_config_db->putDouble("lsm_min_delta_saturation", 0.1);

    // Numerical method parameters
    pqs_config_db->putString("lsm_spatial_derivative_type", "WENO5");
    pqs_config_db->putInteger("time_integration_order", 1);

    // Algorithms database
    boost::shared_ptr<tbox::Database> algorithms_config_db =
        pqs_config_db->putDatabase("Algorithms");

    // Slightly Compressible Model database
    boost::shared_ptr<tbox::Database> scm_config_db =
        algorithms_config_db->putDatabase("SlightlyCompressibleModel");
    scm_config_db->putDouble("reference_pressure", 1.0);
    scm_config_db->putDouble("bulk_modulus", 1.0);
    scm_config_db->putDouble("target_volume", 0.5);
    scm_config_db->putDouble("surface_tension", 0.1);

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

    // --- Initialize PQS objects

    // Solver
    solver = NULL;

    // TestPoreInitModule (implements PoreInitStrategy)
    pore_init_strategy = boost::shared_ptr<pqs::PoreInitStrategy>(
        new TestPoreInitModule());

    // TestInterfaceInitModule (implements InterfaceInitStrategy)
    interface_init_strategy = boost::shared_ptr<pqs::InterfaceInitStrategy>(
        new TestInterfaceInitModule());
} // pqsTest::pqsTest()

pqsTest::~pqsTest() {
    // Clean up SAMRAI objects
    patch_hierarchy.reset();
    grid_geometry.reset();

    // Clean up PQS objects
    if (solver) {
        delete solver;
        solver = NULL;
    }
    pore_init_strategy.reset();
    interface_init_strategy.reset();
    config_db.reset();

    // Shutdown SAMRAI
    tbox::SAMRAIManager::shutdown();
    if (s_num_tests_remaining > 0) {
        s_num_tests_remaining--;
    } else {
        tbox::SAMRAIManager::finalize();
        tbox::SAMRAI_MPI::finalize();
    }
} // pqsTest::~pqsTest()

void pqsTest::initializeGeometryAndHierarchy(
        boost::shared_ptr<tbox::Database> config_db,
        boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
        boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int num_dimensions)
{
    // --- Check arguments

    if (config_db == NULL) {
        PQS_ERROR_STATIC("pqsTest", "initializeGeometryAndHierarchy",
                  "'config_db' must not be NULL");
    }
    if ((num_dimensions != 2) && (num_dimensions != 3)) {
        PQS_ERROR_STATIC("pqsTest", "initializeGeometryAndHierarchy",
                  std::string("'num_dimensions' (= ") +
                  std::to_string(num_dimensions) +
                  std::string(") must be equal to 2 or 3"));
    }

    // --- Preparations

    boost::shared_ptr<tbox::Database> samrai_config_db =
        config_db->getDatabase("SAMRAI");

    // Geometry database
    boost::shared_ptr<tbox::Database> geometry_config_db =
        samrai_config_db->putDatabase("Geometry");

    // PatchHierarchy database
    boost::shared_ptr<tbox::Database> patch_hierarchy_config_db =
        samrai_config_db->putDatabase("PatchHierarchy");

    // --- Geometry configuration

    // Problem dimension
    geometry_config_db->putInteger("dim", num_dimensions);
    const tbox::Dimension dim(geometry_config_db->getInteger("dim"));

    if (num_dimensions == 2) {
        // Physical bounds
        double x_lo[2] = {-1.0, -1.0};
        double x_up[2] = {1.0, 1.0};
        geometry_config_db->putDoubleArray("x_lo", x_lo, 2);
        geometry_config_db->putDoubleArray("x_up", x_up, 2);

        // Box
        int box_lower[2] = {0, 0};
        int box_upper[2] = {19, 19};
        const tbox::Dimension dim(geometry_config_db->getInteger("dim"));
        tbox::DatabaseBox domain_boxes(dim, box_lower, box_upper);
        geometry_config_db->putDatabaseBoxArray("domain_boxes",
                                                &domain_boxes, 1);
    } else {
        // Physical bounds
        double x_lo[3] = {-1.0, -1.0, -1.0};
        double x_up[3] = {1.0, 1.0, 1.0};
        geometry_config_db->putDoubleArray("x_lo", x_lo, 3);
        geometry_config_db->putDoubleArray("x_up", x_up, 3);

        // Box
        int box_lower[3] = {0, 0, 0};
        int box_upper[3] = {19, 19, 19};
        tbox::DatabaseBox domain_boxes(dim, box_lower, box_upper);
        geometry_config_db->putDatabaseBoxArray("domain_boxes",
                                                &domain_boxes, 1);
    }

    // --- PatchHierarchy configuration

    int max_levels = 3;
    patch_hierarchy_config_db->putInteger("max_levels", max_levels);
    patch_hierarchy_config_db->putDatabase("ratio_to_coarser");
    for (int ln=1; ln <= max_levels; ln++) {
        std::string level_name = std::string("level_") + std::to_string(ln);

        if (num_dimensions == 2) {
            int ratio_to_coarser[2] = {2, 2};
            patch_hierarchy_config_db->putIntegerArray(level_name,
                                                       ratio_to_coarser, 2);
        } else {
            int ratio_to_coarser[3] = {2, 2, 2};
            patch_hierarchy_config_db->putIntegerArray(level_name,
                                                       ratio_to_coarser, 3);
        }
    }

    // --- Initialize Geometry and PatchHierarchy

    // Geometry
    grid_geometry = boost::shared_ptr<geom::CartesianGridGeometry>(
        new geom::CartesianGridGeometry(dim, "TestCartesianGeometry",
                                        geometry_config_db));

    // PatchHierarchy
    patch_hierarchy = boost::shared_ptr<hier::PatchHierarchy>(
        new hier::PatchHierarchy("TestPatchHierarchy", grid_geometry,
                                 patch_hierarchy_config_db));

} // pqsTest::configureGeometryAndHierarchy()

} // pqsTests namespace
