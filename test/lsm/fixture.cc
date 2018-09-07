/*! \file fixture.cc
 *
 * \brief
 * Implementation of fixture for unit tests for PQS::lsm classes.
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
#include "PQS/utilities/error.h"
#include "PQS/pqs/Solver.h"

// PQS test
#include "fixture.h"
#include "TestInterfaceInitModule.h"
#include "TestPoreInitModule.h"


// --- Fixtures

namespace lsmTests {

// Static data members
int lsmTest::s_num_tests = 1;
int lsmTest::s_num_tests_remaining = lsmTest::s_num_tests;

// Methods
lsmTest::lsmTest() {

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

    config_db = shared_ptr<tbox::Database>(
        new tbox::MemoryDatabase("Configuration Parameters"));

    // ------ PQS configuration

    // PQS database
    shared_ptr<tbox::Database> pqs_config_db =
        config_db->putDatabase("PQS");

    // Physical parameters
    pqs_config_db->putDouble("initial_curvature", 0.5);
    pqs_config_db->putDouble("final_curvature", 1.0);
    pqs_config_db->putDouble("curvature_step", 0.1);

    // Level set method parameters
    pqs_config_db->putInteger("lsm_t_max", 1.0);
    pqs_config_db->putInteger("lsm_max_iterations", 5);
    pqs_config_db->putDouble("lsm_min_delta_phi", 0.1);
    pqs_config_db->putDouble("lsm_min_delta_saturation", 0.1);

    // Numerical method parameters
    pqs_config_db->putString("lsm_spatial_derivative_type", "WENO5");
    pqs_config_db->putInteger("time_integration_order", 1);

    // Algorithms database
    shared_ptr<tbox::Database> algorithms_config_db =
        pqs_config_db->putDatabase("Algorithms");

    // Prescribed Curvature Model database
    shared_ptr<tbox::Database> pcm_config_db =
        algorithms_config_db->putDatabase("PrescribedCurvatureModel");
    pcm_config_db->putDouble("pressure", 1.0);
    pcm_config_db->putDouble("surface_tension", 0.1);

    // Slightly Compressible Model database
    shared_ptr<tbox::Database> scm_config_db =
        algorithms_config_db->putDatabase("SlightlyCompressibleModel");
    scm_config_db->putDouble("pressure", 1.0);
    scm_config_db->putDouble("bulk_modulus", 1.0);
    scm_config_db->putDouble("target_volume", 0.5);
    scm_config_db->putDouble("surface_tension", 0.1);

    // ------ SAMRAI configuration

    // SAMRAI database
    shared_ptr<tbox::Database> samrai_config_db =
        config_db->putDatabase("SAMRAI");

    // BoxGenerator database
    shared_ptr<tbox::Database> box_generator_config_db =
        samrai_config_db->putDatabase("BoxGenerator");

    // LoadBalancer database
    shared_ptr<tbox::Database> load_balancer_config_db =
        samrai_config_db->putDatabase("LoadBalancer");

    // GriddingAlgorithm database
    shared_ptr<tbox::Database> gridding_algorithm_config_db =
        samrai_config_db->putDatabase("GriddingAlgorithm");

    // --- Initialize PQS objects

    // Solver
    solver = shared_ptr<pqs::Solver>(NULL);

    // TestInterfaceInitModule (implements InterfaceInitStrategy)
    interface_init_strategy = shared_ptr<pqs::InterfaceInitStrategy>(NULL);

    // TestPoreInitModule (implements PoreInitStrategy)
    pore_init_strategy = shared_ptr<pqs::PoreInitStrategy>(
        new TestPoreInitModule());

} // lsmTest::lsmTest()

lsmTest::~lsmTest() {
    // Clean up SAMRAI objects
    patch_hierarchy.reset();
    grid_geometry.reset();

    // Clean up PQS objects
    solver.reset(); 
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
} // lsmTest::~lsmTest()

void lsmTest::constructSolver(const int num_dimensions, const double radius)
{
    // --- Check arguments

    if (config_db == NULL) {
        PQS_ERROR(this, "constructSolver",
                  "'config_db' must not be NULL");
    }
    if ((num_dimensions != 2) && (num_dimensions != 3)) {
        PQS_ERROR(this, "constructSolver",
                  std::string("'num_dimensions' (= ") +
                  std::to_string(num_dimensions) +
                  std::string(") must be equal to 2 or 3"));
    }
    if (radius <= 0) {
        PQS_ERROR(this, "constructSolver",
                  std::string("'radius' (= ") +
                  std::to_string(radius) +
                  std::string(") must be greater than zero"));
    }

    // --- Preparations

    shared_ptr<tbox::Database> samrai_config_db =
        config_db->getDatabase("SAMRAI");

    // Geometry database
    shared_ptr<tbox::Database> geometry_config_db =
        samrai_config_db->putDatabase("Geometry");

    // PatchHierarchy database
    shared_ptr<tbox::Database> patch_hierarchy_config_db =
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
        int box_upper[2] = {99, 99};
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
        int box_upper[3] = {99, 99, 99};
        tbox::DatabaseBox domain_boxes(dim, box_lower, box_upper);
        geometry_config_db->putDatabaseBoxArray("domain_boxes",
                                                &domain_boxes, 1);
    }

    // --- PatchHierarchy configuration

    int max_levels = 1;  // TODO: add multi-level tests
    patch_hierarchy_config_db->putInteger("max_levels", max_levels);
    shared_ptr<tbox::Database> ratio_to_coarser_db =
            patch_hierarchy_config_db->putDatabase("ratio_to_coarser");
    for (int ln=1; ln <= max_levels; ln++) {
        std::string level_name = std::string("level_") + std::to_string(ln);

        if (num_dimensions == 2) {
            int ratio_to_coarser[2] = {2, 2};
            ratio_to_coarser_db->putIntegerArray(level_name,
                                                       ratio_to_coarser, 2);
        } else {
            int ratio_to_coarser[3] = {2, 2, 2};
            ratio_to_coarser_db->putIntegerArray(level_name,
                                                       ratio_to_coarser, 3);
        }
    }

    // --- Initialize Geometry and PatchHierarchy

    // Geometry
    grid_geometry = shared_ptr<geom::CartesianGridGeometry>(
        new geom::CartesianGridGeometry(dim, "TestCartesianGeometry",
                                        geometry_config_db));

    // PatchHierarchy
    patch_hierarchy = shared_ptr<hier::PatchHierarchy>(
        new hier::PatchHierarchy("TestPatchHierarchy", grid_geometry,
                                 patch_hierarchy_config_db));

    // --- Construct TestInterfaceInitModule object

    interface_init_strategy = shared_ptr<pqs::InterfaceInitStrategy>(
            new TestInterfaceInitModule(num_dimensions, radius));

    // --- Construct PQS::pqs::Solver object

    solver = shared_ptr<pqs::Solver>(
            new pqs::Solver(config_db,
                            pore_init_strategy,
                            interface_init_strategy));

} // lsmTest::constructSolver()

} // lsmTests namespace
