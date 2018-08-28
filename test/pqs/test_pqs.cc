/*! \file test_pqs.cc
 *
 * \brief
 * Unit tests for PQS::pqs classes.
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
#include <memory>
#include <stddef.h>
#include <sstream>

// Google Test
#include "gtest/gtest.h"
#include "gtest/gtest-message.h"
#include "gtest/gtest-test-part.h"

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/Database.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/Solver.h"

// PQS test
#include "fixture.h"

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations


// --- Test cases

namespace pqsTests {

// Test case: validate structure of configuration database
TEST_F(pqsTest, test_config_db_structure)
{
    // Configure and initialize geometry and PatchHierarchy
    int num_dimensions = 3;
    pqsTest::initializeGeometryAndHierarchy(config_db,
                                            grid_geometry,
                                            patch_hierarchy,
                                            num_dimensions);

    // Sub-databases
    EXPECT_TRUE(config_db->isDatabase("PQS"));
    EXPECT_TRUE(config_db->isDatabase("SAMRAI"));

    // --- PQS database

    shared_ptr<tbox::Database> pqs_config_db =
        config_db->getDatabase("PQS");
    EXPECT_TRUE(pqs_config_db->isDouble("initial_curvature"));
    EXPECT_TRUE(pqs_config_db->isDouble("final_curvature"));
    EXPECT_TRUE(pqs_config_db->isDouble("curvature_increment"));

    // ------ Algorithms database

    EXPECT_TRUE(pqs_config_db->isDatabase("Algorithms"));
    shared_ptr<tbox::Database> algorithms_config_db =
        pqs_config_db->getDatabase("Algorithms");

    // PrescribedCurvatureModel database
    EXPECT_TRUE(algorithms_config_db->isDatabase("PrescribedCurvatureModel"));

    // SlightlyCompressibleModel database
    EXPECT_TRUE(algorithms_config_db->isDatabase("SlightlyCompressibleModel"));

    // --- SAMRAI database

    shared_ptr<tbox::Database> samrai_config_db =
        config_db->getDatabase("SAMRAI");
    EXPECT_TRUE(samrai_config_db->isDatabase("Geometry"));
    EXPECT_TRUE(samrai_config_db->isDatabase("PatchHierarchy"));
    EXPECT_TRUE(samrai_config_db->isDatabase("BoxGenerator"));
    EXPECT_TRUE(samrai_config_db->isDatabase("LoadBalancer"));
    EXPECT_TRUE(samrai_config_db->isDatabase("GriddingAlgorithm"));

    // Geometry database
    shared_ptr<tbox::Database> geometry_config_db =
        samrai_config_db->getDatabase("Geometry");
    EXPECT_TRUE(geometry_config_db->keyExists("dim"));
    EXPECT_TRUE(geometry_config_db->keyExists("x_lo"));
    EXPECT_TRUE(geometry_config_db->keyExists("x_up"));
    EXPECT_TRUE(geometry_config_db->keyExists("domain_boxes"));
}

// Test case: test constructor for PQS::pqs::Solver class with PatchHierarchy
TEST_F(pqsTest, test_Solver_Solver_with_patch_hierarchy)
{
    // --- Preparations

    // Configure and initialize geometry and PatchHierarchy
    int num_dimensions = 3;
    pqsTest::initializeGeometryAndHierarchy(config_db,
                                            grid_geometry,
                                            patch_hierarchy,
                                            num_dimensions);

    // Configuration databases
    shared_ptr<tbox::Database> pqs_config_db =
        config_db->getDatabase("PQS");

    shared_ptr<tbox::Database> samrai_config_db =
        config_db->getDatabase("SAMRAI");

    // --- Exercise functionality

    // Construct PQS::pqs::Solver object
    solver = new pqs::Solver(config_db,
                             pore_init_strategy,
                             interface_init_strategy,
                             patch_hierarchy);

    // --- Check results

    // Status code
    EXPECT_NE(solver, (pqs::Solver*) NULL);

    // Solver parameters
    EXPECT_EQ(solver->getInitialCurvature(),
              pqs_config_db->getDouble("initial_curvature"));
    EXPECT_EQ(solver->getFinalCurvature(),
              pqs_config_db->getDouble("final_curvature"));
    EXPECT_EQ(solver->getStepCount(), 0);

    // Solver state
    EXPECT_GE(solver->getPoreSpacePatchDataId(), 0);
    EXPECT_GE(solver->getInterfacePatchDataId(), 0);

    // PatchHierarchy configuration
    shared_ptr<hier::PatchHierarchy> patch_hierarchy =
        solver->getPatchHierarchy();
    EXPECT_TRUE(patch_hierarchy);  // check that pointer is not NULL

    int expected_num_patch_levels = samrai_config_db->
        getDatabase("PatchHierarchy")->getInteger("max_levels");
    EXPECT_EQ(patch_hierarchy->getNumberOfLevels(), expected_num_patch_levels);
}

// Test case: test constructor for PQS::pqs::Solver class without PatchHierarchy
TEST_F(pqsTest, test_Solver_Solver_without_patch_hierarchy)
{
    // --- Preparations

    // Configure and initialize geometry and PatchHierarchy
    int num_dimensions = 3;
    pqsTest::initializeGeometryAndHierarchy(config_db,
                                            grid_geometry,
                                            patch_hierarchy,
                                            num_dimensions);

    // Configuration databases
    shared_ptr<tbox::Database> pqs_config_db =
        config_db->getDatabase("PQS");

    shared_ptr<tbox::Database> samrai_config_db =
        config_db->getDatabase("SAMRAI");

    // --- Exercise functionality

    // Construct PQS::pqs::Solver object
    solver = new pqs::Solver(config_db, pore_init_strategy,
                             interface_init_strategy);

    // --- Check results

    // Status code
    EXPECT_NE(solver, (pqs::Solver*) NULL);

    // Solver parameters
    EXPECT_EQ(solver->getInitialCurvature(),
              pqs_config_db->getDouble("initial_curvature"));
    EXPECT_EQ(solver->getFinalCurvature(),
              pqs_config_db->getDouble("final_curvature"));
    EXPECT_EQ(solver->getStepCount(), 0);

    // Solver state
    EXPECT_GE(solver->getPoreSpacePatchDataId(), 0);
    EXPECT_GE(solver->getInterfacePatchDataId(), 0);

    // PatchHierarchy configuration
    shared_ptr<hier::PatchHierarchy> patch_hierarchy =
        solver->getPatchHierarchy();
    EXPECT_TRUE(patch_hierarchy);  // check that pointer is not NULL

    int expected_num_patch_levels = samrai_config_db->
        getDatabase("PatchHierarchy")->getInteger("max_levels");
    EXPECT_EQ(patch_hierarchy->getNumberOfLevels(), expected_num_patch_levels);
}

} // pqsTests namespace
