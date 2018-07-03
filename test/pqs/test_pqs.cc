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
#include <sstream>
#include <stddef.h>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// Google Test
#include "gtest/gtest.h"
#include "gtest/gtest-message.h"
#include "gtest/gtest-test-part.h"

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MemoryDatabase.h"

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
TEST_F(pqsTests, config_db_structure)
{
    // --- Validate structure of configuration database

    // Sub-databases
    EXPECT_TRUE(config_db->isDatabase("PQS"));
    EXPECT_TRUE(config_db->isDatabase("Geometry"));
    EXPECT_TRUE(config_db->isDatabase("SAMRAI"));

    // PQS database
    boost::shared_ptr<tbox::Database> pqs_config_db =
        config_db->getDatabase("PQS");
    EXPECT_TRUE(pqs_config_db->isDouble("initial_curvature"));

    // Geometry database
    boost::shared_ptr<tbox::Database> geometry_config_db =
        config_db->getDatabase("Geometry");
    EXPECT_TRUE(geometry_config_db->keyExists("x_lo"));
    EXPECT_TRUE(geometry_config_db->keyExists("x_up"));
    EXPECT_TRUE(geometry_config_db->keyExists("domain_boxes"));

    // SAMRAI database
    boost::shared_ptr<tbox::Database> samrai_config_db =
        config_db->getDatabase("SAMRAI");
    EXPECT_TRUE(samrai_config_db->isDatabase("BoxGenerator"));
    EXPECT_TRUE(samrai_config_db->isDatabase("LoadBalancer"));
    EXPECT_TRUE(samrai_config_db->isDatabase("GriddingAlgorithm"));
}

// Test case: construct PQS::pqs::Solver object
TEST_F(pqsTests, Solver_Solver)
{
    // --- Preparations


    // --- Exercise functionality

    // Construct PQS::pqs::Solver object
    pqs::Solver *solver =
        new pqs::Solver(config_db, patch_hierarchy, data_init_strategy);

    // --- Check results

    // Status code
    EXPECT_NE(solver, (pqs::Solver*) NULL);

    // Solver parameters
    boost::shared_ptr<tbox::Database> pqs_config_db =
        config_db->getDatabase("PQS");
    EXPECT_EQ(solver->getCurvature(),
              pqs_config_db->getDouble("initial_curvature"));
    EXPECT_EQ(solver->getStep(), 0);

    // PatchHierarchy configuration
    // TODO
}

} // pqsTests namespace

