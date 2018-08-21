/*! \file test_pqs_3d.cc
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

// Google Test
#include "gtest/gtest.h"

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
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

// Test case: PQS::pqs::Solver::equilibrateInterface()
TEST_F(pqsTest, test_Solver_equilibrateInterface_3d)
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

    // Construct PQS::pqs::Solver object
    solver = new pqs::Solver(config_db, pore_init_strategy,
                             interface_init_strategy);

    // mean curvature
    double curvature = 1.0;

    // --- Exercise functionality

    solver->equilibrateInterface(curvature);

    // --- Check results

    // TODO
}

// Test case: PQS::pqs::Solver::advanceInterface()
TEST_F(pqsTest, test_Solver_advanceInterface_3d)
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

    // TODO
}

} // pqsTests namespace
