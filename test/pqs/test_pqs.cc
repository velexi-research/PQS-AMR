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
#include "SAMRAI/tbox/MemoryDatabase.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/DataInitStrategy.h"
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

// Test case: construct PQS::pqs::Solver object
TEST_F(pqsTests, Solver_Constructor)
{
    // --- Preparations

    // Construct configuration parameter database
    boost::shared_ptr<tbox::MemoryDatabase> config_db =
        boost::shared_ptr<tbox::MemoryDatabase>(
            new tbox::MemoryDatabase("Solver"));

    config_db->putDatabase("GridManager");

    // Construct DataInitStrategy object
    boost::shared_ptr<pqs::DataInitStrategy> data_init_strategy =
        boost::shared_ptr<pqs::DataInitStrategy>(
            new pqs::DataInitStrategy());

    // --- Exercise functionality

    // Construct PQS::pqs::Solver object
    pqs::Solver *solver =
        new pqs::Solver(config_db, patch_hierarchy, data_init_strategy);

    // --- Check results

    // Status code
    EXPECT_NE(solver, (pqs::Solver*) NULL);
}

} // pqsTests namespace

