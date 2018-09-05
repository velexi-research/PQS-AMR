/*! \file test_lsm_3d.cc
 *
 * \brief
 * Unit tests for PQS::lsm classes.
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

// Google Test
#include "gtest/gtest.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/lsm/kernels_3d.h"

// PQS test
#include "fixture.h"

// Namespaces
using namespace std;

// Class/type declarations


// --- Test cases

namespace lsmTests {

// Test case: compute volume (no control volume)
TEST_F(lsmTest, test_compute_volume)
{
    // --- Preparations

    // Configure and initialize geometry and PatchHierarchy
    num_dimensions = 3;
    N = 10;
    lsmTest::setUpGrid();
    lsmTest::setUpGridFunction();

    // --- Exercise functionality

    // TODO

    // --- Check results

    // TODO
}

} // lsmTests namespace
