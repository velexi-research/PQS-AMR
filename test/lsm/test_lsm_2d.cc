/*! \file test_lsm_2d.cc
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
#include <cmath>

// Google Test
#include "gtest/gtest.h"
#include "gtest/gtest-message.h"
#include "gtest/gtest-test-part.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/lsm/kernels_2d.h"

// PQS test
#include "fixture.h"
#include "kernels/kernels_2d.h"


// Namespaces
using namespace std;

// Class/type declarations


// --- Test cases

namespace lsmTests {

// Test case: compute area (no control volume)
TEST_F(lsmTest, test_compute_area)
{
    // --- Preparations

    num_dimensions = 2;

    // --- Exercise functionality and check results

    // ------ circles

    double radius;
    double area, expected_area;
    double eps;

    // radius = 0.25, N = 200
    N = 200;
    lsmTest::setUpGrid();
    lsmTest::setUpGridFunction();

    radius = 0.25;
    expected_area = M_PI * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_CIRCLE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    area = PQS_LSM_2D_AREA_PHI_LESS_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                              ib_lo, ib_hi,
                                              dx, &eps);
    EXPECT_NEAR((area - expected_area)/expected_area, 0, eps);

    lsmTest::tearDownGridFunction();
    lsmTest::tearDownGrid();

    // radius = 0.5, N = 100
    N = 100;
    lsmTest::setUpGrid();
    lsmTest::setUpGridFunction();

    radius = 0.5;
    expected_area = M_PI * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_CIRCLE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    area = PQS_LSM_2D_AREA_PHI_LESS_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                              ib_lo, ib_hi,
                                              dx, &eps);
    EXPECT_NEAR((area - expected_area)/expected_area, 0, eps);

    lsmTest::tearDownGridFunction();
    lsmTest::tearDownGrid();

    // radius = 0.75, N = 500
    N = 500;
    lsmTest::setUpGrid();
    lsmTest::setUpGridFunction();

    radius = 0.75;
    expected_area = M_PI * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_CIRCLE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    area = PQS_LSM_2D_AREA_PHI_LESS_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                              ib_lo, ib_hi,
                                              dx, &eps);
    EXPECT_NEAR((area - expected_area)/expected_area, 0, eps);

    lsmTest::tearDownGridFunction();
    lsmTest::tearDownGrid();

    // --- Clean up
    //
}

} // lsmTests namespace
