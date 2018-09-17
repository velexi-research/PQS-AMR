/*! \file test_lsm_kernels_2d.cc
 *
 * \brief
 * Unit tests for numerical kernels for level set method in 2D
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
#include "PQS/math/kernels/level_set_method_2d.h"

// PQS test
#include "fixture_kernels.h"
#include "kernels/kernels_2d.h"


// Namespaces
using namespace std;

// Class/type declarations


// --- Test cases

namespace lsmTests {

// Test case: compute area phi < 0 (no control volume)
TEST_F(lsmKernelTest, test_compute_area_phi_lt_zero)
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
    lsmKernelTest::setUpGrid();
    lsmKernelTest::setUpGridFunction();

    radius = 0.25;
    expected_area = M_PI * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_CIRCLE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    area = LSM_2D_AREA_PHI_LESS_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                          ib_lo, ib_hi,
                                          dx, &eps);
    EXPECT_NEAR((area - expected_area)/expected_area, 0, eps);

    lsmKernelTest::tearDownGridFunction();
    lsmKernelTest::tearDownGrid();

    // radius = 0.5, N = 100
    N = 100;
    lsmKernelTest::setUpGrid();
    lsmKernelTest::setUpGridFunction();

    radius = 0.5;
    expected_area = M_PI * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_CIRCLE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    area = LSM_2D_AREA_PHI_LESS_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                          ib_lo, ib_hi,
                                          dx, &eps);
    EXPECT_NEAR((area - expected_area)/expected_area, 0, eps);

    lsmKernelTest::tearDownGridFunction();
    lsmKernelTest::tearDownGrid();

    // radius = 0.75, N = 500
    N = 500;
    lsmKernelTest::setUpGrid();
    lsmKernelTest::setUpGridFunction();

    radius = 0.75;
    expected_area = M_PI * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_CIRCLE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    area = LSM_2D_AREA_PHI_LESS_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                          ib_lo, ib_hi,
                                          dx, &eps);
    EXPECT_NEAR((area - expected_area)/expected_area, 0, eps);

    lsmKernelTest::tearDownGridFunction();
    lsmKernelTest::tearDownGrid();
}

// Test case: compute area phi > 0 (no control volume)
TEST_F(lsmKernelTest, test_compute_area_phi_gt_zero)
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
    lsmKernelTest::setUpGrid();
    lsmKernelTest::setUpGridFunction();

    radius = 0.25;
    expected_area = 4 - M_PI * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_CIRCLE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    area = LSM_2D_AREA_PHI_GREATER_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                             ib_lo, ib_hi,
                                             dx, &eps);
    EXPECT_NEAR((area - expected_area)/expected_area, 0, eps);

    lsmKernelTest::tearDownGridFunction();
    lsmKernelTest::tearDownGrid();

    // radius = 0.5, N = 100
    N = 100;
    lsmKernelTest::setUpGrid();
    lsmKernelTest::setUpGridFunction();

    radius = 0.5;
    expected_area = 4 - M_PI * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_CIRCLE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    area = LSM_2D_AREA_PHI_GREATER_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                             ib_lo, ib_hi,
                                             dx, &eps);
    EXPECT_NEAR((area - expected_area)/expected_area, 0, eps);

    lsmKernelTest::tearDownGridFunction();
    lsmKernelTest::tearDownGrid();

    // radius = 0.75, N = 500
    N = 500;
    lsmKernelTest::setUpGrid();
    lsmKernelTest::setUpGridFunction();

    radius = 0.75;
    expected_area = 4 - M_PI * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_CIRCLE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    area = LSM_2D_AREA_PHI_GREATER_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                             ib_lo, ib_hi,
                                             dx, &eps);
    EXPECT_NEAR((area - expected_area)/expected_area, 0, eps);

    lsmKernelTest::tearDownGridFunction();
    lsmKernelTest::tearDownGrid();
}

} // lsmTests namespace
