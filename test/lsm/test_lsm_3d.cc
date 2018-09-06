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
#include <cmath>

// Google Test
#include "gtest/gtest.h"
#include "gtest/gtest-message.h"
#include "gtest/gtest-test-part.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/lsm/kernels_3d.h"

// PQS test
#include "fixture.h"
#include "kernels/kernels_3d.h"


// Namespaces
using namespace std;

// Class/type declarations


// --- Test cases

namespace lsmTests {

// Test case: compute volume phi < 0 (no control volume)
TEST_F(lsmTest, test_compute_volume_phi_lt_zero)
{
    // --- Preparations

    num_dimensions = 3;

    // --- Exercise functionality and check results

    // ------ circles

    double radius;
    double volume, expected_volume;
    double eps;

    // radius = 0.25, N = 200
    N = 200;
    lsmTest::setUpGrid();
    lsmTest::setUpGridFunction();

    radius = 0.25;
    expected_volume = 4 * M_PI / 3 * radius * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_SPHERE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    volume = LSM_3D_VOLUME_PHI_LESS_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                              ib_lo, ib_hi,
                                              dx, &eps);
    EXPECT_NEAR((volume - expected_volume)/expected_volume, 0, eps);

    lsmTest::tearDownGridFunction();
    lsmTest::tearDownGrid();

    // radius = 0.5, N = 100
    N = 100;
    lsmTest::setUpGrid();
    lsmTest::setUpGridFunction();

    radius = 0.5;
    expected_volume = 4 * M_PI / 3 * radius * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_SPHERE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    volume = LSM_3D_VOLUME_PHI_LESS_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                              ib_lo, ib_hi,
                                              dx, &eps);
    EXPECT_NEAR((volume - expected_volume)/expected_volume, 0, eps);

    lsmTest::tearDownGridFunction();
    lsmTest::tearDownGrid();

    // radius = 0.75, N = 250
    N = 250;
    lsmTest::setUpGrid();
    lsmTest::setUpGridFunction();

    radius = 0.75;
    expected_volume = 4 * M_PI / 3 * radius * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_SPHERE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    volume = LSM_3D_VOLUME_PHI_LESS_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                              ib_lo, ib_hi,
                                              dx, &eps);
    EXPECT_NEAR((volume - expected_volume)/expected_volume, 0, eps);

    lsmTest::tearDownGridFunction();
    lsmTest::tearDownGrid();

    // --- Clean up
    //
}

// Test case: compute volume phi > 0 (no control volume)
TEST_F(lsmTest, test_compute_volume_phi_gt_zero)
{
    // --- Preparations

    num_dimensions = 3;

    // --- Exercise functionality and check results

    // ------ circles

    double radius;
    double volume, expected_volume;
    double eps;

    // radius = 0.25, N = 200
    N = 200;
    lsmTest::setUpGrid();
    lsmTest::setUpGridFunction();

    radius = 0.25;
    expected_volume = 8 - 4 * M_PI / 3 * radius * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_SPHERE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    volume = LSM_3D_VOLUME_PHI_GREATER_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                                 ib_lo, ib_hi,
                                                 dx, &eps);
    EXPECT_NEAR((volume - expected_volume)/expected_volume, 0, eps);

    lsmTest::tearDownGridFunction();
    lsmTest::tearDownGrid();

    // radius = 0.5, N = 100
    N = 100;
    lsmTest::setUpGrid();
    lsmTest::setUpGridFunction();

    radius = 0.5;
    expected_volume = 8 - 4 * M_PI / 3 * radius * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_SPHERE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    volume = LSM_3D_VOLUME_PHI_GREATER_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                                 ib_lo, ib_hi,
                                                 dx, &eps);
    EXPECT_NEAR((volume - expected_volume)/expected_volume, 0, eps);

    lsmTest::tearDownGridFunction();
    lsmTest::tearDownGrid();

    // radius = 0.75, N = 250
    N = 250;
    lsmTest::setUpGrid();
    lsmTest::setUpGridFunction();

    radius = 0.75;
    expected_volume = 8 - 4 * M_PI / 3 * radius * radius * radius;
    eps = 1.5 * dx[0];
    SET_PHI_SPHERE(phi, phi_gb_lo, phi_gb_hi,
                   ib_lo, ib_hi,
                   x_lo, dx, &radius);
    volume = LSM_3D_VOLUME_PHI_GREATER_THAN_ZERO(phi, phi_gb_lo, phi_gb_hi,
                                                 ib_lo, ib_hi,
                                                 dx, &eps);
    EXPECT_NEAR((volume - expected_volume)/expected_volume, 0, eps);

    lsmTest::tearDownGridFunction();
    lsmTest::tearDownGrid();

    // --- Clean up
    //
}

} // lsmTests namespace
