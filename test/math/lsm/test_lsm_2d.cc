/*! \file test_lsm_2d.cc
 *
 * \brief
 * Unit tests for level set method functionality in two space dimensions
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
#include <memory>

// Google Test
#include "gtest/gtest.h"
#include "gtest/gtest-message.h"
#include "gtest/gtest-test-part.h"

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/math/LSMToolbox.h"
#include "PQS/pqs/Solver.h"

// PQS test
#include "fixture.h"

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations


// --- Test cases

namespace lsmTests {

// Test case: PQS::math::computeVolume() in 2D
TEST_F(lsmTest, test_LSMToolbox_computeVolume_2d)
{
    // --- Preparations

    double radius;
    double area;
    double expected_area;

    const double eps = 1.5 * 2./100;

    // --- Exercise functionality and check results

    // ------ circle with radius = 0.5

    const int num_dimensions = 2;
    radius = 0.5;
    lsmTest::constructSolver(num_dimensions, radius);

    // interior of circle
    area = math::computeVolume(
            solver->getPatchHierarchy(),
            solver->getInterfacePatchDataId(),
            -1, // phi < 0
            solver->getControlVolumePatchDataId());

    expected_area = M_PI * radius * radius;

    EXPECT_NEAR((area - expected_area)/expected_area, 0, eps);

    // exterior of circle
    area = math::computeVolume(
            solver->getPatchHierarchy(),
            solver->getInterfacePatchDataId(),
            1, // phi > 0
            solver->getControlVolumePatchDataId());

    expected_area = 4 - M_PI * radius * radius;

    EXPECT_NEAR((area - expected_area)/expected_area, 0, eps);
}

} // lsmTests namespace
