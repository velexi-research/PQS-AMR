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
#include <memory>

// Google Test
#include "gtest/gtest.h"
#include "gtest/gtest-message.h"
#include "gtest/gtest-test-part.h"

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/lsm/Toolbox.h"
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

// Test case: PQS::lsm::Toolbox::computeVolume() in 3D
TEST_F(lsmTest, test_Toolbox_computeVolume_3d)
{
    // --- Preparations

    double radius;
    double volume;
    double expected_volume;

    const double eps = 1.5 * 2./100;

    // --- Exercise functionality and check results

    // ------ sphere with radius = 0.5

    const int num_dimensions = 3;
    radius = 0.5;
    lsmTest::constructSolver(num_dimensions, radius);

    // interior of sphere
    volume = lsm::Toolbox::computeVolume(
            solver->getPatchHierarchy(),
            solver->getInterfacePatchDataId(),
            solver->getControlVolumePatchDataId(),
            -1); // phi < 0

    expected_volume = M_PI * 4 / 3 * radius * radius * radius;

    EXPECT_NEAR((volume - expected_volume)/expected_volume, 0, eps);

    // exterior of sphere
    volume = lsm::Toolbox::computeVolume(
            solver->getPatchHierarchy(),
            solver->getInterfacePatchDataId(),
            solver->getControlVolumePatchDataId(),
            1); // phi < 0

    expected_volume = 8 - M_PI * 4 / 3 * radius * radius * radius;

    EXPECT_NEAR((volume - expected_volume)/expected_volume, 0, eps);
}

} // lsmTests namespace
