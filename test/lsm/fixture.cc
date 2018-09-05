/*! \file fixture.cc
 *
 * \brief
 * Implementation of fixture for unit tests for PQS::lsm package
 */

/*
 * ---------------------------------------------------------------------
 * COPYRIGHT/LICENSE. This file is part of the XYZ package. It is
 * subject to the license terms in the LICENSE file found in the
 * top-level directory of this distribution. No part of the XYZ
 * package, including this file, may be copied, modified, propagated,
 * or distributed except according to the terms contained in the
 * LICENSE file.
 * ---------------------------------------------------------------------
 */

// --- Headers, namespaces, and type declarations

// Standard library
#include <memory>
#include <stddef.h>
#include <string>

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/utilities/error.h"

// PQS test
#include "fixture.h"


// --- Fixtures

namespace lsmTests {

// Static data members
int lsmTest::s_num_tests = 0;
int lsmTest::s_num_tests_remaining = lsmTest::s_num_tests;

// Methods
lsmTest::lsmTest() {
    // Initialize num_dimensions and N to invalid values
    num_dimensions = 0;
    N = 0;

    // Initalize pointers to NULL
    x_lo = NULL;
    x_hi = NULL;

    // Initialze grid parameters to NULL
    dx = NULL;
    phi_gb_lo = NULL;
    phi_gb_hi = NULL;
    ib_lo = NULL;
    ib_hi = NULL;

    // Initialize phi to NULL
    phi = NULL;

    // Initialize grid_setup_complete to false
    grid_setup_complete = false;

} // lsmTest::lsmTest()

lsmTest::~lsmTest() {
    // Tear down grid
    tearDownGrid();

    // Tear down grid functdion
    tearDownGridFunction();

} // lsmTest::~lsmTest()

void lsmTest::setUpGrid()
{
    // --- Check prerequisites

    if ((num_dimensions != 2) && (num_dimensions != 3)) {
        PQS_ERROR_STATIC("lsmTest", "setUpGrid",
                  string("'num_dimensions' (= ") +
                  to_string(num_dimensions) +
                  string(") must be equal to 2 or 3"));
    }

    if (N <= 0) {
        PQS_ERROR_STATIC("lsmTest", "setUpGrid",
                  string("'N' (= ") + to_string(N) +
                  string(") must positive"));
    }

    // --- Set up grid

    // Set lower and upper corners of grid
    x_lo = new double[num_dimensions];
    x_hi = new double[num_dimensions];
    for (int i = 0; i < num_dimensions; i++) {
        x_lo[i] = -1.0;
        x_hi[i] = 1.0;
    }

    // Compute dx
    dx = new double[num_dimensions];
    for (int i = 0; i < num_dimensions; i++) {
        dx[i] = 1.0 * (x_hi[i] - x_lo[i]) / N;
    }

    // Compute lower and upper corners of index space
    phi_gb_lo = new int[num_dimensions];
    phi_gb_hi = new int[num_dimensions];
    ib_lo = new int[num_dimensions];
    ib_hi = new int[num_dimensions];
    for (int i = 0; i < num_dimensions; i++) {
        phi_gb_lo[i] = 1;
        phi_gb_hi[i] = N;
        ib_lo[i] = 1;
        ib_hi[i] = N;
    }

    // Set grid_setup_complete to true
    grid_setup_complete = true;

} // lsmTest::setUpGrid()

void lsmTest::tearDownGrid() {
    // Clean up geometry parameters
    if (x_lo) {
        delete [] x_lo;
        x_lo = NULL;
    }
    if (x_hi) {
        delete [] x_hi;
        x_hi = NULL;
    }

    // Clean up grid parameters
    if (dx) {
        delete [] dx;
        dx = NULL;
    }
    if (phi_gb_lo) {
        delete [] phi_gb_lo;
        phi_gb_lo = NULL;
    }
    if (phi_gb_hi) {
        delete [] phi_gb_hi;
        phi_gb_hi = NULL;
    }
    if (ib_lo) {
        delete [] ib_lo;
        ib_lo = NULL;
    }
    if (ib_hi) {
        delete [] ib_hi;
        ib_hi = NULL;
    }

    // Set grid_setup_complete to false
    grid_setup_complete = false;

} // lsmTest::tearDownGrid()

void lsmTest::setUpGridFunction()
{
    // Check prerequisites
    if (!grid_setup_complete) {
        PQS_ERROR_STATIC("lsmTest", "setUpGridFunction",
            "setUpGrid() must be called before setUpGridFunction()");
    }

    // Free pre-existing memory for phi
    tearDownGridFunction();

    // Allocate memory for phi
    int grid_size = 1;
    for (int i = 0; i < num_dimensions; i++) {
        grid_size *= N;
    }
    phi = new double[grid_size];

} // lsmTest::setUpGridFunction()

void lsmTest::tearDownGridFunction()
{
    // Free memory for phi
    if (phi) {
        delete [] phi;
        phi = NULL;
    }
} // lsmTest::tearDownGridFunction()

} // lsmTests namespace
