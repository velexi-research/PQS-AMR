/*! \file fixture.h
 *
 * \brief
 * Fixture for unit tests for PQS::lsm package
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

#ifndef INCLUDED_PQS_lsm_fixture_h
#define INCLUDED_PQS_lsm_fixture_h

// --- Headers, namespaces, and type declarations

// Standard library

// Google Test
#include "gtest/gtest.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep

// Namespaces
using namespace std;

// Class/type declarations


// --- Fixtures

namespace lsmTests {

class lsmTest: public ::testing::Test
{
protected:
    // --- Fixture data

    // geometry
    int num_dimensions; // number of spatial dimensions for grid
    double *x_lo; // lower corner of grid
    double *x_hi; // upper corner of grid

    // grid parameters
    int N; // number of grid cells in each coordinate direction
    double *dx; // dimensions of grid cell in each coordinate direction
    int *phi_gb_lo; // lower corner of ghostbox
    int *phi_gb_hi; // upper corner of ghostbox
    int *ib_lo; // lower corner of interior box
    int *ib_hi; // upper corner of interior box

    // grid functions
    PQS_REAL *phi;

    // state variables
    bool grid_setup_complete;

    // --- Fixture set up and tear down

    // Constructor (set up)
    lsmTest();

    // Destructor (tear down)
    //
    // Free heap memory used by unit tests.
    virtual ~lsmTest();

    // --- Helper methods

    /*
     * Set up numerical grid.
     *
     * Notes
     * -----
     * - num_dimensions, the number of dimensions in the grid, and
     *   N, the number of grid cells in each coordinate direction,
     *   must be set before calling setUpGrid().
     */
    void setUpGrid();

    /*
     * Tear down numerical grid.
     *
     * Notes
     * -----
     * - tearDownGrid() does not free memory used for grid functions.
     */
    void tearDownGrid();

    /*
     * Set up grid function.
     *
     * Notes
     * -----
     * - setUpGrid() must be called before setUpGridFunction().
     */
    void setUpGridFunction();

    /*
     * Tear down grid function. Free memory allocated for grid function.
     */
    void tearDownGridFunction();

    // --- Test management
    //
    // NOTE: used to ensure that per-process initialization and cleanup
    // are not called multiple times.
    static int s_num_tests;
    static int s_num_tests_remaining;
};

} // lsmTests namespace

#endif // INCLUDED_PQS_lsm_fixture_h
