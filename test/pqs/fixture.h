/*! \file fixture.h
 *
 * \brief
 * Fixture for unit tests for PQS::pqs classes.
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

#ifndef INCLUDED_PQS_fixture_h
#define INCLUDED_PQS_fixture_h

// --- Headers, namespaces, and type declarations

// Standard library
#include <memory>

// Google Test
#include "gtest/gtest.h"

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/InterfaceInitStrategy.h"  // IWYU pragma: keep
#include "PQS/pqs/PoreInitStrategy.h"  // IWYU pragma: keep

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations
namespace SAMRAI { namespace geom { class CartesianGridGeometry; } }
namespace SAMRAI { namespace hier { class PatchHierarchy; } }
namespace SAMRAI { namespace tbox { class Database; } }
namespace PQS { namespace pqs { class Solver; } }


// --- Fixtures

namespace pqsTests {

class pqsTest: public ::testing::Test
{
protected:
    // --- Fixture data

    // Configuration database
    shared_ptr<tbox::Database> config_db;

    // SAMRAI::geom::Geometry
    shared_ptr<geom::CartesianGridGeometry> grid_geometry;

    // SAMRAI::hier::PatchHierarchy
    shared_ptr<hier::PatchHierarchy> patch_hierarchy;

    // PQS::pqs::PoreInitStrategy
    shared_ptr<pqs::PoreInitStrategy> pore_init_strategy;

    // PQS::pqs::InterfaceInitStrategy
    shared_ptr<pqs::InterfaceInitStrategy> interface_init_strategy;

    // PQS::pqs::Solver
    pqs::Solver *solver;

    // --- Fixture set up and tear down

    // Constructor (set up)
    pqsTest();

    // Destructor (tear down)
    virtual ~pqsTest();

    // --- Helper methods

    /*
     * Set configuration parameters, initialize geometry, and initialize
     * PatchHierarchy based on 'num_dimensions'.
     */
    static void initializeGeometryAndHierarchy(
            shared_ptr<tbox::Database> config_db,
            shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int num_dimensions);

    // --- Test management
    //
    // NOTE: used to ensure that per-process initialization and cleanup
    // are not called multiple times.
    static int s_num_tests;
    static int s_num_tests_remaining;
};

} // pqsTests namespace

#endif // INCLUDED_PQS_fixture_h
