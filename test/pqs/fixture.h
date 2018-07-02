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

// --- Headers, namespaces, and type declarations

// Google Test
#include "gtest/gtest.h"

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/DataInitStrategy.h"  // IWYU pragma: keep

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations
namespace SAMRAI { namespace geom { class CartesianGridGeometry; } }
namespace SAMRAI { namespace hier { class PatchHierarchy; } }
namespace SAMRAI { namespace tbox { class MemoryDatabase; } }


// --- Fixtures

namespace pqsTests {

class pqsTests: public ::testing::Test
{
protected:
    // --- Fixture data

    // Configuration database
    boost::shared_ptr<tbox::MemoryDatabase> config_db;

    // SAMRAI::geom::Geometry
    boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry;

    // SAMRAI::hier::PatchHierarchy
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy;

    // PQS::pqs::DataInitStrategy
    boost::shared_ptr<pqs::DataInitStrategy> data_init_strategy;

    // --- Fixture set up and tear down

    // Constructor (set up)
    pqsTests();

    // Destructor (tear down)
    virtual ~pqsTests();
};

} // pqsTests namespace
