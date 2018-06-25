/*! \file GridManager.cc
 *
 * \brief
 * Implementation file for GridManager class.
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

// Standard library headers
#include <cstddef>
#include <sstream>

// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI headers
#include "SAMRAI/mesh/TagAndInitializeStrategy.h"

// PQS headers
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/GridManager.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchHierarchy; } }
namespace SAMRAI { namespace tbox { class Database; } }


// --- Implementation for PQS::pqs::GridManager methods

namespace PQS {
namespace pqs {

// --- Implementation of public PQS::pqs::GridManager methods

// Constructor
GridManager::GridManager(
    boost::shared_ptr<tbox::Database> config_db,
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy):
    mesh::TagAndInitializeStrategy(d_object_name)
{
    // Check parameters
    if (patch_hierarchy == NULL) {
        // TODO
    }

    d_patch_hierarchy = patch_hierarchy;
}

void GridManager::printClassData(ostream& os) const
{
    os << endl
       << "===================================" << endl;
    os << "PQS::GridManager" << endl;

    os << "Object Pointers" << endl;
    os << "---------------" << endl;
    os << "(GridManager*) this = " << (GridManager*) this << endl;
    os << "d_patch_hierarchy = " << d_patch_hierarchy.get() << endl;

    os << "===================================" << endl << endl;
}

// --- Implementation of private PQS::pqs::GridManager methods

// Copy constructor
GridManager::GridManager(const GridManager& rhs):
    mesh::TagAndInitializeStrategy(d_object_name)
{}

} // PQS::pqs namespace
} // PQS namespace
