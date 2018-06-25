/*! \file DataInitStrategy.cc
 *
 * \brief
 * Implementation file for DataInitStrategy class.
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
#include <cstddef>                         // for NULL
#include <sstream>                         // for operator<<, basic_ostream

// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>  // for shared_ptr, operator==

// SAMRAI headers

// PQS headers
#include "PQS/PQS_config.h"                // IWYU pragma: keep
#include "PQS/pqs/DataInitStrategy.h"                // for DataInitStrategy

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchHierarchy; } }


// --- Implementation for PQS::pqs::DataInitStrategy methods

namespace PQS {
namespace pqs {

// Constructor
DataInitStrategy::DataInitStrategy(boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy) {
    // Check parameters
    if (patch_hierarchy == NULL) {
        // TODO
    }

    d_patch_hierarchy = patch_hierarchy;
}

// printClassData()
void DataInitStrategy::printClassData(ostream& os) const
{
    os << endl
       << "===================================" << endl;
    os << "PQS::DataInitStrategy" << endl;

    os << "Object Pointers" << endl;
    os << "---------------" << endl;
    os << "(DataInitStrategy*) this = " << (DataInitStrategy*) this << endl;
    os << "d_patch_hierarchy = " << d_patch_hierarchy.get() << endl;

    os << "===================================" << endl << endl;
}

} // PQS::pqs namespace
} // PQS namespace
