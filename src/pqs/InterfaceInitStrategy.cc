/*! \file InterfaceInitStrategy.cc
 *
 * \brief
 * Implementation file for InterfaceInitStrategy class.
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
#include <sstream>

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/InterfaceInitStrategy.h"

// Class/type declarations


// --- Class implementation

namespace PQS {
namespace pqs {

// --- Implementation of public methods

void InterfaceInitStrategy::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::InterfaceInitStrategy::printClassData..." << endl;
    os << "(InterfaceInitStrategy*) this = " << (InterfaceInitStrategy*) this
       << endl;
}

} // PQS::pqs namespace
} // PQS namespace
