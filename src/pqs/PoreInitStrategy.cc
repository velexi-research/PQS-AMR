/*! \file PoreInitStrategy.cc
 *
 * \brief
 * Implementation file for PoreInitStrategy class.
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
#include "PQS/pqs/PoreInitStrategy.h"

// Class/type declarations


// --- Class implementation

namespace PQS {
namespace pqs {

// --- Public methods

void PoreInitStrategy::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::PoreInitStrategy::printClassData..." << endl;
    os << "(PoreInitStrategy*) this = " << (PoreInitStrategy*) this << endl;
}

} // PQS::pqs namespace
} // PQS namespace
