/*! \file Algorithms.cc
 *
 * \brief
 * Implementation file for Algorithms class.
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

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/Algorithms.h"

// Class/type declarations


// --- Class implementation

namespace PQS {
namespace pqs {


// --- Implementation of public methods

void Algorithms::computePrescribedCurvatureModelRHS()
{
} // Algorithms::computePrescribedCurvatureModelRHS()

void Algorithms::computeSlightlyCompressibleModelRHS()
{
} // Algorithms::computeSlightlyCompressibleModelRHS()

void Algorithms::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::Algorithms::printClassData..." << endl;
    os << "(Algorithms*) this = " << (Algorithms*) this << endl;
} // Algorithms::printClassData()

} // PQS::pqs namespace
} // PQS namespace
