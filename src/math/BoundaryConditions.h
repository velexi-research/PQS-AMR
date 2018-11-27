/*! \file BoundaryConditions.h
 *
 * \brief
 * Header file for functions to fill boundary data
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

#ifndef INCLUDED_PQS_math_BoundaryConditions_h
#define INCLUDED_PQS_math_BoundaryConditions_h

/*! PQS::math functions
 *
 * \brief
 * TODO: add description
 */

// --- Headers, namespaces, and type declarations

// Standard headers
#include <ostream>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/PatchHierarchy.h"

// PQS headers
#include "PQS/PQS_config.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations


// --- PQS::math functions

namespace PQS {
namespace math {

/*!
 * Fill boundary data using linear extrapolation from the interior
 * of the domain.
 *
 * Parameters
 * ----------
 * patch_hierarchy: PatchHierarchy to perform computation on
 *
 * u_id: PatchData id of variable to fill boundary data for
 */
void fillBdryDataLinearExtrapolation(
        const shared_ptr<hier::PatchHierarchy> patch_hierarchy,
        const int u_id);

}  // PQS::math namespace
}  // PQS namespace

#endif // INCLUDED_PQS_math_BoundaryConditions_h
