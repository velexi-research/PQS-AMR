/*! \file TimeIntegration.cc
 *
 * \brief
 * Implementation file for TimeIntegration class.
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

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/hier/PatchHierarchy.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/utilities.h"
#include "PQS/math/TimeIntegration.h"

// Class/type declarations


// --- Class implementation

namespace PQS {
namespace math {

// --- Implementation of public methods

void TimeIntegration::forwardEulerStep(
            boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_next_id,
            const int u_current_id,
            const int rhs_id,
            const double dt)
{
    // Check arguments
    if (patch_hierarchy == NULL) {
        PQS_ERROR_STATIC("forwardEulerStep",
                         "'patch_hierarchy' must not be NULL");
    }
    if (u_next_id < 0) {
        PQS_ERROR_STATIC("forwardEulerStep",
                         "'u_next_id' must be non-negative");
    }
    if (u_current_id < 0) {
        PQS_ERROR_STATIC("forwardEulerStep",
                         "'u_current_id' must be non-negative");
    }
    if (rhs_id < 0) {
        PQS_ERROR_STATIC("forwardEulerStep",
                         "'rhs_id' must be non-negative");
    }
    if (dt <= 0) {
        PQS_ERROR_STATIC("forwardEulerStep",
                         "'dt' must be positive");
    }

    // Advance solution by time step dt
}

} // PQS::math namespace
} // PQS namespace
