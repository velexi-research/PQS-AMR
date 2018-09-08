/*! \file Toolbox.h
 *
 * \brief
 * Header file for math toolbox
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

#ifndef INCLUDED_PQS_math_Toolbox_h
#define INCLUDED_PQS_math_Toolbox_h

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
#include "SAMRAI/mesh/GriddingAlgorithm.h"

// PQS headers
#include "PQS/PQS_config.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations


// --- PQS::math functions

namespace PQS {
namespace math {

//! @{
/*!
 ************************************************************************
 *
 * @name Utility functions
 *
 ************************************************************************/

/*!
 * Compute max norm of (u - v).
 *
 * Parameters
 * ----------
 * patch_hierarchy: PatchHierarchy on which to perform computation
 *
 * u_id: PatchData id for 'u'
 *
 * v_id: PatchData id for 'v'
 *
 * control_volume_id: PatchData id for control volume
 *
 * Return value
 * ------------
 * max norm of (u - v)
 */
PQS_REAL maxNormDiff(
        const shared_ptr<hier::PatchHierarchy> patch_hierarchy,
        const int u_id,
        const int v_id,
        const int control_volume_id);

/*!
 * Compute the volume of one of two following regions:
 * region with phi < 0 or region with phi > 0.
 *
 * Parameters
 * ----------
 * patch_hierarchy: PatchHierarchy on which to compute the volume
 *
 * phi_id: PatchData id for phi
 *
 * control_volume_id: PatchData id for control volume
 *
 * region_indicator: integer indicating which region to integrate over
 *      - region_indicator >0:  integration region = {x | phi(x) > 0}
 *      - region_indicator <=0: integration region = {x | phi(x) <= 0}
 *
 * Return value
 * ------------
 * volume of the specified domain that is enclosed by the zero level set
 *
 * Notes
 * -----
 * - When phi is a signed distance function, the smoothed Heaviside
 *   function used compute the integral has a width approximately equal
 *   to the three grid cells widths.
 */
PQS_REAL computeVolume(
        const shared_ptr<hier::PatchHierarchy> patch_hierarchy,
        const int phi_id,
        const int control_volume_id,
        const int region_indicator);

}  // PQS::math namespace
}  // PQS namespace

#endif
