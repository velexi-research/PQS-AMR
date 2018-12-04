/*! \file boundary_conditions_2d.h
 *
 * \brief
 * Header files for numerical kernels for filling boundary data for
 * partial differnetial equations in two space dimensions
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

#ifndef INCLUDED_PQS_math_boundary_conditions_2d_h
#define INCLUDED_PQS_math_boundary_conditions_2d_h

#include "PQS/PQS_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief
 * @ref boundary_conditions_2d.h provides support for filling boundary
 * data for partial differential equations in two space dimensions.
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in             name in
 *      C/C++ code          Fortran code
 *      ----------          ------------
 */
#define PQS_MATH_2D_FILL_EDGE_BDRY_DATA_LINEAR_EXTRAPOLATION \
                            pqsmath2dfilledgebdrydatalinearextrapolation_

#define PQS_MATH_2D_FILL_NODE_BDRY_DATA_LINEAR_EXTRAPOLATION \
                            pqsmath2dfillnodebdrydatalinearextrapolation_


/*!
 * See documentation in boundary_conditions_2d.f.in
 */
void PQS_MATH_2D_FILL_EDGE_BDRY_DATA_LINEAR_EXTRAPOLATION(
    PQS_REAL *u,
    const int *u_gb_lo,
    const int *u_gb_hi,
    const int *interior_box_lo,
    const int *interior_box_hi,
    const int *bdry_box_lo,
    const int *bdry_box_hi,
    const int *bdry_location_idx);

/*!
 * See documentation in boundary_conditions_2d.f.in
 */
void PQS_MATH_2D_FILL_NODE_BDRY_DATA_LINEAR_EXTRAPOLATION(
    PQS_REAL *u,
    const int *u_gb_lo,
    const int *u_gb_hi,
    const int *interior_box_lo,
    const int *interior_box_hi,
    const int *bdry_box_lo,
    const int *bdry_box_hi,
    const int *bdry_location_idx);

#ifdef __cplusplus
}
#endif

#endif
