/*! \file BoundaryConditions.cc
 *
 * \brief
 * Implementation file for functions to fill boundary data
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
#include <memory>
#include <string>
#include <vector>

// SAMRAI headers
#include "SAMRAI/appu/CartesianBoundaryDefines.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS headers
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/math/BoundaryConditions.h"  // IWYU pragma: keep
#include "PQS/math/kernels/boundary_conditions_2d.h"
#include "PQS/math/kernels/boundary_conditions_3d.h"
#include "PQS/utilities/error.h"
#include "PQS/utilities/macros.h"

// Class/type declarations


// --- Implementation for PQS::math toolbox functions

namespace PQS {
namespace math {

void fillBdryDataLinearExtrapolation(
        const shared_ptr<hier::PatchHierarchy> patch_hierarchy,
        const int u_id)
{
    // --- Check arguments

    const int dim = patch_hierarchy->getDim().getValue();
    if ((dim != 2) && (dim != 3)) {
        PQS_ERROR_STATIC("math", "fillBdryDataLinearExtrapolation",
                         string("Invalid dimension (=") + to_string(dim) +
                         string("for patch_hierarchy. Valid dimensions: 2, 3"));
    }

    // --- Fill boundary data

    // loop over PatchHierarchy fill boundary data on each Patch
    // by calling Fortran subroutines
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    for (int level_num = 0; level_num < num_levels; level_num++) {

        shared_ptr<hier::PatchLevel> level =
                patch_hierarchy->getPatchLevel(level_num);

        for (hier::PatchLevel::Iterator pi(level->begin());
                pi!=level->end(); pi++) {

            // loop over patches
            shared_ptr<hier::Patch> patch = *pi;
            if (patch==NULL) {
                PQS_ERROR_STATIC("math", "fillBdryDataLinearExtrapolation",
                                 "Null patch pointer");
            }

            // get PatchData and PatchGeometry
            shared_ptr< pdat::CellData<PQS_REAL> > u_data =
                    SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                            patch->getPatchData(u_id));

            const shared_ptr<geom::CartesianPatchGeometry> pgeom(
                    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                            patch->getPatchGeometry()));

            // get boundary boxes
            vector<hier::BoundaryBox> bdry_boxes;
            if (dim == 3) {
                bdry_boxes = pgeom->getCodimensionBoundaries(Bdry::FACE3D);

            } else if (dim == 2) {
                bdry_boxes = pgeom->getCodimensionBoundaries(Bdry::EDGE2D);
            }

            // skip Patch if it does not touch physical boundaries
            if (bdry_boxes.size() == 0) {
                continue;
            }

            // get index space ranges
            hier::Box u_ghostbox = u_data->getGhostBox();
            const hier::IntVector u_ghostbox_lower = u_ghostbox.lower();
            const hier::IntVector u_ghostbox_upper = u_ghostbox.upper();
            PQS_INT_VECT_TO_INT_ARRAY(u_ghostbox_lo, u_ghostbox_lower);
            PQS_INT_VECT_TO_INT_ARRAY(u_ghostbox_hi, u_ghostbox_upper);

            // interior box
            hier::Box interior_box = patch->getBox();
            const hier::IntVector interior_box_lower = interior_box.lower();
            const hier::IntVector interior_box_upper = interior_box.upper();
            PQS_INT_VECT_TO_INT_ARRAY(interior_box_lo, interior_box_lower);
            PQS_INT_VECT_TO_INT_ARRAY(interior_box_hi, interior_box_upper);

            // get pointer to data
            PQS_REAL* u = u_data->getPointer();

            // fill boundary data
            const hier::IntVector& ghostcell_width(u_data->getGhostCellWidth());
            for (int i = 0; i < static_cast<int>(bdry_boxes.size()); i++) {

                const int bdry_location_idx =
                    bdry_boxes[i].getLocationIndex();

                // get boundary box
                hier::Box bdry_box(pgeom->getBoundaryFillBox(
                        bdry_boxes[i], interior_box, ghostcell_width));
                const hier::IntVector bdry_box_lower = bdry_box.lower();
                const hier::IntVector bdry_box_upper = bdry_box.upper();
                PQS_INT_VECT_TO_INT_ARRAY(bdry_box_lo, bdry_box_lower);
                PQS_INT_VECT_TO_INT_ARRAY(bdry_box_hi, bdry_box_upper);

                if (dim == 3) {
                    PQS_MATH_3D_FILL_BDRY_DATA_LINEAR_EXTRAPOLATION(
                        u, u_ghostbox_lo, u_ghostbox_hi,
                        interior_box_lo, interior_box_hi,
                        bdry_box_lo, bdry_box_hi,
                        &bdry_location_idx);
                } else if (dim == 2) {
                    PQS_MATH_2D_FILL_BDRY_DATA_LINEAR_EXTRAPOLATION(
                        u, u_ghostbox_lo, u_ghostbox_hi,
                        interior_box_lo, interior_box_hi,
                        bdry_box_lo, bdry_box_hi,
                        &bdry_location_idx);
                }
            }

        } // end loop over patches in level
    } // end loop over levels in hierarchy
}

} // PQS::math namespace
} // PQS namespace
