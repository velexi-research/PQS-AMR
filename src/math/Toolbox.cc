/*! \file Toolbox.cc
 *
 * \brief
 * Implementation file for math toolbox
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

// External package headers
#include "mpi.h"

// SAMRAI headers
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS headers
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/math/Toolbox.h"
#include "PQS/math/level_set_method_2d.h"
#include "PQS/math/level_set_method_3d.h"
#include "PQS/math/utilities_2d.h"
#include "PQS/math/utilities_3d.h"
#include "PQS/utilities/error.h"
#include "PQS/utilities/macros.h"

// Class/type declarations


// --- Implementation for PQS::math toolbox functions

namespace PQS {
namespace math {

PQS_REAL maxNormDiff(
        const shared_ptr<hier::PatchHierarchy> patch_hierarchy,
        const int u_id,
        const int v_id,
        const int control_volume_id)
{
    // --- Check arguments

    const int dim = patch_hierarchy->getDim().getValue();
    if ((dim != 2) && (dim != 3)) {
        PQS_ERROR_STATIC("math", "maxNormDiff",
                         string("Invalid dimension (=") + to_string(dim) +
                         string("for patch_hierarchy. Valid dimensions: 2, 3"));
    }

    // --- Compute max norm of (u - v)

    PQS_REAL max_norm_diff = 0.0;

    // loop over PatchHierarchy and compute the integral on each Patch
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
                PQS_ERROR_STATIC("math", "maxNormDiff",
                                 "Null patch pointer");
            }

            // get pointers to data and index space ranges
            shared_ptr< pdat::CellData<PQS_REAL> > u_data =
                    SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                            patch->getPatchData( u_id ));
            shared_ptr< pdat::CellData<PQS_REAL> > control_volume_data =
                    SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                        patch->getPatchData( control_volume_id ));

            hier::Box u_ghostbox = u_data->getGhostBox();
            const hier::IntVector u_ghostbox_lower = u_ghostbox.lower();
            const hier::IntVector u_ghostbox_upper = u_ghostbox.upper();
            PQS_INT_VECT_TO_INT_ARRAY(u_ghostbox_lo, u_ghostbox_lower);
            PQS_INT_VECT_TO_INT_ARRAY(u_ghostbox_hi, u_ghostbox_upper);

            shared_ptr< pdat::CellData<PQS_REAL> > v_data =
                    SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                            patch->getPatchData( v_id ));

            hier::Box v_ghostbox = v_data->getGhostBox();
            const hier::IntVector v_ghostbox_lower = v_ghostbox.lower();
            const hier::IntVector v_ghostbox_upper = v_ghostbox.upper();
            PQS_INT_VECT_TO_INT_ARRAY(v_ghostbox_lo, v_ghostbox_lower);
            PQS_INT_VECT_TO_INT_ARRAY(v_ghostbox_hi, v_ghostbox_upper);

            hier::Box control_volume_ghostbox =
                    control_volume_data->getGhostBox();
            const hier::IntVector control_volume_ghostbox_lower =
                    control_volume_ghostbox.lower();
            const hier::IntVector control_volume_ghostbox_upper =
                    control_volume_ghostbox.upper();
            PQS_INT_VECT_TO_INT_ARRAY(control_volume_ghostbox_lo,
                                      control_volume_ghostbox_lower);
            PQS_INT_VECT_TO_INT_ARRAY(control_volume_ghostbox_hi,
                                      control_volume_ghostbox_upper);

            // interior box
            hier::Box interior_box = patch->getBox();
            const hier::IntVector interior_box_lower = interior_box.lower();
            const hier::IntVector interior_box_upper = interior_box.upper();
            PQS_INT_VECT_TO_INT_ARRAY(interior_box_lo, interior_box_lower);
            PQS_INT_VECT_TO_INT_ARRAY(interior_box_hi, interior_box_upper);

            PQS_REAL* u = u_data->getPointer();
            PQS_REAL* v = v_data->getPointer();
            PQS_REAL* control_volume = control_volume_data->getPointer();
            PQS_REAL max_norm_diff_on_patch = 0.0;

            if ( dim == 3 ) {
                max_norm_diff_on_patch = PQS_MATH_3D_MAX_NORM_DIFF(
                        u, u_ghostbox_lo, u_ghostbox_hi,
                        v, v_ghostbox_lo, v_ghostbox_hi,
                        interior_box_lo, interior_box_hi);
            } else if ( dim == 2 ) {
                max_norm_diff_on_patch = PQS_MATH_2D_MAX_NORM_DIFF(
                        u, u_ghostbox_lo, u_ghostbox_hi,
                        v, v_ghostbox_lo, v_ghostbox_hi,
                        interior_box_lo, interior_box_hi);
            }

            if (max_norm_diff < max_norm_diff_on_patch) {
                max_norm_diff = max_norm_diff_on_patch;
            }

        } // end loop over patches in level
    } // end loop over levels in hierarchy

    int status =
            tbox::SAMRAI_MPI::getSAMRAIWorld().AllReduce(&max_norm_diff, 1,
                                                         MPI_MAX);

    if (status != 0) {
        PQS_ERROR_STATIC("math", "maxNormDiff",
                         string("AllReduce error code=") + to_string(status));
    }

    return max_norm_diff;
}


PQS_REAL computeVolume(
        const shared_ptr<hier::PatchHierarchy> patch_hierarchy,
        const int phi_id,
        const int control_volume_id,
        const int region_indicator)
{
    // --- Check arguments

    const int dim = patch_hierarchy->getDim().getValue();
    if ((dim != 2) && (dim != 3)) {
        PQS_ERROR_STATIC("math::LSMToolbox", "computeVolume",
                         string("Invalid dimension (=") + to_string(dim) +
                         string("for patch_hierarchy. Valid dimensions: 2, 3"));
    }

    // --- Compute volume

    PQS_REAL volume = 0.0;

    // loop over PatchHierarchy and compute the integral on each Patch
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
                PQS_ERROR_STATIC("math::LSMToolbox", "computeVolume",
                                 "Null patch pointer");
            }

            // get dx and compute epsilon
            shared_ptr<geom::CartesianPatchGeometry> patch_geom =
                    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                            patch->getPatchGeometry());
            const double* dx = patch_geom->getDx();

            double max_dx = dx[0];
            for (int i = 1; i < dim; i++) {
                if (max_dx < dx[i]) max_dx = dx[i];
            }

            double epsilon = 1.5 * max_dx;

            // get pointers to data and index space ranges
            shared_ptr< pdat::CellData<PQS_REAL> > phi_data =
                    SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                            patch->getPatchData( phi_id ));
            shared_ptr< pdat::CellData<PQS_REAL> > control_volume_data =
                    SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                        patch->getPatchData( control_volume_id ));

            hier::Box phi_ghostbox = phi_data->getGhostBox();
            const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
            const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();
            PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_lo, phi_ghostbox_lower);
            PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_hi, phi_ghostbox_upper);

            hier::Box control_volume_ghostbox =
                    control_volume_data->getGhostBox();
            const hier::IntVector control_volume_ghostbox_lower =
                    control_volume_ghostbox.lower();
            const hier::IntVector control_volume_ghostbox_upper =
                    control_volume_ghostbox.upper();
            PQS_INT_VECT_TO_INT_ARRAY(control_volume_ghostbox_lo,
                                      control_volume_ghostbox_lower);
            PQS_INT_VECT_TO_INT_ARRAY(control_volume_ghostbox_hi,
                                      control_volume_ghostbox_upper);

            // interior box
            hier::Box interior_box = patch->getBox();
            const hier::IntVector interior_box_lower = interior_box.lower();
            const hier::IntVector interior_box_upper = interior_box.upper();
            PQS_INT_VECT_TO_INT_ARRAY(interior_box_lo, interior_box_lower);
            PQS_INT_VECT_TO_INT_ARRAY(interior_box_hi, interior_box_upper);

            PQS_REAL* phi = phi_data->getPointer();
            PQS_REAL* control_volume = control_volume_data->getPointer();
            PQS_REAL volume_on_patch = 0.0;
            int control_volume_sgn = 1;

            if ( dim == 3 ) {
                if (region_indicator > 0) {
                    // integrate over region {x | phi(x) > 0}
                    volume_on_patch = LSM_3D_VOLUME_PHI_GREATER_THAN_ZERO(
                            phi, phi_ghostbox_lo, phi_ghostbox_hi,
                            interior_box_lo, interior_box_hi,
                            dx, &epsilon);
                } else {
                    // integrate over region {x | phi(x) < 0}
                    volume_on_patch = LSM_3D_VOLUME_PHI_LESS_THAN_ZERO(
                            phi, phi_ghostbox_lo, phi_ghostbox_hi,
                            interior_box_lo, interior_box_hi,
                            dx, &epsilon);
                }
            } else if ( dim == 2 ) {
                if (region_indicator > 0) {
                    // integrate over region {x | phi(x) > 0}
                    volume_on_patch = LSM_2D_AREA_PHI_GREATER_THAN_ZERO(
                            phi, phi_ghostbox_lo, phi_ghostbox_hi,
                            interior_box_lo, interior_box_hi,
                            dx, &epsilon);
                } else {
                    // integrate over region {x | phi(x) < 0}
                    volume_on_patch = LSM_2D_AREA_PHI_LESS_THAN_ZERO(
                            phi, phi_ghostbox_lo, phi_ghostbox_hi,
                            interior_box_lo, interior_box_hi,
                            dx, &epsilon);
                }
            }

            volume += volume_on_patch;

        } // end loop over patches in level
    } // end loop over levels in hierarchy

    int status =
            tbox::SAMRAI_MPI::getSAMRAIWorld().AllReduce(&volume,1, MPI_SUM);

    if (status != 0) {
        PQS_ERROR_STATIC("math::LSMToolbox", "computeVolume",
                         string("AllReduce error code=") + to_string(status));
    }

    return volume;
}

} // PQS::math namespace
} // PQS namespace
