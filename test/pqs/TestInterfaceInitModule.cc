/*! \file TestInterfaceInitModule.cc
 *
 * \brief
 * Implementation file for concrete subclass of pqs::InterfaceInitStrategy to
 * use for testing.
 */

/*
 * ---------------------------------------------------------------------
 * COPYRIGHT/LICENSE. This file is part of the XYZ package. It is
 * subject to the license terms in the LICENSE file found in the
 * top-level directory of this distribution. No part of the XYZ
 * package, including this file, may be copied, modified, propagated,
 * or distributed except according to the terms contained in the
 * LICENSE file.
 * ---------------------------------------------------------------------
 */

// --- Headers, namespaces, and type declarations

// Standard library
#include <cmath>
#include <memory>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep

// PQS test
#include "TestInterfaceInitModule.h"

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations


// --- Fixtures

namespace pqsTests {

void TestInterfaceInitModule::initializeInterface(
        hier::Patch& patch, int phi_id)
{
    // --- Preparations

    // Get geometry
    shared_ptr<geom::CartesianPatchGeometry> patch_geom =
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                    patch.getPatchGeometry());
    const double* dx = patch_geom->getDx();
    const double* X_lower = patch_geom->getXLower();
    const int dim = patch_geom->getDim().getValue();

    // Get patch box
    hier::Box patch_box = patch.getBox();
    const hier::IntVector patch_box_lower = patch_box.lower();
    const hier::IntVector patch_box_upper = patch_box.upper();

    // Get patch data
    shared_ptr< pdat::CellData<PQS_REAL> > phi_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch.getPatchData(phi_id));

    const hier::IntVector ghostbox_lower = phi_data->getGhostBox().lower();
    const hier::IntVector ghostbox_upper = phi_data->getGhostBox().upper();

    PQS_REAL* phi = phi_data->getPointer();

    // --- Initialize phi to signed distance to circle of radius 0.25

    if (dim == 2) {
        for (int j = patch_box_lower[1] ; j <= patch_box_upper[1]; j++) {
            for (int i = patch_box_lower[0] ; i <= patch_box_upper[0]; i++) {
                double x = X_lower[0] + dx[0] * (i - patch_box_lower[0] + 0.5);
                double y = X_lower[1] + dx[1] * (j - patch_box_lower[1] + 0.5);
                int idx = (i - ghostbox_lower[0]) +
                          (j - ghostbox_lower[1]) *
                              (ghostbox_upper[0] - ghostbox_lower[0]);
                phi[idx] = sqrt(x*x + y*y) - 0.25;
            }
        }
    } else if (dim == 3) {
        for (int k = patch_box_lower[2] ; k <= patch_box_upper[2]; k++) {
            for (int j = patch_box_lower[1] ; j <= patch_box_upper[1]; j++) {
                for (int i = patch_box_lower[0] ; i <= patch_box_upper[0]; i++) {
                    double x = X_lower[0] + dx[0] * (i - patch_box_lower[0] + 0.5);
                    double y = X_lower[1] + dx[1] * (j - patch_box_lower[1] + 0.5);
                    double z = X_lower[2] + dx[2] * (k - patch_box_lower[2] + 0.5);
                    int idx = (i - ghostbox_lower[0]) +
                              (j - ghostbox_lower[1]) *
                                  (ghostbox_upper[0] - ghostbox_lower[0]) +
                              (k - ghostbox_lower[2]) *
                                  (ghostbox_upper[1] - ghostbox_lower[1]) *
                                  (ghostbox_upper[0] - ghostbox_lower[0]);
                    phi[idx] = sqrt(x*x + y*y + z*z) - 0.25;
                }
            }
        }
    }
}

} // pqsTests namespace
