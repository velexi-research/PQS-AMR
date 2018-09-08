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
#include <iosfwd>
#include <memory>
#include <string>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/utilities/error.h"
#include "PQS/utilities/macros.h"

// PQS test
#include "TestInterfaceInitModule.h"
#include "kernels/kernels_2d.h"
#include "kernels/kernels_3d.h"

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations


// --- Fixtures

namespace lsmTests {

TestInterfaceInitModule::TestInterfaceInitModule(
        const int num_dimensions, const double radius)
{
    // --- Check arguments
    if ((num_dimensions != 2) && (num_dimensions != 3)) {
        PQS_ERROR(this, "setPhi",
                  std::string("'num_dimensions' (= ") +
                  std::to_string(num_dimensions) +
                  std::string(") must be equal to 2 or 3"));
    }
    if (radius <= 0) {
        PQS_ERROR(this, "setPhi",
                  std::string("'radius' (= ") +
                  std::to_string(radius) +
                  std::string(") must be greater than zero"));
    }

    // --- Initialize data members

    d_num_dimensions = num_dimensions;
    d_radius = radius;
}

void TestInterfaceInitModule::initializeInterface(
        hier::Patch& patch, int phi_id)
{
    // --- Get geometry parameters

    shared_ptr<geom::CartesianPatchGeometry> patch_geom =
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                    patch.getPatchGeometry());
    const double* dx = patch_geom->getDx();
    const double* x_lo = patch_geom->getXLower();

    // --- Get pointers to data and index space ranges

    shared_ptr< pdat::CellData<PQS_REAL> > phi_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch.getPatchData(phi_id));

    // phi
    hier::Box phi_ghostbox = phi_data->getGhostBox();
    const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
    const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_lo, phi_ghostbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_hi, phi_ghostbox_upper);

    PQS_REAL* phi = phi_data->getPointer();

    // fill box
    hier::Box fillbox = patch.getBox();
    const hier::IntVector fillbox_lower = fillbox.lower();
    const hier::IntVector fillbox_upper = fillbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(fillbox_lo, fillbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(fillbox_hi, fillbox_upper);

    // --- Set phi

    if (d_num_dimensions == 2) {
        SET_PHI_CIRCLE(phi, phi_ghostbox_lo, phi_ghostbox_hi,
                       fillbox_lo, fillbox_hi,
                       x_lo, dx, &d_radius);
    } else if (d_num_dimensions == 3) {
        SET_PHI_SPHERE(phi, phi_ghostbox_lo, phi_ghostbox_hi,
                       fillbox_lo, fillbox_hi,
                       x_lo, dx, &d_radius);
    }
}

} // lsmTests namespace
