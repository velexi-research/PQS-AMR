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
#include <memory>

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

    // Get patch box
    hier::Box patch_box = patch.getBox();

    // Get patch data
    shared_ptr< pdat::CellData<PQS_REAL> > phi_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch.getPatchData(phi_id));

    PQS_REAL* phi = phi_data->getPointer();

    // --- Initialize phi

    // Set phi to 1 everywhere except for single edge of box
    phi_data->fill(1.0);

    const hier::IntVector patch_box_lower = patch_box.lower();
    const hier::IntVector patch_box_upper = patch_box.upper();
    for (int i = patch_box_lower[0] ; i <= patch_box_upper[0]; i++) {
        phi[i] = 0.0;
    }
}

} // pqsTests namespace
