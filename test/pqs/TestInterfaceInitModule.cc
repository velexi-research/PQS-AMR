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
    shared_ptr< pdat::CellData<PQS_REAL> > phi_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch.getPatchData(phi_id));

    phi_data->fill(1.0);
}

} // pqsTests namespace
