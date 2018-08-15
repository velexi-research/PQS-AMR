/*! \file PoreInitModule.cc
 *
 * \brief
 * Concrete implementation of pqs::PoreInitStrategy to use for TODO example.
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

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep

// PQS test
#include "PoreInitModule.h"

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchData; } }


// --- Class implementation

void PoreInitModule::initializePoreSpace(
        hier::Patch& patch, int psi_id)
{
    boost::shared_ptr< pdat::CellData<PQS_REAL> > psi_data =
            BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch.getPatchData(psi_id));

    psi_data->fill(1.0);
}
