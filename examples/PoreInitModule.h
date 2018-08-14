/*! \file PoreInitModule.h
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

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/PoreInitStrategy.h"

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations
namespace SAMRAI { namespace hier { class Patch; } }


// --- Class declaration

class PoreInitModule: public pqs::PoreInitStrategy
{
public:
    /*
     * Set value at all grid points to 1.0.
     */
    virtual void initializePoreSpace(hier::Patch& patch, int psi_id);
};
