/*! \file TestInterfaceInitModule.h
 *
 * \brief
 * Header file for concrete subclass of pqs::InterfaceInitStrategy to use for
 * testing.
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
#include "PQS/pqs/InterfaceInitStrategy.h"

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations
namespace SAMRAI { namespace hier { class Patch; } }


// --- Fixtures

namespace pqsTests {

class TestInterfaceInitModule: public pqs::InterfaceInitStrategy
{
public:
    /*
     * Set value at all grid points to 1.0.
     */
    virtual void initializeInterface(hier::Patch& patch, int phi_id);
};

} // pqsTests namespace
