/*! \file TestDataInitModule.h
 *
 * \brief
 * Concrete implementation of pqs::DataInitStrategy to use for testing.
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

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/DataInitStrategy.h"

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations


// --- Fixtures

namespace pqsTests {

class TestDataInitModule: public pqs::DataInitStrategy
{
public:
    /*
     * No-op. Only used for testing.
     */
    virtual void initializePoreSpace(hier::Patch& patch, int psi_id) {};

    /*
     * No-op. Only used for testing.
     */
    virtual void initializeInterface(hier::Patch& patch, int phi_id) {};
};

} // pqsTests namespace
