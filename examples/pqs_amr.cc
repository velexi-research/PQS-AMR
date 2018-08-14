/*! \file pqs_amr.cc
 *
 * \brief
 * Example PQS AMR application.
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

// --- Headers, namespaces, and type declarations

// Standard

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/utilities/application.h"

// Class/type declarations
namespace SAMRAI { namespace tbox { class Database; } }

// Namespaces
using namespace std;
using namespace PQS;
using namespace SAMRAI;


// --- Main program

int main(int argc, char *argv[])
{

    // Initialize PQS simulation
    boost::shared_ptr<tbox::Database> config_db = initialize_pqs(argc, argv);

    // Run PQS simulation
    // run();

    // Shutdown PQS simulation
    shutdown_pqs();

    return(0);
}
