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

// PQS Application
#include "InterfaceInitModule.h"
#include "PoreInitModule.h"

// Class/type declarations
namespace SAMRAI { namespace tbox { class Database; } }
namespace PQS { namespace pqs { class InterfaceInitStrategy; } }
namespace PQS { namespace pqs { class PoreInitStrategy; } }

// Namespaces
using namespace std;
using namespace PQS;
using namespace SAMRAI;


// --- Main program

int main(int argc, char *argv[])
{

    // Initialize PQS simulation
    boost::shared_ptr<tbox::Database> config_db = initialize_pqs(argc, argv);

    // Construct fluid-fluid interface and pore interface initialization
    // modules

    // TestPoreInitModule (implements PoreInitStrategy)
    boost::shared_ptr<pqs::PoreInitStrategy> pore_init_strategy =
            boost::shared_ptr<pqs::PoreInitStrategy>(
                    new PoreInitModule());

    // TestInterfaceInitModule (implements InterfaceInitStrategy)
    boost::shared_ptr<pqs::InterfaceInitStrategy> interface_init_strategy =
            boost::shared_ptr<pqs::InterfaceInitStrategy>(
                    new InterfaceInitModule());

    // Run PQS simulation
    run_pqs(config_db, pore_init_strategy, interface_init_strategy);

    // Shutdown PQS simulation
    shutdown_pqs();

    return(0);
}
