/*! \file application.h
 *
 * \brief
 * Header file for PQS simulation application utility functions.
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

#ifndef INCLUDED_PQS_application_h
#define INCLUDED_PQS_application_h

// --- Headers, namespaces, and type declarations

// Standard
#include <iosfwd>
#include <iostream>
#include <stdlib.h>
#include <string>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

// PQS
#include "PQS/PQS_config.h"
#include "PQS/pqs/Solver.h"
//#include "PQS/lsm/Toolbox.h"

// Class/type declarations

// Namespaces
using namespace std;
using namespace PQS;
using namespace SAMRAI;


namespace PQS {

// --- Function declarations

/*!
 * Initialize PQS simulation.
 *
 * Parameters
 * ----------
 * argc: number of arguments passed to 'main()'
 *
 * argv: arguments passed to 'main()'
 *
 * Return value
 * ------------
 * database containing configuration parameters for simulation
 */
boost::shared_ptr<tbox::Database> initialize_pqs(int argc, char *argv[]);

/*!
 * Shut down PQS simulation.
 */
void shutdown_pqs();

/*!
 * Run PQS simulation.
 *
 * Parameters
 * ----------
 * config_db: database containing configuration parameters
 *
 * pore_init_strategy: user implementation of PoreInitStrategy strategy
 *      interface
 *
 * interface_init_strategy: user implementation of
 *      InterfaceInitStrategy strategy interface
 */
void run_pqs(
        const boost::shared_ptr<tbox::Database> config_db,
        const boost::shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const boost::shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy);

}  // PQS namespace

#endif  // INCLUDED_PQS_application_h
