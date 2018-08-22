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
shared_ptr<tbox::Database> initialize_pqs(int *argc, char **argv[]);

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
 * pqs_solver: pqs::Solver object
 *
 * visit_data_writer: VisIt data file writer
 */
void run_pqs(
        const shared_ptr<tbox::Database>& config_db,
        const shared_ptr<pqs::Solver>& pqs_solver,
        const shared_ptr<appu::VisItDataWriter>& visit_data_writer = 0);

}  // PQS namespace

#endif  // INCLUDED_PQS_application_h
