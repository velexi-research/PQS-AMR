/*! \file Algorithms.cc
 *
 * \brief
 * Implementation file for Algorithms class.
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

// Standard library
#include <sstream>
#include <stddef.h>
#include <string>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/hier/Patch.h"
//#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/utilities.h"
#include "PQS/pqs/Algorithms.h"

// Class/type declarations


// --- Class implementation

namespace PQS {
namespace pqs {


// --- Implementation of public methods


Algorithms::Algorithms(
        const boost::shared_ptr<tbox::Database>& config_db,
        const int psi_id,
        const int grad_psi_id)
{
    // Check arguments
    if (config_db == NULL) {
        PQS_ERROR(this, "Algorithms", "'config_db' must not be NULL");
    }
    verifyConfigurationDatabase(config_db);

    if (psi_id < 0) {
        PQS_ERROR(this, "Algorithms", "'psi_id' must be non-negative");
    }

    // Set data members
    d_psi_id = psi_id;
    d_grad_psi_id = grad_psi_id;

} // Algorithms::Algorithms()

double Algorithms::computePrescribedCurvatureModelRHS(
        const boost::shared_ptr<hier::Patch> patch,
        const int phi_id) const
{
    patch->getPatchData(phi_id);
    return 1.0;
} // Algorithms::computePrescribedCurvatureModelRHS()

double Algorithms::computeSlightlyCompressibleModelRHS(
        const boost::shared_ptr<hier::Patch> patch,
        const int phi_id) const
{
    patch->getPatchData(phi_id);
    return 1.0;
} // Algorithms::computeSlightlyCompressibleModelRHS()

void Algorithms::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::Algorithms::printClassData..." << endl;
    os << "(Algorithms*) this = " << (Algorithms*) this << endl;
} // Algorithms::printClassData()

void Algorithms::verifyConfigurationDatabase(
    const boost::shared_ptr<tbox::Database>& config_db) const
{
    if (config_db == NULL) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'config_db' must not be NULL");
    }

    // --- Verify "SlightlyCompressibleModel" database

    if (!config_db->isDatabase("SlightlyCompressibleModel")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  std::string("'SlightlyCompressibleModel' database ") +
                  std::string("missing from 'config_db'"));
    }
    boost::shared_ptr<tbox::Database> scm_config_db =
            config_db->getDatabase("SlightlyCompressibleModel");

    if (!scm_config_db->isDouble("P_reference")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'P_reference' missing from 'PQS' database");
    }
    if (!scm_config_db->isDouble("V_target")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'V_target' missing from 'PQS' database");
    }
    if (!scm_config_db->isDouble("surface_tension")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'surface_tension' missing from 'PQS' database");
    }
    if (!scm_config_db->isDouble("bulk_modulus")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'bulk_modulus' missing from 'PQS' database");
    }

} // Algorithms::verifyConfigurationDatabase()

void Algorithms::loadConfiguration(
        const boost::shared_ptr<tbox::Database>& config_db)
{
    // TODO
} // Algorithms::loadConfiguration()

} // PQS::pqs namespace
} // PQS namespace
