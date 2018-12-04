/*! \file PoreInitModule.h
 *
 * \brief
 * Header file for concrete subclass of pqs::PoreInitStrategy to use in
 * example application with pore space defined by a simple packing of
 * circles (2D) or spheres (3D).
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
#include <vector>

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
namespace SAMRAI { namespace tbox { class Database; } }


// --- Class declaration

class PoreInitModule: public pqs::PoreInitStrategy
{
public:
    //! @{

    /*!
     ************************************************************************
     *
     * @name Constructor and destructor
     *
     ************************************************************************/

    /*!
     * The constructor loads configuration parameters and allocates
     * memory required by the object.
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     *
     * dim: problem dimension
     *
     * X_lower: lower corner of computational domain
     *
     * X_upper: upper corner of computational domain
     */
    PoreInitModule(
            const shared_ptr<tbox::Database>& config_db,
            const int dim,
            const vector<double>& X_lower,
            const vector<double>& X_upper);

    /*!
     * Default destructor frees memory allocated for the object.
     */
    virtual ~PoreInitModule();

    //! @}

    //! {@

    /*!
     ************************************************************************
     *
     * @name Pore space initialization methods
     *
     ************************************************************************/

    /*
     * Set value of psi to -1.0 at all grid points (i.e., entire
     * region is pore space).
     */
    virtual void initializePoreSpace(hier::Patch& patch, int psi_id);

    //! @}

protected:

    /****************************************************************
     *
     * Data Members
     *
     ****************************************************************/

    // --- Parameters

    int d_dim;  // problem dimension
    double *d_X_lower;  // lower corner of computational domain
    double *d_X_upper;  // upper corner of computational domain

    double d_radius;  // radius of solid phase particles

private:

    /*
     * Verify that configuration database from specified database.
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     *
     * Return value
     * ------------
     * true
     *
     * Notes
     * -----
     * - If 'config_db' is not valid, this method throws an exception
     *   containing error information.
     */
     void verifyConfigurationDatabase(
        const shared_ptr<tbox::Database>& config_db) const;

    /*
     * Load configuration parameters from specified database.
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     */
    void loadConfiguration(const shared_ptr<tbox::Database>& config_db);
};
