/*! \file InterfaceInitModule.h
 *
 * \brief
 * Header file for concrete subclass of pqs::InterfaceInitStrategy to use in
 * example application that has no solid phase.
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
#include <iosfwd>
#include <memory>

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
namespace SAMRAI { namespace tbox { class Database; } }


// --- Class declaration

class InterfaceInitModule: public pqs::InterfaceInitStrategy
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
     */
    InterfaceInitModule(
            const shared_ptr<tbox::Database>& config_db,
            const int dim);

    /*!
     * Default destructor frees memory allocated for the object.
     */
    virtual ~InterfaceInitModule();

    //! @}

    //! {@

    /*!
     ************************************************************************
     *
     * @name Interface initialization methods
     *
     ************************************************************************/

    /*
     * Set phi so that fluid-fluid interface is a perturbed circle/sphere
     * that is fully contained in computational domain.
     */
    virtual void initializeInterface(hier::Patch& patch, int phi_id);

    //! @}

protected:

    /****************************************************************
     *
     * Data Members
     *
     ****************************************************************/

    // --- Parameters

    int d_dim;  // problem dimension

    string d_interface_type;
    double *d_center;
    double d_radius;

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
