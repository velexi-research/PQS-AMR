/*! \file utilities.h
 *
 * \brief
 * Header file for pqs library utility functions and macros.
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

#ifndef INCLUDED_PQS_pqs_utilities_h
#define INCLUDED_PQS_pqs_utilities_h

// --- Headers, namespaces, and type declarations

// Standard library

// Namespaces

// Class/type declarations


// --- Utility functions and macros

/*!
 * Macro for handling errors that arise in the pqs library.
 *
 * Parameters
 * ----------
 * obj: pointer to object that the error originates in
 *
 * message: error message
 */
#define PQS_ERROR(obj, method, message) \
    throw "PQS_ERROR::" + string(typeid(obj).name()) + "::" + \
          string(method) + "::" + string(message)

#endif
