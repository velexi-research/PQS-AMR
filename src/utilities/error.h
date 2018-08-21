/*! \file error.h
 *
 * \brief
 * Header file for error handling in PQS library.
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

#ifndef INCLUDED_PQS_error_h
#define INCLUDED_PQS_error_h

// --- Headers, namespaces, and type declarations

// Standard
#include <string>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

// PQS
#include "PQS/PQS_config.h"

// Namespaces
using namespace SAMRAI;

// Class/type declarations

// --- Macros

/*!
 * Macro for gracefully aborting PQS programs.
 *
 * Parameters
 * ----------
 * message: error message
 */
#define PQS_ABORT(message) \
    tbox::pout << "Program abort called..." << std::endl; \
    tbox::pout << "Error: " << std::string(message) << std::endl; \
    tbox::pout << std::flush; \
    tbox::SAMRAI_MPI::abort();

/*!
 * Macro for raising errors that arise in the PQS library.
 *
 * Parameters
 * ----------
 * obj: pointer to object that the error originates in
 *
 * method: name of method that error originates in
 *
 * message: error message
 */
#define PQS_ERROR(obj, method, message) \
    throw PQS::exception(std::string(typeid(obj).name()) + "::" + \
                         std::string(method) + "::" + std::string(message));

/*!
 * Macro for raising errors in static methods that arise in the PQS library.
 *
 * Parameters
 * ----------
 * class_name: name of class
 *
 * method: name of method that error originates in
 *
 * message: error message
 */
#define PQS_ERROR_STATIC(class_name, method, message) \
    throw PQS::exception(std::string(class_name) + "::" + \
                         std::string(method) + "::" + std::string(message));


// --- PQS::exception Class

namespace PQS {

/*!
 * Custom exception class for PQS library.
 */
class exception : std::exception {

public:
    exception(const std::string message) : std::exception() {
        this->message = message;
    }

    virtual const char* what() const throw() {
        return message.c_str();
    }

private:
    std::string message;
};

}  // PQS namespace

#endif  // INCLUDED_PQS_error_h
