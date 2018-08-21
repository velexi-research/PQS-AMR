/*! \file macros.h
 *
 * \brief
 * Header file for PQS utility macros.
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

#ifndef INCLUDED_PQS_macros_h
#define INCLUDED_PQS_macros_h

// --- Headers, namespaces, and type declarations

// Standard library

// Namespaces

// Class/type declarations


// --- Utility macros

/*!
 * Macro for creating local int array copy of IntVector.
 *
 * Parameters
 * ----------
 * int_array: variable name for integer array
 *
 * int_vector: IntVector to copy
 */
#define PQS_INT_VECT_TO_INT_ARRAY(int_array, int_vector) \
    int int_array[3]; \
    for (int int_array_i=0; int_array_i < int_vector.getDim().getValue(); \
            int_array_i++) { \
        int_array[int_array_i] = int_vector[int_array_i]; \
    }

#endif  // INCLUDED_PQS_macros_h
