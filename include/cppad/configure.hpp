# ifndef CPPAD_CONFIGURE_HPP
# define CPPAD_CONFIGURE_HPP
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
// SPDX-FileCopyrightText: Bradley M. Bell <bradbell@seanet.com>
// SPDX-FileContributor: 2003-24 Bradley M. Bell
// ----------------------------------------------------------------------------

/*!
{xrst_begin configure.hpp dev}
{xrst_spell
   cl
   cmd
   complier
   gettimeofday
   mkstemp
   noexcept
   nullptr
   pragmas
   tmpnam
   unreferenced
   yyyy
   yyyymmdd
}

Preprocessor Symbols Set By CMake Command
#########################################

CPPAD_COMPILER_HAS_CONVERSION_WARN
**********************************
is the compiler a variant of g++ and has conversion warnings
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_COMPILER_HAS_CONVERSION_WARN 1
/* {xrst_code}
{xrst_spell_on}

CPPAD_DISABLE_SOME_MICROSOFT_COMPILER_WARNINGS
**********************************************
This macro is only used to document the pragmas that disables the
follow warnings:

C4100
=====
unreferenced formal parameter.

C4127
=====
conditional expression is constant.

C4723
=====
The second operand in a divide operation evaluated to zero at compile time.

{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_DISABLE_SOME_MICROSOFT_COMPILER_WARNINGS 1
# if _MSC_VER
# pragma warning( disable : 4100 )
# pragma warning( disable : 4127 )
# pragma warning( disable : 4723 )
# endif
# undef CPPAD_DISABLE_SOME_MICROSOFT_COMPILER_WARNINGS
/* {xrst_code}
{xrst_spell_on}

CPPAD_DEBUG_AND_RELEASE
***********************
Starting with 2023-12-24,
this flag is set by the cmake command; see
:ref:`cmake@cppad_debug_and_release` .
Before then, one would add -D CPPAD_DEBUG_AND_RELEASE
when compiling CppAD code.
{xrst_code hpp} */
# define CPPAD_DEBUG_AND_RELEASE 1
/* {xrst_code}

CPPAD_USE_CPLUSPLUS_2011
************************
Deprecated 2020-12-03:
Is it OK to use C++11 features. This is always 1 (for true).
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_USE_CPLUSPLUS_2011 1
/* {xrst_code}
{xrst_spell_on}

CPPAD_USE_CPLUSPLUS_2017
************************
Deprecated 2020-12-03:
Is it OK for CppAD use C++17 features.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_USE_CPLUSPLUS_2017 0
/* {xrst_code}
{xrst_spell_on}

CPPAD_PACKAGE_STRING
********************
cppad-yyyymmdd as a C string where yyyy is year, mm is month, and dd is day.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_PACKAGE_STRING "cppad-20241211"
/* {xrst_code}
{xrst_spell_on}

CPPAD_HAS_ADOLC
***************
Was include_adolc=true on the cmake command line.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_HAS_ADOLC 0
/* {xrst_code}
{xrst_spell_on}

CPPAD_HAS_COLPACK
*****************
Was a colpack_prefix specified on the cmake command line.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_HAS_COLPACK 0
/* {xrst_code}
{xrst_spell_on}

CPPAD_HAS_EIGEN
***************
Was Eigen found and c++14 is supported.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_HAS_EIGEN 1
/* {xrst_code}
{xrst_spell_on}

CPPAD_HAS_IPOPT
***************
Was include_ipopt=true on the cmake command line.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_HAS_IPOPT 0
/* {xrst_code}
{xrst_spell_on}

CPPAD_DEPRECATED
****************
This symbol is not currently being used.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_DEPRECATED 
/* {xrst_code}
{xrst_spell_on}

CPPAD_BOOSTVECTOR
*****************
If this symbol is one, and _MSC_VER is not defined,
we are using boost vector for CPPAD_TESTVECTOR.
It this symbol is zero,
we are not using boost vector for CPPAD_TESTVECTOR.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_BOOSTVECTOR 0
/* {xrst_code}
{xrst_spell_on}

CPPAD_CPPADVECTOR
*****************
If this symbol is one,
we are using CppAD vector for CPPAD_TESTVECTOR.
It this symbol is zero,
we are not using CppAD vector for CPPAD_TESTVECTOR.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_CPPADVECTOR 1
/* {xrst_code}
{xrst_spell_on}

CPPAD_STDVECTOR
***************
If this symbol is one,
we are using standard vector for CPPAD_TESTVECTOR.
It this symbol is zero,
we are not using standard vector for CPPAD_TESTVECTOR.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_STDVECTOR 0
/* {xrst_code}
{xrst_spell_on}

CPPAD_EIGENVECTOR
*****************
If this symbol is one,
we are using Eigen vector for CPPAD_TESTVECTOR.
If this symbol is zero,
we are not using Eigen vector for CPPAD_TESTVECTOR.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_EIGENVECTOR 0
/* {xrst_code}
{xrst_spell_on}

CPPAD_HAS_GETTIMEOFDAY
**********************
If this symbol is one, and _MSC_VER is not defined,
this system supports the gettimeofday function.
Otherwise, this symbol should be zero.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_HAS_GETTIMEOFDAY 1
/* {xrst_code}
{xrst_spell_on}

CPPAD_TAPE_ADDR_TYPE
********************
Is the type used to store address on the tape.
If it is not size_t, then
{xrst_code cpp}
   sizeof(CPPAD_TAPE_ADDR_TYPE) < sizeof( size_t )
{xrst_code}
can be used to conserve memory.
This type must support std::numeric_limits,
the <= operator,
and conversion to size_t.
Make sure that the type chosen returns true for is_pod<CPPAD_TAPE_ADDR_TYPE>
in pod_vector.hpp.
This type is later defined as addr_t in the CppAD namespace.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_TAPE_ADDR_TYPE unsigned int
/* {xrst_code}
{xrst_spell_on}

CPPAD_IS_SAME_TAPE_ADDR_TYPE_SIZE_T
***********************************
Is size_t the type the same as CPPAD_TAPE_ADDR_TYPE.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_IS_SAME_TAPE_ADDR_TYPE_SIZE_T \
   0
/* {xrst_code}
{xrst_spell_off}


CPPAD_TAPE_ID_TYPE
******************
Is the type used to store tape identifiers.
If it is not size_t, then
{xrst_code cpp}
   sizeof(CPPAD_TAPE_ID_TYPE) < sizeof( size_t )
{xrst_code}
can be used to conserve memory.
This type must support std::numeric_limits,
the <= operator,
and conversion to size_t.
Make sure that the type chosen returns true for is_pod<CPPAD_TAPE_ID_TYPE>
in pod_vector.hpp.
This type is later defined as tape_id_t in the CppAD namespace.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_TAPE_ID_TYPE unsigned int
/* {xrst_code}
{xrst_spell_on}

CPPAD_MAX_NUM_THREADS
*********************
Specifies the maximum number of threads that CppAD can support
(must be greater than or equal four).

The user may define CPPAD_MAX_NUM_THREADS before including any of the CppAD
header files.  If it is not yet defined,
{xrst_spell_off}
{xrst_code hpp} */
# ifndef CPPAD_MAX_NUM_THREADS
# define CPPAD_MAX_NUM_THREADS 48
# endif
/* {xrst_code}
{xrst_spell_on}

CPPAD_HAS_MKSTEMP
*****************
if true, mkstemp works in C++ on this system.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_HAS_MKSTEMP 1
/* {xrst_code}
{xrst_spell_on}

CPPAD_HAS_TMPNAM_S
******************
If true, tmpnam_s works in C++ on this system.
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_HAS_TMPNAM_S 0
/* {xrst_code}
{xrst_spell_on}

CPPAD_NULL
**********
Deprecated 2020-12-03:
This preprocessor symbol was used for a null pointer before c++11.
Replace it by ``nullptr`` .

CPPAD_NOEXCEPT
**************
Deprecated 2020-12-03:
This preprocessor symbol was used for no exception before c++11,
replace it by ``noexcept`` .

CPPAD_NDEBUG_NOEXCEPT
=====================
This preprocessor symbol is
``noexcept`` when ``NDEBUG`` is defined.
Otherwise it is empty.

CPPAD_C_COMPILER_CMD
********************
This is the command that runs the C compiler as a C string;
i.e., surrounded by double quotes.
It can be used to run the C compiler; e.g. see for :ref:`create_dll_lib-name` .
{xrst_code hpp} */
# define CPPAD_C_COMPILER_CMD "cc"
/* {xrst_code}

CPPAD_C_COMPILER_GNU_FLAGS
**************************
If true, the C complier uses the same flags as ``gcc``
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_C_COMPILER_GNU_FLAGS 1
/* {xrst_code}
{xrst_spell_on}

CPPAD_C_COMPILER_MSVC_FLAGS
***************************
If true, the C complier uses the same flags as ``cl``
{xrst_spell_off}
{xrst_code hpp} */
# define CPPAD_C_COMPILER_MSVC_FLAGS 0
/* {xrst_code}
{xrst_spell_on}

CPPAD_IS_SAME_UNSIGNED_INT_SIZE_T
*********************************
If true, ``unsigned int`` and ``size_t`` are the same type
{xrst_code hpp} */
# define CPPAD_IS_SAME_UNSIGNED_INT_SIZE_T 0
/* {xrst_code}

CPPAD_PADDING_BLOCK_T
*********************
Is a string used to define an object that pads the block_t structure
so that its size is a multiple of the size of a double.
{xrst_code hpp} */
# define CPPAD_PADDING_BLOCK_T 
/* {xrst_code}

{xrst_end configure.hpp}
*/
// -------------------------------------------------
# define CPPAD_NULL                nullptr
# define CPPAD_NOEXCEPT            noexcept
//
# ifdef NDEBUG
# define CPPAD_NDEBUG_NOEXCEPT     noexcept
# else
# define CPPAD_NDEBUG_NOEXCEPT
# endif
// -------------------------------------------------

# endif
