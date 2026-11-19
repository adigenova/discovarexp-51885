# AC_OPENMP -- adapted from AX_OPENMP (Autoconf Archive), which carries
# the following notice:
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copying and distribution of this file, with or without modification,
#   are permitted in any medium without royalty provided the copyright
#   notice and this notice are preserved.
#
# This file provides the AC_OPENMP macro used by configure.ac (the tree
# historically shipped it only inside the pre-generated configure script,
# so it is recorded here for autoreconf).
#
# Sets OPENMP_CFLAGS and OPENMP_CXXFLAGS.

AC_DEFUN([AC_OPENMP], [
AC_PREREQ([2.59])
AC_ARG_WITH([openmp],
    [AS_HELP_STRING([--with-openmp],
        [use OpenMP (default: yes) if available])],
    [if test "$withval" = no; then
         acx_openmp=no
     else
         acx_openmp=yes
     fi],
    [acx_openmp=yes])
if test "$acx_openmp" = yes; then

AC_CACHE_CHECK([for OpenMP flag of C compiler], [ac_cv_prog_c_openmp], [
    ac_cv_prog_c_openmp=unknown
    # Try each of the likely candidates.
    AC_LANG_PUSH([C])
    for flag in -fopenmp -xopenmp -openmp -mp -omp -qsmp=omp; do
        ac_save_CFLAGS=$CFLAGS
        CFLAGS="$flag $CFLAGS"
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <omp.h>
            ]], [[
return omp_get_num_threads ();
            ]])],[ac_cv_prog_c_openmp=$flag; break],[])
        CFLAGS=$ac_save_CFLAGS
    done
    AC_LANG_POP([C])
])
if test "$ac_cv_prog_c_openmp" != unknown; then
    OPENMP_CFLAGS="$ac_cv_prog_c_openmp"
else
    OPENMP_CFLAGS=""
fi

AC_CACHE_CHECK([for OpenMP flag of C++ compiler], [ac_cv_prog_cxx_openmp], [
    ac_cv_prog_cxx_openmp=unknown
    AC_LANG_PUSH([C++])
    for flag in -fopenmp -xopenmp -openmp -mp -omp -qsmp=omp; do
        ac_save_CXXFLAGS=$CXXFLAGS
        CXXFLAGS="$flag $CXXFLAGS"
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <omp.h>
            ]], [[
return omp_get_num_threads ();
            ]])],[ac_cv_prog_cxx_openmp=$flag; break],[])
        CXXFLAGS=$ac_save_CXXFLAGS
    done
    AC_LANG_POP([C++])
])
if test "$ac_cv_prog_cxx_openmp" != unknown; then
    OPENMP_CXXFLAGS="$ac_cv_prog_cxx_openmp"
else
    OPENMP_CXXFLAGS=""
fi

else
    OPENMP_CFLAGS=""
    OPENMP_CXXFLAGS=""
fi
AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(OPENMP_CXXFLAGS)
])
