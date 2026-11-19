AC_INIT([DiscovarExp], [51885], [crdhelp@broadinstitute.org])
AM_INIT_AUTOMAKE
LT_INIT([disable-shared])
AC_PROG_CXX
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])

AC_LITTLE_ENDIAN_CHECK
AC_GXX_CHECK
AX_PTHREAD
AC_JEMALLOC_LIB
AC_OPENMP
AC_OPENMP_CHECK
AC_LIBZ_LIB

dnl Modern-compiler build (tested with GCC 11 through 16).
dnl -std=gnu++14: the tree still uses a few pre-C++17 idioms
dnl (binary_function, random_shuffle, dynamic exception specs) that are
dnl removed in C++17; under gnu++14 they are merely deprecated.
dnl Flags dropped versus the 2014 (GCC 4.8) build: -fno-nonansi-builtins
dnl (no longer accepted by GCC), -mieee-fp (x86-only, rejected elsewhere),
dnl -Wno-long-long (moot since C++11).
CXXFLAGS="$CXXFLAGS -std=gnu++14 -Wextra -Wall -Wsign-promo -Woverloaded-virtual -Wendif-labels -Wno-unused -Wno-deprecated -Wno-parentheses -Wno-unused-parameter -fno-strict-aliasing -iquote . -ggdb -DNDEBUG"

OS_RELEASE="`uname -r`"
AC_SUBST(CXXFLAGS)
AC_SUBST(OS_RELEASE, [$OS_RELEASE])
AC_SUBST(VERSION, [$VERSION])
AC_OUTPUT
