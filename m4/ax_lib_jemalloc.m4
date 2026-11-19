############################################################
#
# JEMALLOC path
#
# This function finds the jemalloc library.  If a --with-jemalloc=DIR
# is supplied it is used directly; otherwise a system jemalloc is
# detected.  If none is found the build proceeds without it: jemalloc
# is an optional dependency (recommended for its fragmentation
# behavior on large genomes).
#
# input: none
#
# output: JEMALLOC_LIB
#
############################################################
AC_DEFUN([AC_JEMALLOC_LIB],
 [
  JEMALLOC_LIB=""

  AC_ARG_WITH([jemalloc],
        AS_HELP_STRING([--with-jemalloc=DIR],
                [Force directory for location of jemalloc library.]),
        [
        if test -d $withval
        then
            JEMALLOC_LIB="-L$withval -Wl,-rpath,$withval -ljemalloc"
        else
            AC_MSG_ERROR(--with-jemalloc expected directory name)
        fi
        ]
  )

  if test "x$JEMALLOC_LIB" = "x" ; then
      AC_CHECK_LIB([jemalloc],[malloc_stats_print],
          [JEMALLOC_LIB="-ljemalloc"],
          [AC_MSG_NOTICE([jemalloc not found: building without it.  For best memory behavior on large genomes, install jemalloc or use --with-jemalloc=DIR.])])
  fi

  AC_SUBST(JEMALLOC_LIB,[$JEMALLOC_LIB])
 ]
)
