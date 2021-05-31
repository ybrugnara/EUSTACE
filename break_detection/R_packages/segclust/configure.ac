AC_INIT([segclust],0.08)

dnl Copy/Pasted from documentation example.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC
AC_PROG_CPP

dnl Package specific configuration switches
AC_PATH_PROG(GSL_CONFIG,gsl-config)
if test -z "${GSL_CONFIG}" ; then
   AC_MSG_ERROR("GSL library configuration tool 'gsl-config' not found.")
fi
AC_SUBST(GSL_CONFIG)

AC_CONFIG_FILES(src/Makevars)
AC_OUTPUT()

