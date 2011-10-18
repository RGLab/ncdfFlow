#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.67)
AC_INIT([ncdfFlow],0.10.1)
AC_CONFIG_SRCDIR([src/])

 : ${R_HOME=`R RHOME`}
     if test -z "${R_HOME}"; then
       echo "could not determine R_HOME"
       exit 1
     fi
     CC=`"${R_HOME}/bin/R" CMD config CC`
     CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
     CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
     LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

# Checks for programs.
AC_PROG_CC

AC_ARG_WITH([netcdf],[AS_HELP_STRING([--with-netcdf], [Prefix to the location of the netcdf installation])],[

if test -e $with_netcdf/include -a -e $with_netcdf/lib -a -d $with_netcdf/include -a -d $with_netcdf/lib ; then
	AC_MSG_RESULT($with_netcdf)
		CFLAGS="-I$with_netcdf/include $CFLAGS"
	CPPFLAGS="-I$with_netcdf/include $CPPFLAGS"
	LDFLAGS="-L$with_netcdf/lib $LDFLAGS"
	LIBS="-lnetcdf"
else
AC_CHECK_LIB([netcdf],[nc_def_var_deflate],[],AC_MSG_ERROR([netCDF library not found!]))
fi

],[
AC_CHECK_LIB([netcdf],[nc_def_var_deflate],[],AC_MSG_ERROR([netCDF library not found!]))
])
# Checks for libraries.
#AC_CHECK_LIB([netcdf],[nc_def_var_deflate],[],AC_MSG_ERROR([netCDF library not found!]))

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_HEADERS([netcdf.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_HEADER_STDBOOL

# Checks for library functions.
AC_SUBST(LIBS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
