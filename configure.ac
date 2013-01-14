AC_INIT([ncdfFlow],1.5.10,wjiang2@fhcrc.org)
AC_ARG_WITH([bundled-netcdf],
    [AC_HELP_STRING([--with-bundled-netcdf=yes/no]
    				, [whether to build ncdfFlow with the bundled netcdf library]
    				)
    ],
    [NC_BOUNDLED="${withval}"],
    [NC_BOUNDLED="yes"]
)

# library path for HDF5
AC_ARG_WITH([hdf5],
    [AC_HELP_STRING([--with-hdf5=DIR], [location of hdf5 libraries])],
    [hdf5_path="${withval}"]
)
if test -n "$hdf5_path"; then
	HDF5_LIBS="-L${hdf5_path}/lib64"
else
	HDF5_LIBS="-L/usr/lib64 -L/usr/local/lib64"
fi
HDF5_LIBS="${HDF5_LIBS} -lhdf5 -lhdf5_hl"

: ${R_HOME=`R RHOME`}
		if test -z "${R_HOME}"; then
			echo "Could not determine R_HOME"
			exit 1
		fi
		CC=`${R_HOME}/bin/R CMD config CC`
		CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
		CPPFAGS=`${R_HOME}/bin/R CMD config CPPFLAGS`
		LDFLAGS=`${R_HOME}/bin/R CMD config LDFLAGS`
		AC_PROG_CC

	
if test ${NC_BOUNDLED} = "yes"; then

	PKG_DIR=`pwd`
	cd src
	BASENETCDFNAME="trimmed_netcdf-4.2.1.1"
	NETCDFTGZNAME="${BASENETCDFNAME}.tar.gz"
	if test -d ./${BASENETCDFNAME};
	        then
	                echo 'found ' $BASENETCDFNAME ' header sources and tar archive;using what is there.'
	        else
	                echo "untarring netcdf include tree...";
	                gunzip -dc ${NETCDFTGZNAME} | tar xf -;
	fi;

    AC_MSG_NOTICE([compiling netcdf library...])
    
    BUILD_DIR="${BASENETCDFNAME}/build"
    mkdir "${BUILD_DIR}"
    cd "${BUILD_DIR}"
    CPICFLAGS=`${R_HOME}/bin/R CMD config CPICFLAGS`
    `${CC} ${CFLAGS} ${CPICFLAGS} -c -I../ -I../include ${CPPFLAGS} ../liblib/*.c ../libdispatch/*.c ../libsrc/*.c ../libsrc4/*.c`
    cd "${PKG_DIR}"
	NETCDF_LIBS="${BUILD_DIR}/*.o"

else

	AC_ARG_WITH([nc-config],
	    [AC_HELP_STRING([--with-nc-config=PATH], [location of nc-config])],
	    [NC_CONFIG=$withval])
	if test -z "${NC_CONFIG}"; then
	   AC_PATH_PROG(NC_CONFIG, nc-config)
	else
	  if ! test -x "${NC_CONFIG}"; then
	    NC_CONFIG=""
	    AC_MSG_WARN(['NC_CONFIG' is set to an invalid value and will be ignored])
	  fi
	fi

	if test -n "${NC_CONFIG}"; then
	
	  NETCDF_CPPFLAGS=`${NC_CONFIG} --cflags`
	  NETCDF_LIBS=`${NC_CONFIG} --libs`
	  NC4=`${NC_CONFIG} --has-nc4`
		if  test $NC4 != "yes"; then
			echo "Need netcdf version 4 linked against hdf5"
			exit 1
		fi
	else

		AC_ARG_WITH([netcdf-include],
		    [AC_HELP_STRING([--with-netcdf-include=DIR], 
				   [location of netcdf header files])],
		    [netcdf_include_path=$withval])
		if test -n "${netcdf_include_path}"; then
		   NETCDF_CPPFLAGS="-I${netcdf_include_path}"
		else
			if test -n ${NETCDF_INCLUDE}; then
			   NETCDF_CPPFLAGS="-I${NETCDF_INCLUDE}"
			fi
		fi
		AC_CHECK_HEADERS(netcdf.h, , AC_MSG_ERROR(netcdf header netcdf.h not found))
	
		AC_ARG_WITH([netcdf-lib],
		    [AC_HELP_STRING([--with-netcdf-lib=DIR], [location of netcdf libraries])],
		    [netcdf_lib_path=$withval])
		if test -n "${netcdf_lib_path}"; then
		   NETCDF_LIBS="-L${netcdf_lib_path}"
		else
			if test -n "${NETLIB_LIB}"; then
		   		NETCDF_LIBS="-L${NETLIB_LIB}"
		   	fi
		fi
	fi
fi

PKG_CPPFLAGS="${NETCDF_CPPFLAGS}"
PKG_LIBS="${NETCDF_LIBS} ${HDF5_LIBS} ${LDFLAGS}"

AC_MSG_NOTICE([Using the following compilation and linking flags for ncdfFlow])
AC_MSG_NOTICE([	PKG_CPPFLAGS=${PKG_CPPFLAGS}])
AC_MSG_NOTICE([	PKG_LIBS=${PKG_LIBS}])
AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)

AC_OUTPUT(src/Makevars)
