AC_INIT([ncdfFlow],1.5.10,wjiang2@fhcrc.org)


#: ${R_HOME=`R RHOME`}
#		if test -z "${R_HOME}"; then
#			echo "Could not determine R_HOME"
#			exit 1
#		fi
#		CC=`${R_HOME}/bin/R CMD config CC`
#		CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
#		CPPFAGS=`${R_HOME}/bin/R CMD config CPPFLAGS`
#		LDFLAGS=`${R_HOME}/bin/R CMD config LDFLAGS`
#		AC_PROG_CC


#	AC_MSG_NOTICE([using extenal hdf5 library])
	
	AC_ARG_WITH([hdf5],
	    [AC_HELP_STRING([--with-hdf5=DIR], [prefix of hdf5 libraries])],
	    [HDF5_PREFIX="${withval}"]
	)
	HDF5_FILE_TOTEST="libhdf5.so"
	#if prefix is set by user then search for libhdf5 under that path
	if test -n "${HDF5_PREFIX}"; then
			for HDF5_LIB_PATH in lib lib64
			do
			HDF5_LIB_PATH_ABS=${HDF5_PREFIX}/${HDF5_LIB_PATH}
			AC_MSG_NOTICE([search for hdf5 library in: ${HDF5_LIB_PATH_ABS}])
			if test -e ${HDF5_LIB_PATH_ABS}/${HDF5_LIB_TOTEST}; then
				HDF5_LIBS="-L${HDF5_LIB_PATH_ABS}"
				HDF5_LIBS="${HDF5_LIBS} -lhdf5 -lhdf5_hl"
				HDF5_CFLAGS="-I${HDF5_PREFIX}/include"
				AC_MSG_NOTICE([found hdf5 library in: ${HDF5_LIB_PATH_ABS}])
				break;
			fi
			done
	# if prefix not provided,then search some default standard location			
	else  
		for HDF5_PREFIX in "/usr" "/usr/local"
		do
			for HDF5_LIB_PATH in lib lib64
			do
			HDF5_LIB_PATH_ABS=${HDF5_PREFIX}/${HDF5_LIB_PATH}
			AC_MSG_NOTICE([search for hdf5 library in: ${HDF5_LIB_PATH_ABS}])
			if test -e ${HDF5_LIB_PATH_ABS}/${HDF5_LIB_TOTEST}; then
				HDF5_LIBS="-L${HDF5_LIB_PATH_ABS}"
				HDF5_LIBS="${HDF5_LIBS} -lhdf5 -lhdf5_hl"
				HDF5_CFLAGS="-I${HDF5_PREFIX}/include"
				AC_MSG_NOTICE([found hdf5 library in: ${HDF5_LIB_PATH_ABS}])
				break;
			fi
			done
		if test -n "${HDF5_LIBS}"; then
			break
		fi	
		done
	fi
	
	if test -z "${HDF5_LIBS}"; then
			AC_MSG_ERROR([hdf5 library not found in "${HDF5_PREFIX}"])
	fi
	

	

AC_MSG_NOTICE([using bundled netcdf4 library])
netcdf_basename=trimmed_netcdf-4.2.1.1
netcdf_tarname="${netcdf_basename}.tar.gz"
cd src
if test -d ${netcdf_basename};then 
	AC_MSG_NOTICE(['found ' ${netcdf_basename} ' header sources and tar archive;using what is there.'])
else 
	AC_MSG_NOTICE([untarring netcdf include tree...])
	gunzip -dc ${netcdf_tarname} | tar xf - 
fi
cd ../

PKG_CPPFLAGS=${HDF5_CFLAGS}
PKG_LIBS="${HDF5_LIBS}"

AC_MSG_NOTICE([Using the following compilation and linking flags for ncdfFlow])
AC_MSG_NOTICE([	PKG_CPPFLAGS=${PKG_CPPFLAGS}])
AC_MSG_NOTICE([	PKG_LIBS=${PKG_LIBS}])

AC_SUBST(netcdf_basename)
AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)


AC_OUTPUT(src/Makevars)
