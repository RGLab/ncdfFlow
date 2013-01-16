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
	    [AC_HELP_STRING([--with-hdf5=DIR], [location of hdf5 libraries])],
	    [hdf5_path="${withval}"]
	)
	if test -n "$hdf5_path"; then
		HDF5_LIBS="-L${hdf5_path}/lib64"
		AC_MSG_NOTICE([hdf5 library path:${HDF5_LIBS}])
	else
		HDF5_LIBS="-L/usr/lib64 -L/usr/local/lib64"
		AC_MSG_NOTICE([trying default hdf5 library path:${HDF5_LIBS}])
	fi
	HDF5_LIBS="${HDF5_LIBS} -lhdf5 -lhdf5_hl"

	

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

PKG_CPPFLAGS=
PKG_LIBS="${HDF5_LIBS}"

AC_MSG_NOTICE([Using the following compilation and linking flags for ncdfFlow])
AC_MSG_NOTICE([	PKG_CPPFLAGS=${PKG_CPPFLAGS}])
AC_MSG_NOTICE([	PKG_LIBS=${PKG_LIBS}])

AC_SUBST(netcdf_basename)
AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)


AC_OUTPUT(src/Makevars)
