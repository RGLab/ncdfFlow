AC_INIT([ncdfFlow],1.5.10,wjiang2@fhcrc.org)


	AC_ARG_WITH([hdf5],
	    [AC_HELP_STRING([--with-hdf5=DIR], [prefix of hdf5 libraries])],
	    [HDF5_PREFIX="${withval}"]
	)
	HDF5_LIB_TOTEST="libhdf5.so"
	#if prefix is set by user then search for libhdf5 under that path
	if test -n "${HDF5_PREFIX}"; then
			for HDF5_LIB_PATH in lib lib64
			do
			HDF5_LIB_PATH_ABS=${HDF5_PREFIX}/${HDF5_LIB_PATH}
			AC_MSG_NOTICE([search for hdf5 library in: ${HDF5_LIB_PATH_ABS}])
			if test -e ${HDF5_LIB_PATH_ABS}/${HDF5_LIB_TOTEST}; then
				HDF5_LIBS="-L${HDF5_LIB_PATH_ABS}"
				HDF5_LIBS="${HDF5_LIBS} -lhdf5"
				HDF5_CFLAGS="-I${HDF5_PREFIX}/include"
				AC_MSG_NOTICE([found hdf5 library in: ${HDF5_LIB_PATH_ABS}])
				break;
			fi
			done
				
	else  
		# if prefix not provided, try pkg-config
		AC_MSG_NOTICE([No directory was specified for --with-hdf5. Trying to find hdf5 using pkg-config.])
		if test -z "${PKG_CONFIG}"; then
			AC_PATH_PROG(PKG_CONFIG, pkg-config)
	    fi
	    if ! test -z "${PKG_CONFIG}" ; then
            HDF_CONFIG="${PKG_CONFIG} hdf5"
            HDF5_CFLAGS="`${HDF_CONFIG} --cflags`"
		fi
		
		if test -z "${HDF5_CFLAGS}" ; then
        	AC_MSG_NOTICE([pkg-config was not able to find the hdf5 library.])
        	AC_MSG_NOTICE([Trying some common locations."])
			#still not detected, then search some default standard location
			for HDF5_PREFIX in "/usr/local" #skip searching/usr since it is standard system path
			do
				for HDF5_LIB_PATH in lib lib64 lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
				do
				HDF5_LIB_PATH_ABS=${HDF5_PREFIX}/${HDF5_LIB_PATH}
				AC_MSG_NOTICE([search for hdf5 library in: ${HDF5_LIB_PATH_ABS}])
				if test -e ${HDF5_LIB_PATH_ABS}/${HDF5_LIB_TOTEST}; then
					HDF5_LIBS="-L${HDF5_LIB_PATH_ABS}"
					HDF5_LIBS="${HDF5_LIBS} -lhdf5"
					HDF5_CFLAGS="-I${HDF5_PREFIX}/include"
					AC_MSG_NOTICE([found hdf5 library in: ${HDF5_LIB_PATH_ABS}])
					break;
				fi
				done
			if test -n "${HDF5_LIBS}"; then
				break
			fi	
			done
		else
			HDF5_LIBS="`${HDF_CONFIG} --libs`"
		fi
	fi
	
	#fail to detect hdf5 lib, use the standard system path (which can be ommited) 
	if test -z "${HDF5_LIBS}"; then
			HDF5_LIBS="-lhdf5"
			HDF5_CFLAGS="" 
	fi

# add mpi header path when hdf5 is built with parallel enabled
# otherwise the H5public.h will fail to find mpi.h	
if test -z "${H5CC}"; then
			AC_PATH_PROG(H5CC, h5cc)
fi
if ! test -z "${H5CC}" ; then
	if test -n "`h5cc -showconfig|grep 'Parallel HDF5'|grep yes`"; then
		MPI_CPPFLAGS="`mpicc -showme:compile`"
	fi
fi
	


PKG_CPPFLAGS="${HDF5_CFLAGS} ${MPI_CPPFLAGS}"
PKG_LIBS="${HDF5_LIBS}"

AC_MSG_NOTICE([PKG_CPPFLAGS = ${PKG_CPPFLAGS}])
AC_MSG_NOTICE([PKG_LIBS = ${PKG_LIBS}])

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)


AC_OUTPUT(src/Makevars)
