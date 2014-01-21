AC_INIT([ncdfFlow],1.5.10,wjiang2@fhcrc.org)


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
				HDF5_LIBS="${HDF5_LIBS} -lhdf5"
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
	fi
	
	if test -z "${HDF5_LIBS}"; then
			AC_MSG_ERROR([hdf5 library not found in "${HDF5_PREFIX}"])
	fi
	

	


PKG_CPPFLAGS=${HDF5_CFLAGS}
PKG_LIBS="${HDF5_LIBS}"


AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)


AC_OUTPUT(src/Makevars)
