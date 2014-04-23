#include "wrappers.h"

void ERR(int e){
	char * errmsg ="";
//	error("Error from C API:");
//	errmsg = HEstring(e);
	REprintf("hdf Error: %s\n", errmsg);
}
/*
 * metaSize and compress are currently not used
 */
SEXP createFile(SEXP _fileName, SEXP _nEvent, SEXP _nChannel, SEXP _nSample, SEXP _dim) {

	SEXP k = allocVector(LGLSXP,1); //create logical scalar for return value
    int nDim = INTEGER(_dim)[0];

    int retval;
    if(nDim == 3)
    {
    	int nSample = INTEGER(_nSample)[0];
    	int nChannel = INTEGER(_nChannel)[0];
		int nEvent = INTEGER(_nEvent)[0];
		retval = _createFile3d(translateChar(STRING_ELT(_fileName, 0)), nSample, nChannel, nEvent);
    }
    else
    {

    	retval = _createFile2d(translateChar(STRING_ELT(_fileName, 0)));
    }

    if(retval < 0)
    {
    	ERR(retval);
    	LOGICAL(k)[0] = FALSE; //set return value as FALSE
    }

    LOGICAL(k)[0] = TRUE; //set return value as TRUE
    return(k);
}

/*
 * inline _writeSlice and _writeSlice2d code
 */
SEXP writeSlice(SEXP _fileName, SEXP _mat, SEXP _chIndx, SEXP _sampleIndx, SEXP _sampleName) {

	SEXP k = allocVector(LGLSXP,1);//create logical scalar for return value
	const char * fName = translateChar(STRING_ELT(_fileName, 0));
	double *mat = REAL(_mat);
	int * chIndx = INTEGER(_chIndx);
	int chCount = length(_chIndx);

    int nRow, nCol, sample;
    SEXP Rdim = getAttrib(_mat, R_DimSymbol);
    nRow = INTEGER(Rdim)[0];
    nCol = INTEGER(Rdim)[1];
    int nEvents = nRow;
	/*
	 * Open the file and the dataset.
	 */
	herr_t      status;
	hid_t  file, dataset,dataspace, memspace;         /* handles */
	file = H5Fopen(fName, H5F_ACC_RDWR, H5P_DEFAULT);//open file
	status = H5Lexists(file, DATASETNAME3d, H5P_DEFAULT);

	int is3d;
	dataset = -1;

	if(status == TRUE){
		dataset = H5Dopen2(file, DATASETNAME3d, H5P_DEFAULT);
		dataspace = H5Dget_space(dataset);    /* dataspace handle */
		int nDim = H5Sget_simple_extent_ndims(dataspace);
		is3d = nDim == 3;
	}
	else
		is3d = 0;

	if(is3d)
	{

		int sampleIndx = INTEGER(_sampleIndx)[0];

		hid_t attrID;
		/*
		 * Define the memory dataspace.
		 */
		hsize_t 	dimsm[2]; //dimenstions
		dimsm[0] = chCount;
		dimsm[1] = nEvents;
		memspace = H5Screate_simple(2,dimsm,NULL);


		/*
		 * Define hyperslab in the dataset.
		 */
		hsize_t      count[3];              /* size of the hyperslab in the file */
		hsize_t      offset[3];             /* hyperslab offset in the file */
		hsize_t      count_in[2];          /* size of the hyperslab in memory */
		hsize_t      offset_in[2];         /* hyperslab offset in memory */

		/*
		 * write subsets
		 */
		unsigned i;
		sampleIndx = sampleIndx -1;//convert from R to C indexing
		for(i = 0; i < chCount; i++){
			int colStart = chIndx[i] -1; //convert from R to C indexing
			offset[0] = sampleIndx;//start from sampleIndx-th sample
			offset[1] = colStart; //start from colStart-th channel
			offset[2] = 0; //start from the first event

			count[0]  = 1;//get one sample
			count[1]  = 1;//get one channel
			count[2]  = nEvents; //get all events


			status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,
												count, NULL);
			/*
			 * Define memory hyperslab.
			 */
			offset_in[0] = i;//start from ith column
			offset_in[1] = 0;//start from 0th event


			count_in[0]  = 1;//one channel
			count_in[1]  = nEvents; //all events

			status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_in, NULL,
					count_in, NULL);
			/*
			 * write data to hyperslab in the file from memory .
			 */
			status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, mat);

		}


		 /*
		  * get eCount attribute
		   */
		  hsize_t dims[3];
		  status  = H5Sget_simple_extent_dims(dataspace, dims, NULL); //get dimensions of datset
		  unsigned nSample = dims[0];//get total number of samples
		  if(sampleIndx >= nSample)
		  		error("writeSlice error!sample index exceeds the boundary.");
		  unsigned * eCount = (unsigned *) malloc(sizeof(unsigned) * nSample);
		  attrID = H5Aopen(dataset, "eventCount", H5P_DEFAULT);
		  status = H5Aread(attrID, H5T_NATIVE_UINT32, eCount);
		  //update the eCount for current sample
		  eCount[sampleIndx] = nEvents;
		  /*
		   * write back to hdf
		   */
		  status = H5Awrite(attrID, H5T_NATIVE_UINT32, eCount);

		  free(eCount);
		  H5Aclose(attrID);
	}
	else
	{
		const char * sampleName = translateChar(STRING_ELT(_sampleName, 0));
		/*
		 * Open the file and the dataset.
		 */
		if(dataset>0)
		{
			//close it if it was previously opened for checking dimension
			H5Dclose(dataset);
			H5Sclose(dataspace);
		}
		/*
		 * check if dataset already exists
		 */
		status = H5Lexists(file, sampleName, H5P_DEFAULT);
		if(status == FALSE)
		{
			/* Create the data space for the 2d mat. */
			hsize_t dims[2];
			dims[0] = chCount;
			dims[1] = nEvents;
			dataspace = H5Screate_simple(2, dims, NULL);

			hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);

			//set it to use chunking
			hsize_t		chunk_dims[2] = {1, nEvents};
			H5Pset_chunk(dcpl_id, 2, chunk_dims);

			/* Create the 2d mat. */
			dataset = H5Dcreate2(file, sampleName, H5T_IEEE_F32LE_g, dataspace,
								  H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
			H5Dclose(dcpl_id);
		}
		else
		{
			dataset = H5Dopen2(file, sampleName, H5P_DEFAULT);
			dataspace = H5Dget_space(dataset);    /* dataspace handle */
		}

		/*
		 * Define the memory dataspace.
		 */
		hsize_t 	dimsm[2]; //dimenstions
		dimsm[0] = chCount;
		dimsm[1] = nEvents;
		memspace = H5Screate_simple(2,dimsm,NULL);

		/*
		 * Define hyperslab in the dataset.
		 */
		hsize_t      count[2];              /* size of the hyperslab in the file */
		hsize_t      offset[2];             /* hyperslab offset in the file */
		hsize_t      count_in[2];          /* size of the hyperslab in memory */
		hsize_t      offset_in[2];         /* hyperslab offset in memory */

		/*
		 * write subsets
		 */
		unsigned i;

		for(i = 0; i < chCount; i++){
			int colStart = chIndx[i] -1; //convert from R to C indexing
			offset[0] = colStart; //start from colStart-th channel
			offset[1] = 0; //start from the first event

			count[0]  = 1;//get one channel
			count[1]  = nEvents; //get all events


			status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,
												count, NULL);
			/*
			 * Define memory hyperslab.
			 */
			offset_in[0] = i;//start from ith column
			offset_in[1] = 0;//start from 0th event


			count_in[0]  = 1;//one channel
			count_in[1]  = nEvents; //all events

			status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_in, NULL,
					count_in, NULL);
			/*
			 * write data to hyperslab in the file from memory .
			 */
			status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, mat);

		}

	}

	/*
	 * Close/release resources.
	 */
	H5Dclose(dataset);

	H5Sclose(dataspace);
	H5Sclose(memspace);

	H5Fclose(file);



	if(status < 0)
	{
		ERR(status);
		LOGICAL(k)[0] = FALSE; //set return value as FALSE
	}

	LOGICAL(k)[0] = TRUE; //set return value as TRUE
    return(k);

}
/*
 * The reason we have to inline the code from _readSlice is two fold:
 * 1. in order to use the memory buffer dynamically allocated by R, which is transient within the call
 *    and its size if decided by the events number stored in hdf
 * 2. we want to query the hdf format in order to dispatch to the different logic of IO
 * Both requires the hdf query, which causes disk IO, thus we do it here directly so that hdf only needs to be opened once
 */
SEXP readSlice(SEXP _fileName, SEXP _chIndx, SEXP _sampleIndx, SEXP _sampleName, SEXP _colnames) {


	SEXP k = allocVector(LGLSXP,1);//create logical scalar for return value

	/*
	 * convert R arguments to C type
	 */
    SEXP ans, dnms;
    const char * fName = translateChar(STRING_ELT(_fileName, 0));
    int * chnlIndx = INTEGER(_chIndx);
	int chCount = length(_chIndx);
    /*
     * determine the dataset format
     */
	hid_t       file, dataset,dataspace, memspace;         /* handles */
	hsize_t 	dimsm[2]; //dimenstions
	herr_t      status;
	unsigned nEvents;
	file = H5Fopen(fName, H5F_ACC_RDONLY, H5P_DEFAULT);
	status = H5Lexists(file, DATASETNAME3d, H5P_DEFAULT);

	int is3d;
	dataset = -1;

	if(status == TRUE){
		dataset = H5Dopen2(file, DATASETNAME3d, H5P_DEFAULT);
		dataspace = H5Dget_space(dataset);    /* dataspace handle */
		int nDim = H5Sget_simple_extent_ndims(dataspace);
		is3d = nDim == 3;
	}
	else
		is3d = 0;

	if(is3d)
	{
		/*
		 * read data from 3d mat
		 */
		int sampleIndx = INTEGER(_sampleIndx)[0];
		sampleIndx = sampleIndx -1;//convert from R to C indexing


		/*
		 * get the total number of events for the current sample
		 */
		hsize_t dims[3];
		hid_t attrID;
		status  = H5Sget_simple_extent_dims(dataspace, dims, NULL); //get dimensions of datset
		unsigned nSample = dims[0];//get total number of samples
		if(sampleIndx >= nSample)
			error("readSlice error!sample index exceeds the boundary.");
		unsigned * eCount = (unsigned *) malloc(sizeof(unsigned) * nSample);
		attrID = H5Aopen(dataset, "eventCount", H5P_DEFAULT);
		status = H5Aread(attrID, H5T_NATIVE_UINT32, eCount);
		nEvents = eCount[sampleIndx];
		free(eCount);
		H5Aclose(attrID);

		/*
		 * these two lines is the reason for the _readSlice to be inline code
		 * because we need to open hdf file to get events info
		 *
		 */
		PROTECT(ans = allocVector(REALSXP, nEvents * chCount));
		double *data_out = REAL(ans);

		/*
		 * Define the memory dataspace.
		 */
		dimsm[0] = chCount;
		dimsm[1] = nEvents;
		memspace = H5Screate_simple(2,dimsm,NULL);


		/*
		 * Define hyperslab in the dataset.
		 */
		hsize_t      count[3];              /* size of the hyperslab in the file */
		hsize_t      offset[3];             /* hyperslab offset in the file */
		hsize_t      count_out[2];          /* size of the hyperslab in memory */
		hsize_t      offset_out[2];         /* hyperslab offset in memory */

		unsigned i;
		for(i = 0; i < chCount; i++){
			int colStart = chnlIndx[i] - 1;//convert from R to C indexing
			offset[0] = sampleIndx;//start from sampleIndx-th sample
			offset[1] = colStart; //start from colStart-th channel
			offset[2] = 0; //start from the first event

			count[0]  = 1;//get one sample
			count[1]  = 1;//get one channel
			count[2]  = nEvents; //get all events


			status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,
												count, NULL);


			/*
			 * Define memory hyperslab.
			 */
			offset_out[0] = i;//start from ith column
			offset_out[1] = 0;//start from 0th event

			count_out[0]  = 1;//one channel
			count_out[1]  = nEvents; //all events
			status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL,
						 count_out, NULL);

			/*
			 * Read data from hyperslab in the file into the hyperslab in
			 * memory .
			 */
			status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
					 H5P_DEFAULT, data_out);

		}

	}
	else
	{
		/*
		 * read 2d format
		 */
		const char * sampleName = translateChar(STRING_ELT(_sampleName, 0));
		/*
		 * Open the file and the dataset.
		 */
		if(dataset>0)
		{
			//close it if it was previously opened for checking dimension
			H5Dclose(dataset);
			H5Sclose(dataspace);
		}

		status = H5Lexists(file, sampleName, H5P_DEFAULT);
		if(status == TRUE)
		{

			dataset = H5Dopen2(file, sampleName, H5P_DEFAULT);
			dataspace = H5Dget_space(dataset);    /* dataspace handle */


			/*
			 * get the total number of events for the current sample
			 */
			hsize_t dims[2];

			status  = H5Sget_simple_extent_dims(dataspace, dims, NULL); //get dimensions of datset
			nEvents = dims[1];


			PROTECT(ans = allocVector(REALSXP, nEvents * chCount));
			double *data_out = REAL(ans);

			/*
			 * Define the memory dataspace.
			 */
			dimsm[0] = chCount;
			dimsm[1] = nEvents;
			memspace = H5Screate_simple(2,dimsm,NULL);


			/*
			 * Define hyperslab in the dataset.
			 */
			hsize_t      count[2];              /* size of the hyperslab in the file */
			hsize_t      offset[2];             /* hyperslab offset in the file */
			hsize_t      count_out[2];          /* size of the hyperslab in memory */
			hsize_t      offset_out[2];         /* hyperslab offset in memory */

			unsigned i;
			for(i = 0; i < chCount; i++){
				int colStart = chnlIndx[i] - 1;//convert from R to C indexing
				offset[0] = colStart; //start from colStart-th channel
				offset[1] = 0; //start from the first event

				count[0]  = 1;//get one channel
				count[1]  = nEvents; //get all events


				status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,
													count, NULL);


				/*
				 * Define memory hyperslab.
				 */
				offset_out[0] = i;//start from ith column
				offset_out[1] = 0;//start from 0th event

				count_out[0]  = 1;//one channel
				count_out[1]  = nEvents; //all events
				status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL,
							 count_out, NULL);

				/*
				 * Read data from hyperslab in the file into the hyperslab in
				 * memory .
				 */
				status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
						 H5P_DEFAULT, data_out);

			}
		}
		else
		{
			nEvents = 0;
			PROTECT(ans = allocVector(REALSXP, nEvents * chCount));
			double *data_out = REAL(ans);
		}

	}
	/*
	 * Close/release resources.
	 */

	H5Dclose(dataset);
	H5Sclose(dataspace);
	H5Sclose(memspace);
	H5Fclose(file);
    
	/*
	 * construct R object
	 */
    PROTECT(dnms = allocVector(INTSXP, 2));
    INTEGER(dnms)[0] = nEvents;
    INTEGER(dnms)[1]=  chCount;
    setAttrib(ans,R_DimSymbol, dnms);
    /*
     * attach column names
     */
    SEXP dimnames;
//    const char * fName = translateChar(STRING_ELT(_colnames, 0));
    PROTECT(dimnames = allocVector(VECSXP, 2));
//    VECTOR(dimnames)[0] = getAttrib(x, R_NamesSymbol);
    SET_VECTOR_ELT(dimnames ,1, _colnames);
    setAttrib(ans, R_DimNamesSymbol, dimnames);
    UNPROTECT(3);//another PROTECT statement is within _readSlice call
    return(ans);
}



