#include "wrappers.h"

void ERR(int e){
	char * errmsg ="";
//	error("Error from C API:");
//	errmsg = HEstring(e);
	REprintf("hdf Error: %s\n", errmsg);
}

SEXP createFile(SEXP _fileName, SEXP _nEvent, SEXP _nChannel, SEXP _nSample, SEXP _metaSize,SEXP _compress) {

	SEXP k = allocVector(LGLSXP,1); //create logical scalar for return value


    int nSample = INTEGER(_nSample)[0];
    int nChannel = INTEGER(_nChannel)[0];
    int nEvent = INTEGER(_nEvent)[0];
//    int compress = LOGICAL(_compress)[0];

    int retval = _createFile(translateChar(STRING_ELT(_fileName, 0)), nSample, nChannel, nEvent);
    if(retval < 0)
    {
    	ERR(retval);
    	LOGICAL(k)[0] = FALSE; //set return value as FALSE
    }

    LOGICAL(k)[0] = TRUE; //set return value as TRUE
    return(k);
}


SEXP writeSlice(SEXP _fileName, SEXP _mat, SEXP _chIndx, SEXP _sample) {

	SEXP k = allocVector(LGLSXP,1);//create logical scalar for return value

    int retval, nRow, nCol, sample;
    SEXP Rdim = getAttrib(_mat, R_DimSymbol);
    nRow = INTEGER(Rdim)[0];
    nCol = INTEGER(Rdim)[1];
    sample = INTEGER(_sample)[0];


    double *mat = REAL(_mat);
    int * chIndx = INTEGER(_chIndx);
	int chCount = length(_chIndx);
    
	retval = _writeSlice(translateChar(STRING_ELT(_fileName, 0)), mat, nRow, chIndx, chCount, sample);
	if(retval < 0)
	{
		ERR(retval);
		LOGICAL(k)[0] = FALSE; //set return value as FALSE
	}

	LOGICAL(k)[0] = TRUE; //set return value as TRUE
    return(k);

}

SEXP readSlice(SEXP _fileName, SEXP _chIndx, SEXP _sample ) {


	SEXP k = allocVector(LGLSXP,1);//create logical scalar for return value

	/*
	 * convert R arguments to C type
	 */
    SEXP ans, dnms;
    int sampleIndx = INTEGER(_sample)[0];
    int * chnlIndx = INTEGER(_chIndx);
    int chCount = length(_chIndx);
    const char * fName = translateChar(STRING_ELT(_fileName, 0));


	/*
	 * the code bellow is copied from _readSlice in order to use the memory buffer
	 * dynamically allocated by R, which is transient within the call
	 */

    /*
	 * Open the file and the dataset.
	 */
	hid_t       file, dataset,dataspace, memspace;         /* handles */
	hsize_t 	dimsm[2]; //dimenstions
	herr_t      status;
	file = H5Fopen(fName, H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset = H5Dopen(file, DATASETNAME, H5P_DEFAULT);
	dataspace = H5Dget_space(dataset);    /* dataspace handle */

	sampleIndx = sampleIndx -1;//convert from R to C indexing

	/*
	 * get the total number of events for the current sample
	 */
	unsigned nEvents;
	hid_t dataset_ecountID = H5Dopen(file, "/eventCount", H5P_DEFAULT);//open ecount ds
	hid_t dataspace_ecount = H5Dget_space(dataset_ecountID); //get ds space for ecount
	//select single element slab from dataset
	hsize_t off_ecount = sampleIndx;
	hsize_t count_ecount = 1 ;
	status = H5Sselect_hyperslab(dataspace_ecount, H5S_SELECT_SET, &off_ecount, NULL, &count_ecount, NULL);
	//define memory space(single-element space) for ecount
	hsize_t dimsm_ecount = 1;
	hid_t memspace_ecount = H5Screate_simple(1, &dimsm_ecount,NULL);
	status = H5Dread(dataset_ecountID, H5T_NATIVE_UINT32
						,memspace_ecount ,dataspace_ecount, H5P_DEFAULT, &nEvents);

	/*
	 * these two lines is the reason for the _readSlice to be inline code
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


	/*
	 * Close/release resources.
	 */

	H5Dclose(dataset_ecountID);
	H5Sclose(dataspace_ecount);
	H5Dclose(dataset);
	H5Sclose(dataspace);
	H5Sclose(memspace);
	H5Sclose(memspace_ecount);
	H5Fclose(file);



//    retval = _readSlice(fName, chnlIndx, chCount, sampleIndx, data_out);
//    if(retval < 0)
//	{
//		ERR(retval);
//		LOGICAL(k)[0] = FALSE; //set return value as FALSE
//	}
    

    PROTECT(dnms = allocVector(INTSXP, 2));
    INTEGER(dnms)[0] = nEvents;
    INTEGER(dnms)[1]=  chCount;
    setAttrib(ans,R_DimSymbol, dnms);
    UNPROTECT(2);//another PROTECT statement is within _readSlice call
    return(ans);
}



