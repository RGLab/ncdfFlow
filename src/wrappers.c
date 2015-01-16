#include "wrappers.h"

#define MSG_SIZE       1024

herr_t my_hdf5_error_handler(unsigned n, const H5E_error2_t *err_desc, void *client_data)
{
	char                maj[MSG_SIZE];
	char                min[MSG_SIZE];

	const int		indent = 4;

	if(H5Eget_msg(err_desc->maj_num, NULL, maj, MSG_SIZE)<0)
		return -1;

	if(H5Eget_msg(err_desc->min_num, NULL, min, MSG_SIZE)<0)
		return -1;

	REprintf("%*s error #%03d: in %s(): line %u\n",
		 indent, "", n, err_desc->func_name, err_desc->line);
	REprintf("%*smajor: %s\n", indent*2, "", maj);
	REprintf("%*sminor: %s\n", indent*2, "", min);

   return 0;
}

/*
 * customize the printing function so that it print to R error console
 * also raise the R error once the error stack printing is done
 */
herr_t custom_print_cb(hid_t estack, void *client_data)
{
	hid_t estack_id = H5Eget_current_stack();//copy stack before it is corrupted by my_hdf5_error_handler
	H5Ewalk2(estack_id, H5E_WALK_DOWNWARD, my_hdf5_error_handler, client_data);
	H5Eclose_stack(estack_id);
	error("hdf Error");
    return 0;

}


herr_t _createFile2d(const char * fName){
	/* Create a new file using default properties. */
	hid_t file_id = H5Fcreate(fName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* Close the file. */
	herr_t status = H5Fclose(file_id);

	return status;
}

/*
 * create hdf file for flow data
 *
 * Thus use dataset instead of attribute to store events count
 * because we want to avoid dynamically allocating 1d array for it
 *
 * In order to keep it back-compatible, we have to stick to the attribute for eCount
 */
herr_t _createFile3d(const char * fName, unsigned nSample, unsigned nChnl, unsigned nEvt, unsigned nRatio){
//	hid_t       file_id, dataset_id, dataspace_id, dataspace_ecount_id,dataset_ecount_id;  /* identifiers */
	hid_t       file_id, dataset_id, dataspace_id, dataspace_attr_id, attribute_id;  /* identifiers */
	hsize_t     dims[3], dim_attr;
	herr_t      status;

	/* Create a new file using default properties. */
	file_id = H5Fcreate(fName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* Create the data space for the 3d mat. */
	dims[0] = nSample;
	dims[1] = nChnl;
	dims[2] = nEvt;
	dataspace_id = H5Screate_simple(3, dims, NULL);

	hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);

	//set it to use chunking
	hsize_t		chunk_dims[3] = {1, 1, nEvt};
	H5Pset_chunk(dcpl_id, 3, chunk_dims);
	//set it to use compression (zlib)
	status = H5Pset_deflate (dcpl_id, nRatio);

	/* Create the 3d mat. */
	dataset_id = H5Dcreate2(file_id, DATASETNAME3d, H5T_IEEE_F32LE_g, dataspace_id,
						  H5P_DEFAULT, dcpl_id, H5P_DEFAULT);


	/* Create the data space for the attribute. */
	dim_attr = nSample;
	dataspace_attr_id = H5Screate_simple(1, &dim_attr, NULL);
	attribute_id = H5Acreate2 (dataset_id, "eventCount", H5T_STD_U32LE, dataspace_attr_id,
	               H5P_DEFAULT, H5P_DEFAULT);
//	dataspace_ecount_id = H5Screate_simple(1, &dim_attr, NULL);
//
//	/* Create the 1d ds for eventCount */
//	dataset_ecount_id = H5Dcreate2(file_id, "/eventCount", H5T_STD_U32LE, dataspace_ecount_id,
//							  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);




	/* End access to the dataset and release resources used by it. */
	status = H5Dclose(dataset_id);
	status = H5Aclose(attribute_id);
//	status = H5Dclose(dataset_ecount_id);

	H5Pclose(dcpl_id);

	/* Terminate access to the data space. */
	status = H5Sclose(dataspace_id);
	status = H5Sclose(dataspace_attr_id);
//	status = H5Sclose(dataspace_ecount_id);

	/* Close the file. */
	status = H5Fclose(file_id);

	return status;
}

/*
 * metaSize and compress are currently not used
 */
SEXP createFile(SEXP _fileName, SEXP _nEvent, SEXP _nChannel, SEXP _nSample, SEXP _dim, SEXP _ratio) {

	H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t)custom_print_cb, NULL);

	SEXP k = allocVector(LGLSXP,1); //create logical scalar for return value
    int nDim = INTEGER(_dim)[0];

    int retval;
    if(nDim == 3)
    {
    	int nSample = INTEGER(_nSample)[0];
    	int nChannel = INTEGER(_nChannel)[0];
		int nEvent = INTEGER(_nEvent)[0];
		//compression
		int nCmpRatio = INTEGER(_ratio)[0];
		retval = _createFile3d(translateChar(STRING_ELT(_fileName, 0)), nSample, nChannel, nEvent, nCmpRatio);
    }
    else
    {

    	retval = _createFile2d(translateChar(STRING_ELT(_fileName, 0)));
    }

    LOGICAL(k)[0] = retval >= 0;
    return(k);
}

/*
 * inline _writeSlice and _writeSlice2d code
 */
SEXP writeSlice(SEXP _fileName, SEXP _mat, SEXP _chIndx, SEXP _sampleIndx, SEXP _ratio) {

	H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t)custom_print_cb, NULL);

	SEXP k = allocVector(LGLSXP,1);//create logical scalar for return value
	const char * fName = translateChar(STRING_ELT(_fileName, 0));
	double *mat = REAL(_mat);
	int * chIndx = INTEGER(_chIndx);
	int chCount = length(_chIndx);


    int nRow;
    SEXP Rdim = getAttrib(_mat, R_DimSymbol);
    nRow = INTEGER(Rdim)[0];

    int nEvents = nRow;
    int sampleIndx = INTEGER(_sampleIndx)[0];
    sampleIndx = sampleIndx -1;//convert from R to C indexing
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
		//convert index to string to be used as dataset name
		char * sampleName = (char *)malloc(sizeof(char)*MAXLEN);
		snprintf(sampleName, MAXLEN, "%d", sampleIndx);
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

			//compression
			int nCmpRatio = INTEGER(_ratio)[0];
			status = H5Pset_deflate (dcpl_id, nCmpRatio);

			/* Create the 2d mat. */
			dataset = H5Dcreate2(file, sampleName, H5T_IEEE_F32LE_g, dataspace,
								  H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
			H5Pclose(dcpl_id);
		}
		else
		{
			dataset = H5Dopen2(file, sampleName, H5P_DEFAULT);
			dataspace = H5Dget_space(dataset);    /* dataspace handle */
		}
		free(sampleName);
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




	LOGICAL(k)[0] = status >= 0;
	return(k);

}

//_readSlice is inlined in readFrame function


