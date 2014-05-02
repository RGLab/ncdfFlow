#include "hdfFlow.h"
herr_t _createFile2d(const char * fName){
	/* Create a new file using default properties. */
	hid_t file_id = H5Fcreate(fName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if(file_id<0)
		return(file_id);
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
 * write one sample to hdf
 *
 * each channel  is stored as a data chunk(size= events),
 * which is indexed in file and fast to access as a whole unit(slice),
 * so the C API for accessing individual events are not provided due to the inefficient IO.
*/
herr_t _writeSlice(const char * fName, double * mat, unsigned nEvents, unsigned * chnlIndx, unsigned chCount,  unsigned sampleIndx){
	/*
	 * Open the file and the dataset.
	 */
	hid_t  file, dataset,dataspace, memspace;         /* handles */
	hid_t attrID;
	herr_t      status;
	file = H5Fopen(fName, H5F_ACC_RDWR, H5P_DEFAULT);//open file
	dataset = H5Dopen2(file, DATASETNAME3d, H5P_DEFAULT);//open dataset
	dataspace = H5Dget_space(dataset);    /* dataspace handle */

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
		int colStart = chnlIndx[i] -1; //convert from R to C indexing
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


//	/*
//	 * update the eCount for current sample
//	 */
//	hid_t dataset_ecountID = H5Dopen2(file, "/eventCount", H5P_DEFAULT);//open ecount ds
//	hid_t dataspace_ecount = H5Dget_space(dataset_ecountID); //get ds space for ecount
//	//select single element slab from dataset
//	hsize_t off_ecount = sampleIndx;
//	hsize_t count_ecount =1 ;
//	status = H5Sselect_hyperslab(dataspace_ecount, H5S_SELECT_SET, &off_ecount, NULL, &count_ecount, NULL);
//	//define memory space(single-element space) for ecount
//	hsize_t dimsm_ecount = 1;
//	hid_t memspace_ecount = H5Screate_simple(1, &dimsm_ecount,NULL);
//	status = H5Dwrite(dataset_ecountID, H5T_NATIVE_UINT32,memspace_ecount ,dataspace_ecount, H5P_DEFAULT, &nEvents);

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


	/*
	 * Close/release resources.
	 */


//	H5Dclose(dataset_ecountID);
	H5Dclose(dataset);
//	H5Sclose(dataspace_ecount);
	H5Sclose(dataspace);
	H5Sclose(memspace);
//	H5Sclose(memspace_ecount);
	H5Fclose(file);


	return status;

}
/*
 * For the simplicity reason, we allocated data_out outside of routine
 * when tested in in pure C environment so that we don't have to deal with R C API: allocVector
 *
 * But when called from R wrapper, we need to do in inline
 * because the size of this memory buffer is yet to be determined
 * by reading total number of events from hdf file. (And we don't want to open
 * hdf file twice)
 *
 *
 */
herr_t _readSlice(const char * fName, unsigned * chnlIndx, unsigned chCount, unsigned sampleIndx, double * data_out)
{
	/*
	 * Open the file and the dataset.
	 */
    hid_t       file, dataset,dataspace, memspace;         /* handles */
    hsize_t 	dimsm[2]; //dimenstions
    herr_t      status;
    file = H5Fopen(fName, H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset = H5Dopen2(file, DATASETNAME3d, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);    /* dataspace handle */

    sampleIndx = sampleIndx -1;//convert from R to C indexing

    /*
     * get the total number of events for the current sample
     */
//    unsigned nEvents;
//    hid_t dataset_ecountID = H5Dopen2(file, "/eventCount", H5P_DEFAULT);//open ecount ds
//	hid_t dataspace_ecount = H5Dget_space(dataset_ecountID); //get ds space for ecount
//	//select single element slab from dataset
//	hsize_t off_ecount = sampleIndx;
//	hsize_t count_ecount = 1 ;
//	status = H5Sselect_hyperslab(dataspace_ecount, H5S_SELECT_SET, &off_ecount, NULL, &count_ecount, NULL);
//	//define memory space(single-element space) for ecount
//	hsize_t dimsm_ecount = 1;
//	hid_t memspace_ecount = H5Screate_simple(1, &dimsm_ecount,NULL);
//	status = H5Dread(dataset_ecountID, H5T_NATIVE_UINT32
//						,memspace_ecount ,dataspace_ecount, H5P_DEFAULT, &nEvents);

    /*
	 * get the total number of events for the current sample
	 */
   hsize_t dims[3];
   hid_t attrID;
   status  = H5Sget_simple_extent_dims(dataspace, dims, NULL); //get dimensions of datset
   unsigned nSample = dims[0];//get total number of samples
   unsigned * eCount = (unsigned *) malloc(nSample * sizeof(unsigned));
   attrID = H5Aopen(dataset, "eventCount", H5P_DEFAULT);
   status = H5Aread(attrID, H5T_NATIVE_UINT32, eCount);
   unsigned nEvents = eCount[sampleIndx];
   free(eCount);
   H5Aclose(attrID);

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

//	H5Dclose(dataset_ecountID);
//	H5Sclose(dataspace_ecount);
	H5Dclose(dataset);
	H5Sclose(dataspace);
	H5Sclose(memspace);
//	H5Sclose(memspace_ecount);
	H5Fclose(file);



	return status;
}

