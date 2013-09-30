#include "hdfFlow.h"

/*
 * create hdf file for flow data
 */
herr_t _createFile(const char * fName, unsigned nSample, unsigned nChnl, unsigned nEvt){
	hid_t       file_id, dataset_id, dataspace_id, dataspace_attr_id, attribute_id;  /* identifiers */
	hsize_t     dims[3], dim_attr;
	herr_t      status;

	/* Create a new file using default properties. */
	file_id = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* Create the data space for the 3d mat. */
	dims[0] = nSample;
	dims[1] = nChnl;
	dims[2] = nEvt;
	dataspace_id = H5Screate_simple(3, dims, NULL);

	/* Create the 3d mat. */
	dataset_id = H5Dcreate2(file_id, DATASETNAME, H5T_IEEE_F64LE_g, dataspace_id,
						  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//set it to use chunking
	hsize_t		chunk_dims[3] = {1, 1, nEvt};
	hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(dcpl_id, 3, chunk_dims);

	/* Create the data space for the attribute. */
	dim_attr = nSample;
	dataspace_attr_id = H5Screate_simple(1, &dim_attr, NULL);

	/* Create a eventCount attribute. */
	attribute_id = H5Acreate2 (dataset_id, "eventCount", H5T_STD_U32LE, dataspace_attr_id,
							 H5P_DEFAULT, H5P_DEFAULT);

	/* Close the attribute. */
	status = H5Aclose(attribute_id);

	/* End access to the dataset and release resources used by it. */
	status = H5Dclose(dataset_id);


	/* Terminate access to the data space. */
	status = H5Sclose(dataspace_id);
	status = H5Sclose(dataspace_attr_id);

	/* Close the file. */
	status = H5Fclose(file_id);

	return status;
}

/*
 * write one sample to hdf
 *
 * each channel  is stored as a data chunk(size= events),
 * which is indexed in file and fast to access as a whole unit(slice),
 * so the C API for accessing events are not provided due to the inefficient IO.
*/
herr_t _writeSlice(const char * fName, double * mat, unsigned nEvents, unsigned * chnlIndx, unsigned chCount,  unsigned sampleIndx){
	/*
	 * Open the file and the dataset.
	 */
	hid_t       file, dataset,dataspace, memspace, attrID;         /* handles */
	hsize_t 	dims[3], dimsm[2]; //dimenstions
	herr_t      status;
	file = H5Fopen(fName, H5F_ACC_RDWR, H5P_DEFAULT);//open file
	dataset = H5Dopen(file, DATASETNAME, H5P_DEFAULT);//open dataset
	dataspace = H5Dget_space(dataset);    /* dataspace handle */

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
	hsize_t      count_in[3];          /* size of the hyperslab in memory */
	hsize_t      offset_in[2];         /* hyperslab offset in memory */

	/*
	 * write subset
	 */
	unsigned i;
	for(i = 0; i < chCount; i++){
		int colStart = chnlIndx[i];
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
		 * Read data from hyperslab in the file into the hyperslab in
		 * memory .
		 */
		status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
				 H5P_DEFAULT, mat);

	}

	/*
	 * get eCount attribute
	 */
	status  = H5Sget_simple_extent_dims(dataspace, dims, NULL); //get dimensions of datset
	unsigned nSample = dims[0];//get total number of samples
	unsigned * eCount = (unsigned *) malloc(nSample * sizeof(unsigned));
	attrID = H5Aopen(dataset, "eventCount", H5P_DEFAULT);
	status = H5Aread(attrID, H5T_NATIVE_UINT32, eCount);
	//update the eCount for current sample
	eCount[sampleIndx] = nEvents;
	/*
	 * write back to hdf
	 */
	status = H5Awrite(attrID, H5T_NATIVE_UINT32, eCount);

	/*
	 * Close/release resources.
	 */

	H5Aclose(attrID);
	H5Dclose(dataset);
	H5Sclose(dataspace);
	H5Sclose(memspace);
	H5Fclose(file);

	free(eCount);
	return status;

}

herr_t _readSlice(const char * fName, unsigned * chnlIndx, unsigned chCount, unsigned sampleIndx, double * data_out) {
	/*
	 * Open the file and the dataset.
	 */
    hid_t       file, dataset,dataspace, memspace, attrID;         /* handles */
    hsize_t 	dims[3], dimsm[2]; //dimenstions
    herr_t      status;
    file = H5Fopen(fName, H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset = H5Dopen(file, DATASETNAME, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);    /* dataspace handle */

    /*
     * get the total number of events for the current sample
     */
    status  = H5Sget_simple_extent_dims(dataspace, dims, NULL); //get dimensions of datset
    unsigned nSample = dims[0];//get total number of samples
    unsigned * eCount = (unsigned *) malloc(nSample * sizeof(unsigned));
    attrID = H5Aopen(dataset, "eventCount", H5P_DEFAULT);
    status = H5Aread(attrID, H5T_NATIVE_UINT32, eCount);
    unsigned nEvents = eCount[sampleIndx];

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
		int colStart = chnlIndx[i];
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

	H5Aclose(attrID);
	H5Dclose(dataset);
	H5Sclose(dataspace);
	H5Sclose(memspace);
	H5Fclose(file);

	free(eCount);

	return status;
}

