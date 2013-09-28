#include "hdfFlow.h"

/*
 * create hdf file for flow data
 */
herr_t _createFile(const char * fName, unsigned nSample, unsigned nChnl, unsigned nEvt){
	hid_t       file_id, dataset_id, dataspace_id, dataspace_attr_id, attribute_id;  /* identifiers */
	hsize_t     dims[3], dim_attr;
	herr_t      status;

	/* Create a new file using default properties. */
	file_id = H5Fcreate("test.h5", H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

	/* Create the data space for the 3d mat. */
	dims[0] = nSample;
	dims[1] = nChnl;
	dims[2] = nEvt;
	dataspace_id = H5Screate_simple(3, dims, NULL);

	/* Create the 3d mat. */
	dataset_id = H5Dcreate2(file_id, "/exprsMat", H5T_IEEE_F64LE_g, dataspace_id,
						  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	/* Create the data space for the attribute. */
	dim_attr = nSample;
	dataspace_attr_id = H5Screate_simple(1, &dim_attr, NULL);

	/* Create a eventCount attribute. */
	attribute_id = H5Acreate2 (dataset_id, "eventCount", H5T_STD_I32BE, dataspace_attr_id,
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
void _writeSlice(const char * fName, double * mat, unsigned chnlIndx,  unsigned sampleIndx){

}

double * _readSlice(const char * fName, unsigned chnlIndx, unsigned sampleIndx) {
	double * mat;
	return mat;
}

