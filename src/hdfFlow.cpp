#include "hdfFlow.h"
#include <boost/lexical_cast.hpp>


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
	Rcpp::stop("hdf Error");
    return 0;

}


herr_t _createFile2d(const char * fName, H5F_libver_t libver_lowbound){
	/* Create a new file using default properties. */
	hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_libver_bounds(fapl, libver_lowbound, H5F_LIBVER_LATEST);
	hid_t file_id = H5Fcreate(fName, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);

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

// [[Rcpp::export]]
bool createFile(std::string filename, int nEvent, int nChannel, int nSample, int nDim, int nCompressionRatio, bool is_libver_earliest = false)
{

	H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t)custom_print_cb, NULL);
	int retval;
    if(nDim == 3)
		retval = _createFile3d(filename.c_str(), nSample, nChannel, nEvent, nCompressionRatio);
    else
    	retval = _createFile2d(filename.c_str(), is_libver_earliest?H5F_LIBVER_EARLIEST:H5F_LIBVER_LATEST);

    return retval >= 0;

}
/*
 * We separate the open and read file so that multi-sample read operations can share the same file id.
 *
 */
// [[Rcpp::export]]
void open_hdf(std::string filename, unsigned flags, hid_t & fileid, hid_t & dataset, hid_t & dataspace, bool & is3d)
{
  H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t)custom_print_cb, NULL);

  /*
   * determine the dataset format
   */
  herr_t      status;

  fileid = H5Fopen(filename.c_str(), flags, H5P_DEFAULT);
  status = H5Lexists(fileid, DATASETNAME3d, H5P_DEFAULT);


  dataset = -1;

  if(status == TRUE){
    dataset = H5Dopen2(fileid, DATASETNAME3d, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);    /* dataspace handle */
  int nDim = H5Sget_simple_extent_ndims(dataspace);
  is3d = nDim == 3;
  }
  else
    is3d = false;

}
/*
 * inline _writeSlice and _writeSlice2d code
 */
// [[Rcpp::export]]
bool writeSlice(std::string filename, Rcpp::NumericMatrix data, std::vector<int> chIndx, int sampleIndx, int nRatio) {

	H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t)custom_print_cb, NULL);

	double *mat = REAL(data.get__());

	int chCount = chIndx.size();


    unsigned nEvents = data.rows();

    sampleIndx--;//convert from R to C indexing
	/*
	 * Open the file and the dataset.
	 */

	herr_t      status;
	hid_t  fileid, dataset,dataspace, memspace;
	bool is3d;
	open_hdf(filename, H5F_ACC_RDWR, fileid, dataset, dataspace, is3d);

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
		int i;

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
		  int nSample = dims[0];//get total number of samples
		  if(sampleIndx >= nSample)
		  		Rcpp::stop("writeSlice error!sample index exceeds the boundary.");
		  unsigned * eCount = new unsigned[nSample];
		  attrID = H5Aopen(dataset, "eventCount", H5P_DEFAULT);
		  status = H5Aread(attrID, H5T_NATIVE_UINT32, eCount);
		  //update the eCount for current sample
		  eCount[sampleIndx] = nEvents;
		  /*
		   * write back to hdf
		   */
		  status = H5Awrite(attrID, H5T_NATIVE_UINT32, eCount);

		  delete [] eCount;
		  H5Aclose(attrID);
	}
	else
	{
		//convert index to string to be used as dataset name
		std::string sampleName = boost::lexical_cast<std::string>(sampleIndx);
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
		status = H5Lexists(fileid, sampleName.c_str(), H5P_DEFAULT);
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
			status = H5Pset_deflate (dcpl_id, nRatio);

			/* Create the 2d mat. */
			dataset = H5Dcreate2(fileid, sampleName.c_str(), H5T_IEEE_F32LE_g, dataspace,
								  H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
			H5Pclose(dcpl_id);
		}
		else
		{
			dataset = H5Dopen2(fileid, sampleName.c_str(), H5P_DEFAULT);
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
		int i;

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

	H5Fclose(fileid);




	return status >= 0;


}

/*
 * these 4 APIs (open_hdf, get_event_number, readSlice, close_hdf) used to belong to a single routine
 * They are separated for the purpose of:
 * 1. to be used in two environment: R or pure c++ (exposed to other R packages through R_RegisterCCallable)
 * 2. open file once and readSlice multiple times
 *
 * So they need to be used together to perform the complete and clean read.
 * e.g. readSlice is responsible to close dataset and dataspace opened by open_hdf and get_event_number
 */


/*
 * We separate the open and read file so that multi-sample read operations can share the same file id.
 *
 */

// [[Rcpp::export]]
unsigned get_event_number(hid_t fileid, hid_t & dataset, hid_t & dataspace, unsigned sampleIndx, bool is3d)
{
	unsigned nEvents;
	if(is3d)
	{

		/*
		 * get the total number of events for the current sample
		 */
		hsize_t dims[3];
		hid_t attrID;
		H5Sget_simple_extent_dims(dataspace, dims, NULL); //get dimensions of datset
		unsigned nSample = dims[0];//get total number of samples
		if(sampleIndx >= nSample)
			Rcpp::stop("readSlice error!sample index exceeds the boundary.");
		unsigned * eCount = (unsigned *) malloc(sizeof(unsigned) * nSample);
		attrID = H5Aopen(dataset, "eventCount", H5P_DEFAULT);
		H5Aread(attrID, H5T_NATIVE_UINT32, eCount);
		nEvents = eCount[sampleIndx];
		free(eCount);
		H5Aclose(attrID);
	}
	else
	{
		/*
		 * convert index to string to be used as dataset name
		 * because dataset can not be renamed once created in hdf
		 */
		std::string sampleName = boost::lexical_cast<std::string>(sampleIndx);
		/*
		 * Open the file and the dataset.
		 */
		if(dataset>0)
		{
			//close it if it was previously opened for checking dimension
			H5Dclose(dataset);
			H5Sclose(dataspace);
		}

		herr_t status = H5Lexists(fileid, sampleName.c_str(), H5P_DEFAULT);
		if(status == TRUE)
		{

			dataset = H5Dopen2(fileid, sampleName.c_str(), H5P_DEFAULT);
			dataspace = H5Dget_space(dataset);    /* dataspace handle */


			/*
			 * get the total number of events for the current sample
			 */
			hsize_t dims[2];

			status  = H5Sget_simple_extent_dims(dataspace, dims, NULL); //get dimensions of datset
			nEvents = dims[1];
		}
		else
			nEvents = 0;
	}
	return nEvents;
}
/*
 * since events number is stored in hdf we have to separate get_event_number call from readSlice so that we can separate the logic of
 *  dynamic allocation of memory buffer either through Rcpp (Rcpp::NumericVector mat(nEvents * nCh))
 *  or pure C ( double * mat = new double[nEvents * nCh]) thus we have the option to use them in the pure c++ code evnironment
 * Thus, dataset and dataspace are assumed to be already opened by get_event_number call and there is no validity check in readSlice
 * and up to caller to do this check carefully
 */

// [[Rcpp::plugins(hdf5)]]
// [[Rcpp::depends(BH,RcppArmadillo)]]
void readSlice_cpp(hid_t fileid, hid_t dataset, hid_t dataspace
								, std::vector<unsigned> chIndx
								, unsigned sampleIndx
								, unsigned nEvents
								, double * data_out
								, bool is3d
								)
{

//	H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t)custom_print_cb, NULL);

	hid_t   memspace;         /* handles */
	hsize_t 	dimsm[2]; //dimenstions
	herr_t      status;
	int chCount = chIndx.size();

   	/*
	 * read data from 3d mat
	 */
	if(is3d)
	{

		/*
		 * get the total number of events for the current sample
		 */
		hsize_t dims[3];
		status  = H5Sget_simple_extent_dims(dataspace, dims, NULL); //get dimensions of datset
		unsigned nSample = dims[0];//get total number of samples
		if(sampleIndx >= nSample)
			Rcpp::stop("readSlice error!sample index exceeds the boundary.");


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

		int i;
		for(i = 0; i < chCount; i++){
			int colStart = chIndx.at(i);
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

			int i;
			for(i = 0; i < chCount; i++){
				int colStart = chIndx.at(i);
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
	H5Dclose(dataset);
	H5Sclose(dataspace);
	H5Sclose(memspace);
}

// [[Rcpp::export]]
void close_hdf(hid_t fileid)
{

	H5Fclose(fileid);

}

/*
 * the reason for his wrapper is that Rcpp does not support SEXP to double * conversion
 */
// [[Rcpp::export]]
void readSlice(hid_t fileid, hid_t dataset, hid_t dataspace
								, std::vector<unsigned> chIndx
								, unsigned sampleIndx
								, unsigned nEvents
								, Rcpp::NumericVector data_out
								, bool is3d
								)
{
	double * data = REAL(data_out.get__());
	readSlice_cpp(fileid, dataset, dataspace, chIndx, sampleIndx, nEvents, data, is3d);
}

