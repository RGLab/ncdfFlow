#include <Rcpp.h>
#include <Rinternals.h>
#include "wrappers.h"
#include <boost/lexical_cast.hpp>

typedef std::vector<std::string> strVec;
typedef std::vector<int> intVec;
typedef std::vector<unsigned> uintVec;
typedef std::vector<bool> boolVec;

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
	Rcpp::stop("hdf Error");
    return 0;

}


// [[Rcpp::plugins(hdf5)]]
// [[Rcpp::depends(BH)]]
Rcpp::NumericVector readSlice_cpp(std::string fName
								, std::vector<unsigned> chIndx
								, unsigned sampleIndx
//								, Rcpp::StringVector colnames
								)
{

	H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t)custom_print_cb, NULL);


    Rcpp::NumericVector res;

	int chCount = chIndx.size();
    /*
     * determine the dataset format
     */
	hid_t       file, dataset,dataspace, memspace;         /* handles */
	hsize_t 	dimsm[2]; //dimenstions
	herr_t      status;
	unsigned nEvents;
	file = H5Fopen(fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
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
	/*
	 * read data from 3d mat
	 */
	if(is3d)
	{

		/*
		 * get the total number of events for the current sample
		 */
		hsize_t dims[3];
		hid_t attrID;
		status  = H5Sget_simple_extent_dims(dataspace, dims, NULL); //get dimensions of datset
		unsigned nSample = dims[0];//get total number of samples
		if(sampleIndx >= nSample)
			Rcpp::stop("readSlice error!sample index exceeds the boundary.");
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

		res = Rcpp::NumericVector(nEvents * chCount);
		double *data_out = REAL(res.get__());


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
		H5Dclose(dataset);
		H5Sclose(dataspace);
		H5Sclose(memspace);
	}
	else
	{
		/*
		 * read 2d format
		 */

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

		status = H5Lexists(file, sampleName.c_str(), H5P_DEFAULT);
		if(status == TRUE)
		{

			dataset = H5Dopen2(file, sampleName.c_str(), H5P_DEFAULT);
			dataspace = H5Dget_space(dataset);    /* dataspace handle */


			/*
			 * get the total number of events for the current sample
			 */
			hsize_t dims[2];

			status  = H5Sget_simple_extent_dims(dataspace, dims, NULL); //get dimensions of datset
			nEvents = dims[1];



			res = Rcpp::NumericVector(nEvents * chCount);
			double *data_out = REAL(res.get__());

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

			H5Sclose(dataspace);
			H5Sclose(memspace);
			H5Dclose(dataset);
		}
		else
		{
			nEvents = 0;

			res = Rcpp::NumericVector(nEvents * chCount);

		}

	}



	H5Fclose(file);


//	//set dims
//	Rcpp::IntegerVector dims(2);
//	dims[0] = nEvents;
//	dims[1]=  chCount;
//	res.attr("dim") = dims;
//	// attach column names
//	Rcpp::List dimnms = Rcpp::List::create(R_NilValue, colnames);
//	res.attr("dimnames") = dimnms;
//

    return(res);
}




// [[Rcpp::export]]
Rcpp::S4 readFrame(Rcpp::S4 x
					, Rcpp::RObject i_obj
					, Rcpp::RObject j_obj
					, bool useExpr
					)
{
	/*
	 * parse i index (sample name)
	 */
	std::string sampleName;
	unsigned short i_type = i_obj.sexp_type();

	if(i_type  == STRSXP)
	{
		sampleName = Rcpp::as<std::string>(i_obj.get__());
	}
	else if(i_type == REALSXP || i_type == INTSXP)
	{
		unsigned s_ind = Rcpp::as<unsigned>(i_obj.get__());
		s_ind = s_ind - 1;
		//get local sample vec
		Rcpp::S4 fspd = x.slot("phenoData");
		Rcpp::DataFrame fsdata = fspd.slot("data");
		Rcpp::CharacterVector sn = fsdata["name"];
		sampleName =  sn(s_ind);
	}
	else
		Rcpp::stop("unsupported i type!");


	/*
	 * parse j index (channel names)
	 */
	Rcpp::Environment frEnv = x.slot("frames");
	Rcpp::S4 frObj = frEnv.get(sampleName);
	Rcpp::S4 fr = Rcpp::clone(frObj);

	  //get local channel names
	  Rcpp::StringVector colnames = x.slot("colnames");

	  Rcpp::StringVector ch_selected;
	 /*
	  * subset by j if applicable
	  */
	 int j_type = j_obj.sexp_type();
	 //creating j index used for subsetting colnames and pdata
	 Rcpp::IntegerVector j_indx;

	 if(j_type == STRSXP)//when character vector
	 {
		 ch_selected = Rcpp::StringVector(j_obj.get__());
		 unsigned nCol = ch_selected.size();
		 j_indx = Rcpp::IntegerVector(nCol);
		 //match ch_selected to colnames
		for(unsigned i = 0 ; i < nCol; i ++)
		{
			const Rcpp::internal::string_proxy<STRSXP> &thisCh = ch_selected(i);
			Rcpp::StringVector::iterator match_id = std::find(colnames.begin(), colnames.end(), thisCh);
			if(match_id == colnames.end()){
				std::string strCh = Rcpp::as<std::string>(thisCh);
				Rcpp::stop("j subscript out of bounds: " + strCh);
			}else
			{
				j_indx(i) = match_id - colnames.begin();
			}
		}
	 }
	 else if(j_type == NILSXP)//j is set to NULL in R when not supplied
	 {
		 ch_selected = colnames;
	 }
	 else if(j_type == LGLSXP)
	 {
		 Rcpp::LogicalVector j_val(j_obj.get__());
		 //convert to integer vector
		 unsigned nSel = Rcpp::sum(j_val);
		 unsigned thisCount = 0;
		 j_indx = Rcpp::IntegerVector(nSel);
		 for(unsigned i = 0; i < j_val.size(); i++){
			 if(j_val(i)){
				 j_indx(thisCount++) = i;
			 }

		 }

		 ch_selected = colnames[j_indx];
	 }
	 else if(j_type == INTSXP || j_type == REALSXP)
	 {
		 Rcpp::IntegerVector j_val(j_obj.get__());
		 j_indx = j_val - 1; //convert to 0-based index
		 ch_selected = colnames[j_indx];
	 }
	 else
		 Rcpp::stop("unsupported j expression!");
	/*
	 * update annotationDataFrame
	 * we don't update description slot(i.e. keywords) as flowCore does
	 */
	 if(j_type != NILSXP)
	 {
		Rcpp::S4 pheno = fr.slot("parameters");
		Rcpp::DataFrame pData = pheno.slot("data");

		Rcpp::CharacterVector pd_rn = pData.attr("row.names");
		Rcpp::CharacterVector pd_name = pData["name"];
		Rcpp::CharacterVector pd_desc = pData["desc"];
		Rcpp::NumericVector pd_range = pData["range"];
		Rcpp::NumericVector pd_minRange = pData["minRange"];
		Rcpp::NumericVector pd_maxRange = pData["maxRange"];

		Rcpp::List plist = Rcpp::List::create(Rcpp::Named("name") = pd_name[j_indx]
													,Rcpp::Named("desc") = pd_desc[j_indx]
													,Rcpp::Named("range") = pd_range[j_indx]
													,Rcpp::Named("minRange") = pd_minRange[j_indx]
													,Rcpp::Named("maxRange") = pd_maxRange[j_indx]
													);
		plist.attr("class") = "data.frame";
		plist.attr("row.names") = pd_rn[j_indx];
		pheno.slot("data") = plist;
	 }

	 /*
	  * read data from hdf
	  */
	if(useExpr){
		Rcpp::StringVector origChNames = x.slot("origColnames");
		unsigned nCh = ch_selected.size();
		unsigned nOrig = origChNames.size();

		//convert chnames to global channel index
		std::vector<unsigned> chIndx(nCh);

		for(unsigned i = 0; i < nCh; i++){
			Rcpp::String thisCh = ch_selected(i);

			for(unsigned j = 0; j < nOrig; j++)
			{
//				std::string thisOrig = Rcpp::as<std::string>(origChNames[j]);
				if(thisCh == origChNames(j))
				{
					chIndx.at(i) = j;
					break;
				}
			}
//			Rcpp::Rcout << thisCh << ":" << chIndx.at(i) << std::endl;
		}



	    Rcpp::Environment IndiceEnv = x.slot("indices");
	    Rcpp::RObject Indice = IndiceEnv.get(sampleName);
	    bool subByIndice;
	    if(Indice.isNULL())
	    	Rcpp::stop("Invalid sample name '" + sampleName + "'! It is not found in 'indices' slot!");
	    else{


	    	if(Indice.sexp_type() == LGLSXP)//somehow NA is of logical type by default in R
			{
	    		subByIndice = false;
//	    		Rcpp::Rcout <<"no subsetting"<< std::endl;

			}
	    	else if(Indice.sexp_type() == RAWSXP){
	    		subByIndice = true;
//	    		Rcpp::Rcout <<"subsetting"<< std::endl;
	    	}else
	    		Rcpp::stop("Invalid indices type '" + sampleName + "'!It must be raw vector !");

	    }


//	      get sample index
	      unsigned samplePos;
	      //convert Rcpp vector to std vector
	      Rcpp::StringVector origSampVec = x.slot("origSampleVector");
	      unsigned nSample = origSampVec.size();
	      std::vector<std::string> origSampleVector(nSample);
	      for(unsigned i = 0; i < nSample; i++)
	    	  origSampleVector.at(i) = Rcpp::as<std::string>(origSampVec[i]);

	      std::vector<std::string>::iterator it = std::find(origSampleVector.begin(), origSampleVector.end(), sampleName);

	      if(it == origSampleVector.end())
	    	  Rcpp::stop("Invalid sample name '" + sampleName + "'! It is not found in 'origSampleVector' slot!");
	      else
	      {
	    	  samplePos = it - origSampleVector.begin();
//	    	  Rcpp::Rcout << samplePos << std::endl;
	      }
	      Rcpp::String file = x.slot("file");

	      Rcpp::NumericVector mat = readSlice_cpp(file, chIndx, samplePos);


//	      subset data by indices if necessary
	      if(subByIndice){
	    	  /*
	    	   * convert bytes to bool vector
	    	   */
	    	  Rcpp::RawVector bytes(Indice.get__());
	    	  unsigned len = bytes.attr("bitlen");
	    	  Rcpp::LogicalVector indx(len * nCh);

			  unsigned byteIndex, bitIndex;
			  for(unsigned i =0 ; i < len; i++)
			  {
				  byteIndex = i / 8;
				  bitIndex = i % 8;
				  //mark the event index for each col
				  for(unsigned j = 0; j < nCh; j++)
					  indx[i+ j * len] = bytes[byteIndex] & 1 << bitIndex;//IS_SET(bytes, byteIndex, bitIndex);
			  }

			 mat = mat[indx];

	      }

		//set dims
		Rcpp::IntegerVector dims(2);
		dims[0] = mat.size()/nCh;
		dims[1]=  nCh;
		mat.attr("dim") = dims;
		// attach column names
		Rcpp::List dimnms = Rcpp::List::create(R_NilValue, ch_selected);
		mat.attr("dimnames") = dimnms;
		//update exprs slot
		fr.slot("exprs") = mat;

	    }
	return(fr);

}



