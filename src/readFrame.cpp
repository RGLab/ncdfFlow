#include "hdfFlow.h"
typedef std::vector<std::string> strVec;
typedef std::vector<int> intVec;
typedef std::vector<unsigned> uintVec;
typedef std::vector<bool> boolVec;

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
		Rcpp::RObject fspdObj = fspd.slot("data");
		Rcpp::DataFrame fsdata = Rcpp::DataFrame(fspdObj.get__());
		Rcpp::CharacterVector sn = fsdata.attr("row.names");
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
	Rcpp::S4 pheno = fr.slot("parameters");
	Rcpp::RObject pd = pheno.slot("data");
	Rcpp::DataFrame pData = Rcpp::DataFrame(pd.get__());
	//this implicit construction works for g++ and llvm-g++
	//but fails on clang++
//		Rcpp::DataFrame pData =pheno.slot("data");

	Rcpp::CharacterVector pd_rn = pData.attr("row.names");
	Rcpp::CharacterVector pd_name = pData["name"];

	  //get local channel names
	  Rcpp::StringVector colnames = clone(pd_name);
    colnames.attr("class") = R_NilValue;
    // if(colnames.hasAttribute("name"))
      colnames.attr("names") = R_NilValue;
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

	      //open hdf
	      bool is3d;
	      hid_t fileid, dataset, dataspace;
	      open_hdf(file, H5F_ACC_RDONLY, fileid, dataset, dataspace, is3d);

	      //query nrows
	      unsigned nEvents = get_event_number(fileid, dataset, dataspace, samplePos, is3d);
	      //allocate buffer
	      Rcpp::NumericVector mat(nEvents * nCh);
	      if(dataset>0)//make sure the dataset is opened before making readSlice call
	    	readSlice(fileid, dataset, dataspace, chIndx, samplePos, nEvents, mat, is3d);//read data from hdf to buffer
	      close_hdf(fileid);

//	      subset data by indices if necessary
	      if(subByIndice){
	    	  /*
	    	   * convert bytes to bool vector
	    	   */
	    	  Rcpp::RawVector bytes(Indice.get__());
	    	  unsigned len = bytes.attr("bitlen");
	    	  uintVec indx;
			  unsigned byteIndex, bitIndex;
			  for(unsigned i =0 ; i < len; i++)
			  {
				  byteIndex = i / 8;
				  bitIndex = i % 8;
				  if(bytes[byteIndex] & 1 << bitIndex)
					  indx.push_back(i);
			  }

			  //use armadillo to subset mat
			  arma::mat mat1 = arma::mat(mat.begin(), len, nCh, false);
//			  Rcpp::Rcout << indx.size() << std::endl;
			  arma::uvec uindx = Rcpp::as<arma::uvec>(Rcpp::wrap(indx));
			  arma::mat subMat = mat1.rows(uindx);

			  mat = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(subMat));
			  //update nrows
			  nEvents = uindx.size();
	      }


	  //set dims
		Rcpp::IntegerVector dims(2);
		dims[0] = nEvents;
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



