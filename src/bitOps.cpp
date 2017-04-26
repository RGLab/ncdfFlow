#include <Rcpp.h>
//[[Rcpp::export]]
Rcpp::RawVector toBitVec(Rcpp::LogicalVector indx) {
	unsigned nBit = indx.size();
	unsigned nByte = ceil(float(nBit)/8);
	Rcpp::RawVector bytes(nByte);//default are all 0s
	bytes.attr("bitlen") = nBit;

    unsigned byteIndex, bitIndex ;
    for(unsigned i = 0 ; i < nBit; i++) {
        byteIndex = i / 8;
        bitIndex = i % 8;
        if(indx(i) == 1)
			bytes[byteIndex] = bytes[byteIndex] | 1 << bitIndex;
    }
    return bytes;
}

//[[Rcpp::export]]
Rcpp::LogicalVector toLogical(Rcpp::RawVector bytes) {

	unsigned nBit = bytes.attr("bitlen");
//	unsigned nByte = bytes.size();
    Rcpp::LogicalVector ans(nBit);//default are all 0s

    unsigned byteIndex, bitIndex;

    for(unsigned i =0 ; i < nBit; i++) {
        byteIndex = i / 8;
        bitIndex = i % 8;
        if(bytes[byteIndex] & 1 << bitIndex)
        	ans(i) = 1;
    }

    return(ans);
}

