// #include<vector>
// #include "Rcpp.h"
// #include<vector>
// #include<numeric>
// #include <R.h>
// #include "Rmath.h"
#include<Rcpp.h>
#include<Rinternals.h>
#include "filterIntervals.h"
#include "rcppdatawrap.h"

// using std::vector;
// using std::size_t;


//extract the individual vectors into a dataframe
Rcpp::DataFrame extractDataFrameFromIntervalVector(const vector<Interval>& v){
    //create the vectors
    vector<int> start(v.size());
    vector<int> end(v.size());
    vector<double> pvalue(v.size());
    
    for(unsigned int i = 0; i < v.size(); ++i){
        start[i] = v[i].getStart();
        end[i] = v[i].getEnd();
        pvalue[i] = v[i].getPvalue();
    }

    //create the data frame from the separate vectors
    Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::Named("start")=start,
                                                 Rcpp::Named("end")=end,
                                                 Rcpp::Named("pvalue")=pvalue);

    //return the data frame
    return(df);
}




//--------------------------------------------------------------------------//
//create an empty data frame, in case we need to return it
Rcpp::DataFrame createEmptyDataFrame(){
    //create the vectors
    vector<int> start;
    vector<int> end;
    vector<double> pvalue;
    
    //clear the vectors
    start.clear();
    end.clear();
    pvalue.clear();

    Rcpp::DataFrame emptyDf = Rcpp::DataFrame::create(Rcpp::Named("start")=start,
                                                      Rcpp::Named("end")=end,
                                                      Rcpp::Named("pvalue")=pvalue);
    return(emptyDf);
} 





//extract the individual vectors into a dataframe
Rcpp::DataFrame createDataFrameTauLPvalue(const vector<long long>& tau, const vector<long long>& l, const vector<double>& pvalue){

    vector<int> tauint(tau.begin(), tau.end());
    vector<int> lint(l.begin(), l.end());


    //create the data frame from the separate vectors
    Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::Named("tau")=tauint,
                                                 Rcpp::Named("l")=lint,
                                                 Rcpp::Named("pvalue")=pvalue);

    //return the data frame
    return(df);
}


//--------------------------------------------------------------------------//
//create a return list that indicates an error occurred
Rcpp::DataFrame createErrorReturnList(){
    Rcpp::DataFrame errorDf = Rcpp::DataFrame::create(Rcpp::Named("message")="An error occurred while runnig FastCMH - no output. An error message should have been displayed, and the error probably occurred while reading in the input");
    return(errorDf);
} 
