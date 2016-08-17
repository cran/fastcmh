#ifndef FILTERINTERVALS_H
#define FILTERINTERVALS_H

//filterIntervals.h
//      read csv file containing all significant intervals, group these into
//      overlapping clusters, and return the most significant interval per cluster
//
//Dean Bodenham June 2016

#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<iomanip>
#include<numeric>

using std::vector;
using std::size_t;
using std::string;
using std::stringstream;
using std::cout;

//--------------------------------------------------------------------------//
//--------------------------------------------------------------------------//
class Interval{
    public:
        size_t getStart() const;
        size_t getEnd() const;
        double getPvalue() const;
        size_t getLength() const;
        void setStart(size_t);
        void setEnd(size_t);
        void setEnd(size_t, size_t);
        void setPvalue(double);
        bool overlaps(size_t, size_t) const;
        void printInterval() const;
    private:
        size_t start;
        size_t end;
        double pvalue;
};
//--------------------------------------------------------------------------//

vector<Interval> cpp_filterIntervalsFromMemory(vector<long long>,
                                              vector<long long>,
                                              vector<double> );
//--------------------------------------------------------------------------//
#endif
