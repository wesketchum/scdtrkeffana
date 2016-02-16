#ifndef TRKEFF_LINEARLEASTSQUARES_CXX
#define TRKEFF_LINEARLEASTSQUARES_CXX

#include "LinearLeastSquaresFit.hh"

#include "RecoBase/Hit.h"
#include <limits>

void trkeff::LinearLeastSquaresFit::Clear(){
  _sum_x  = 0;
  _sum_y  = 0;
  _sum_xy = 0;
  _sum_x2 = 0;
  _n      = 0;
}


trkeff::LinearLeastSquaresFit::LeastSquaresResult_t
 trkeff::LinearLeastSquaresFit::LinearFit(std::vector<recob::Hit> const& hit_collection,
					  std::vector<size_t> const& hit_indices,
					  bool invert)
 {

  LeastSquaresResult_t result;
  if(hit_indices.size()<=1) return result;
  
  Clear();  
  for(auto const& i_h : hit_indices){
    _n      += 1.0;
    if(!invert){
      _sum_x  += (double)hit_collection[i_h].WireID().Wire;
      _sum_y  += hit_collection[i_h].PeakTime(); //detprop?
      _sum_x2 += (double)hit_collection[i_h].WireID().Wire * (double)hit_collection[i_h].WireID().Wire;
    }
    else{
      _sum_y  += (double)hit_collection[i_h].WireID().Wire;
      _sum_x  += hit_collection[i_h].PeakTime(); //detprop?
      _sum_x2 += (double)hit_collection[i_h].PeakTime() * (double)hit_collection[i_h].PeakTime();
    }
    _sum_xy += (double)hit_collection[i_h].WireID().Wire * hit_collection[i_h].PeakTime();
  }

  result.intercept =  (_sum_y*_sum_x2 - _sum_x*_sum_xy)/(_n*_sum_x2 - _sum_x*_sum_x);
  result.slope     = (_n*_sum_xy - _sum_x*_sum_y)/(_n*_sum_x2 - _sum_x*_sum_x);

  result.chi2 = 0;
  result.npts = hit_indices.size();
  result.max_outlier_value = 0;

  double min_time = std::numeric_limits<double>::max();
  double max_time = std::numeric_limits<double>::min();
  double tmp,diff;
  for(auto const& i_h : hit_indices){

    if(hit_collection[i_h].PeakTime() < min_time){
      min_time = hit_collection[i_h].PeakTime();
      result.min_wire = hit_collection[i_h].WireID();      
    }
    if(hit_collection[i_h].PeakTime() > max_time){
      max_time = hit_collection[i_h].PeakTime();
      result.max_wire = hit_collection[i_h].WireID();      
    }
    
    if(!invert){
      tmp = result.slope*(double)(hit_collection[i_h].WireID().Wire) + result.intercept;
      diff = (hit_collection[i_h].PeakTime() - tmp) / hit_collection[i_h].RMS();
    }
    else{
      tmp = result.slope*(double)hit_collection[i_h].PeakTime() + result.intercept;
      diff = (double)(hit_collection[i_h].WireID().Wire - tmp) / 0.5;
    }
    result.chi2 += diff*diff;
    
    if( std::abs(diff) > std::abs(result.max_outlier_value) ){
      result.max_outlier_value = diff;
      result.max_outlier_index = i_h;
    }
  }
  
  result.bad_result = false;
  return result;
  
}

#endif
