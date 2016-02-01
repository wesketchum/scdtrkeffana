#ifndef TRKEFF_LINEARLEASTSQUARES_CXX
#define TRKEFF_LINEARLEASTSQUARES_CXX

#include "LinearLeastSquaresFit.hh"

#include "RecoBase/Hit.h"
#include "SimpleTypesAndConstants/geo_types.h"

void trkeff::LinearLeastSquaresFit::Clear(){
  _sum_x  = 0;
  _sum_y  = 0;
  _sum_xy = 0;
  _sum_x2 = 0;
  _n      = 0;
}


trkeff::LinearLeastSquaresFit::LeastSquaresResult_t
 trkeff::LinearLeastSquaresFit::LinearFit(std::vector<recob::Hit> const& hit_collection,
					  std::vector<size_t> const& hit_indices)
 {

  LeastSquaresResult_t result;
  if(hit_indices.size()<=1) return result;
  
  Clear();  
  for(auto const& i_h : hit_indices){
    _n      += 1.0;
    _sum_x  += (double)hit_collection[i_h].WireID().Wire;
    _sum_y  += hit_collection[i_h].PeakTime(); //detprop?
    _sum_xy += (double)hit_collection[i_h].WireID().Wire * hit_collection[i_h].PeakTime();
    _sum_x2 += (double)hit_collection[i_h].WireID().Wire * (double)hit_collection[i_h].WireID().Wire;
  }

  result.intercept =  (_sum_y*_sum_x2 - _sum_x*_sum_xy)/(_n*_sum_x2 - _sum_x*_sum_x);
  result.slope     = (_n*_sum_xy - _sum_x*_sum_y)/(_n*_sum_x2 - _sum_x*_sum_x);

  result.chi2 = 0;
  double tmp;
  for(auto const& i_h : hit_indices){
    tmp = result.slope*(double)hit_collection[i_h].WireID().Wire + result.intercept;
    result.chi2 += (hit_collection[i_h].PeakTime() - tmp)*(hit_collection[i_h].PeakTime() - tmp)/tmp;
  }

  result.bad_result = false;
  return result;
  
}

#endif
