#ifndef TRKEFF_LINEARLEASTSQUARES_H
#define TRKEFF_LINEARLEASTSQUARES_H

#include <vector>

#include "SimpleTypesAndConstants/geo_types.h"

namespace recob{ class Hit; }

namespace trkeff{
  class LinearLeastSquaresFit;
}

class trkeff::LinearLeastSquaresFit{
  
public:

  typedef struct LeastSquaresResult{
    double slope;
    double intercept;
    double chi2;
    unsigned int npts;
    std::size_t max_outlier_index;
    double max_outlier_value;

    geo::WireID min_wire;
    geo::WireID max_wire;
    
    bool   bad_result;
    LeastSquaresResult():slope(0),intercept(0),chi2(-999),npts(0),bad_result(true){}
  } LeastSquaresResult_t;
  
  /// Default constructor
  LinearLeastSquaresFit() { Clear(); }

  void Clear();

  LeastSquaresResult_t LinearFit(std::vector<recob::Hit> const&,std::vector<std::size_t> const&,bool invert=false);
  
  /// Default destructor
  virtual ~LinearLeastSquaresFit(){};

 private:

  double _sum_x;
  double _sum_y;
  double _sum_xy;
  double _sum_x2;
  double _n;

  double _slope;
  double _intercept;
  double _chi2;
};

#endif
