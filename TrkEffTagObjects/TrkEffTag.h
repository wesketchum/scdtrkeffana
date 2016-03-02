/** ****************************************************************************
 * @file TrkEffTag.h
 * @brief Definition of basic track efficiency tagging object
 * @author wketchum@fnal.gov
 * 
 * ****************************************************************************/

#ifndef TRKEFF_TRKEFFTAG_H
#define TRKEFF_TRKEFFTAG_H

#include "SimpleTypesAndConstants/geo_types.h"

#include <stdint.h>
#include <vector>
#include <algorithm>

#ifndef __GCCXML__
#include <array>
#endif

namespace trkeff {

  class TrkEffTag{

  public:
  TrkEffTag():
    fChi2(-999),fTheta(-999),fPhi(-999){}
    
#ifndef __GCCXML__
  public:

    typedef struct TrkEffTag_Tree {
      
      TrkEffTag_Tree(){};
      
    TrkEffTag_Tree(TrkEffTag const& tag) :
      x_start(tag.StartPoint()[0]), y_start(tag.StartPoint()[1]), z_start(tag.StartPoint()[2]),
	x_end(tag.EndPoint()[0]), y_end(tag.EndPoint()[1]), z_end(tag.EndPoint()[2]),
	chi2(tag.Chi2()), theta(tag.Theta()), thetaYZ(tag.ThetaYZ()), thetaXZ(tag.ThetaXZ()), phi(tag.Phi()) {}
      
      double x_start;
      double y_start;
      double z_start;
      double x_end;
      double y_end;
      double z_end;
      double chi2;
      double theta;
      double thetaYZ;
      double thetaXZ;
      double phi;
    } TrkEffTag_Tree_t;
    

    TrkEffTag(std::array<double,3> const& start, std::array<double,3> const& end,
	      double const& chi2, std::vector<geo::WireID> const& wires)
      {
	for(size_t i=0; i<3; i++){
	  fStartPoint[i] = start[i];
	  fEndPoint[i] = end[i];
	}
	//fStartPoint = start;
	//fEndPoint   = end;
	fChi2       = chi2;
	fWires      = wires;
	std::sort(fWires.begin(),fWires.end());
	std::unique(fWires.begin(),fWires.end());
	
	CalculateAngles();
      }    
    
    TrkEffTag(double const start[], double const end[],
	      double const& chi2, std::vector<geo::WireID> const& wires)
      {
	for(size_t i=0; i<3; i++){
	  fStartPoint[i] = *(start+i);
	  fEndPoint[i] = *(end+i);
	}

	fChi2       = chi2;
	fWires      = wires;
	std::sort(fWires.begin(),fWires.end());
	std::unique(fWires.begin(),fWires.end());

	CalculateAngles();
      }    

    std::array<double,3>     StartPoint()  const;
    std::array<double,3>     EndPoint()    const;
    double                   Chi2()       const;
    std::vector<geo::WireID> Wires()      const;
    double                   Theta()      const;
    double                   ThetaYZ()    const;
    double                   ThetaXZ()    const;
    double                   Phi()        const;

    TrkEffTag_Tree_t GetRootTreeType() const;
    
#endif // !__GCCXML__
  private:
    double fStartPoint[3];
    double fEndPoint[3];
    double fChi2;

    std::vector<geo::WireID> fWires;
    
    double fTheta;
    double fThetaYZ;
    double fThetaXZ;
    double fPhi;

    void CalculateAngles();
    
  }; // class TrkEffTag()
  
#ifndef __GCCXML__
  inline std::array<double,3>  trkeff::TrkEffTag::StartPoint() const
    { return std::array<double,3>{fStartPoint[0],fStartPoint[1],fStartPoint[2]}; }
  inline std::array<double,3>  trkeff::TrkEffTag::EndPoint()   const
    { return std::array<double,3>{fEndPoint[0],fEndPoint[1],fEndPoint[2]}; }
  
  inline double                   trkeff::TrkEffTag::Chi2()       const { return fChi2;       }
  inline std::vector<geo::WireID> trkeff::TrkEffTag::Wires()      const { return fWires;      }
  inline double                   trkeff::TrkEffTag::Theta()      const { return fTheta;      }
  inline double                   trkeff::TrkEffTag::ThetaYZ()    const { return fThetaYZ;    }
  inline double                   trkeff::TrkEffTag::ThetaXZ()    const { return fThetaXZ;    }
  inline double                   trkeff::TrkEffTag::Phi()        const { return fPhi;        }

  inline trkeff::TrkEffTag::TrkEffTag_Tree_t trkeff::TrkEffTag::GetRootTreeType() const
  { return TrkEffTag_Tree_t(*this); }
#endif // !__GCCXML__
  
} // namespace trkeff


#endif // TRKEFF_TRKEFFTAG_H

////////////////////////////////////////////////////////////////////////
