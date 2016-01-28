#include "TrkEffTagObjects/TrkEffTag.h"

#include <math.h>

void trkeff::TrkEffTag::CalculateAngles(){

  if(fStartPoint[0]==fEndPoint[0] &&
     fStartPoint[1]==fEndPoint[1])
    fPhi=0;				      
  else
    fPhi = std::atan2(fEndPoint[1]-fStartPoint[1],
		      fEndPoint[0]-fStartPoint[0]);

  if(fStartPoint[0]==fEndPoint[0] &&
     fStartPoint[1]==fEndPoint[1] &&
     fStartPoint[2]==fEndPoint[2])
    fTheta=0;
  else
    fTheta = std::atan2(std::sqrt((fEndPoint[1]-fStartPoint[1])*(fEndPoint[1]-fStartPoint[1]) +
				  (fEndPoint[0]-fStartPoint[0])*(fEndPoint[0]-fStartPoint[0])),
			fEndPoint[2]-fStartPoint[2]);

  if(fStartPoint[1]==fEndPoint[1] &&
     fStartPoint[2]==fEndPoint[2])
    fThetaYZ=0;				      
  else
    fThetaYZ = std::atan2(fEndPoint[1]-fStartPoint[1],
			  fEndPoint[2]-fStartPoint[2]);

  if(fStartPoint[0]==fEndPoint[0] &&
     fStartPoint[2]==fEndPoint[2])
    fThetaXZ=0;				      
  else
    fThetaXZ = std::atan2(fEndPoint[2]-fStartPoint[2],
			  fEndPoint[0]-fStartPoint[0]);

}
