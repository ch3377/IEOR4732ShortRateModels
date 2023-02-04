//
// Created by 黄晨 on 2019-09-19.
//

#include "VG.h"
namespace diffusions
{

    pigComp VG::characterFunc(const pigComp &v, const double &T) {

        pigComp partA = exp(pigI*v*(log(S)+(r-q)*T));
        pigComp partB = pow(1.0/(1.0-pigI*v*theta*vmodel+sigma*sigma*v*v*vmodel/2.0),T/v);

        return partA*partB;
    }

    double VG::getSpot() {
        return S;
    }

    double VG::getRiskFree() {
        return r;
    }
}
