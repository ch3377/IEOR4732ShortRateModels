//
// Created by 黄晨 on 2019-09-11.
//

#include "BMS.h"
namespace diffusions
{
    pigComp BMS::characterFunc(const diffusions::pigComp &v, const double &T)
    {
        return exp(pigI*(log(S)+(r-q-sigma*sigma/2)*T)*v-(sigma*sigma*v*v)*T/2.0);
    }

    double BMS::getSpot()
    {
       return S;
    }

    double BMS::getRiskFree() {
        return r;
    }


}