//
// Created by 黄晨 on 2019-10-30.
//

#ifndef PYPIGS_ISHORTRATESFACTORY_H
#define PYPIGS_ISHORTRATESFACTORY_H

#include "Vasicek.h"
#include "CIR.h"
#include "Longstaff.h"

namespace rates
{
    class IShortRatesFactory{
    public:
        static IShortRatesPtr createVasicek(double kappa, double theta, double sigma);
        static IShortRatesPtr createCIR(double k, double theta, double sigma);
        static IShortRatesPtr createLongstaff(double r0, double k, double sigma);
    };
}


#endif //PYPIGS_ISHORTRATESFACTORY_H
