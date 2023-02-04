//
// Created by 黄晨 on 2019-10-30.
//

#include "IShortRatesFactory.h"
namespace rates
{

    IShortRatesPtr IShortRatesFactory::createVasicek(double kappa, double theta, double sigma) {
        return make_shared<Vasicek>(kappa,theta,sigma);
    }

    IShortRatesPtr IShortRatesFactory::createCIR(double k, double theta, double sigma) {
        return  make_shared<CIR>(k,theta,sigma);
    }

    IShortRatesPtr IShortRatesFactory::createLongstaff(double r0, double k, double sigma) {
        return make_shared<Longstaff>(r0,k,sigma);
    }
}
