//
// Created by 黄晨 on 2019-11-13.
//

#include "CIR.h"
#include <math.h>

namespace rates
{
    double CIR::alpha(double t, double T) {

        double sq = sqrt(k*k+2*sigma*sigma);
        double N = sq*exp(k*(T-t)/2);
        double D = sq*cosh(sq*(T-t)/2)+k*sinh(sq*(T-t)/2);

        return ((2*k*theta)/(sigma*sigma))*log(N/D);
    }

    double CIR::beta(double t, double T) {
        double sq = sqrt(k*k+2*sigma*sigma);
        double D = sq/tanh(sq*(T-t)/2)+k;
        return 2/D;
    }

    double CIR::bondPrice(double t, double T, double r) {
        return exp(alpha(t,T)-beta(t,T)*r);
    }

    double CIR::shortRate(double t, double T, double bondP) {
        return (alpha(t,T)-log(bondP))/beta(t,T);
    }

    double CIR::shortNext(double rNow, double step, double randomNorm) {
        return rNow>0? rNow+k*(theta-rNow)*step+sigma*sqrt(rNow)*sqrt(step)*randomNorm + (sigma*sigma/4)*step*(randomNorm*randomNorm-1)
                        : rNow+k*(theta-rNow)*step;
    }

    double CIR::bondVol(double t, double T) {
        //Not feasible for stochastic vol bond dynamic
        return 0;
    }

    double CIR::bondSigma(double t, double T) {
        //Not feasible for stochastic vol bond dynamic
        return 0;
    }


}