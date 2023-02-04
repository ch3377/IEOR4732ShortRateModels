//
// Created by 黄晨 on 2019-11-23.
//

#include "Longstaff.h"
#include <math.h>

namespace rates
{

    double Longstaff::alpha(double t, double T) {
        auto minus = 1 - exp(sigma*sqrt(2)*(T-t));
        auto plus = 1 + exp(sigma*sqrt(2)*(T-t));
        return 0.5*log(2/plus) + ((T-t)/4)*(sigma*sqrt(2) - 2*k*k/(sigma*sigma))-(k*k/(sqrt(2)*sigma*sigma*sigma))*(minus/plus);
    }

    double Longstaff::beta(double t, double T) {
        auto minus = 1 - exp(sigma*(T-t)/sqrt(2));
        auto plus = 1 + exp(sigma*sqrt(2)*(T-t));
        return (2*k/(sigma*sigma))*minus*minus/plus;
    }

    double Longstaff::gamma(double t, double T) {
        auto minus = 1 - exp(sigma*sqrt(2)*(T-t));
        auto plus = 1 + exp(sigma*sqrt(2)*(T-t));
        return (sqrt(2)/sigma)*minus/plus;
    }

    double Longstaff::bondPrice(double t, double T, double r) {
        return exp(alpha(t,T)+beta(t,T)*sqrt(r)+gamma(t,T)*r);
    }

    double Longstaff::shortRate(double t, double T, double bondP) {
        return r0;
    }

    double Longstaff::shortNext(double rNow, double step, double randomNorm) {
        return rNow>0? rNow + k*(sigma*sigma/(4*k)-sqrt(rNow))*step+sigma*sqrt(rNow)*sqrt(step)*randomNorm + (sigma*sigma/4)*step*(randomNorm*randomNorm-1)
                        :rNow + (sigma*sigma/4)*step;
    }

    double Longstaff::bondVol(double t, double T) {
        return 0;
    }

    double Longstaff::bondSigma(double t, double T) {
        return 0;
    }
}