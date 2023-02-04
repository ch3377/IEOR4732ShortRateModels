//
// Created by 黄晨 on 2019-10-30.
//

#include "Vasicek.h"
#include <math.h>
namespace rates
{

    double Vasicek::bondVol(double t, double T) {
        return -sig*(1.0-exp(-kappa*(T-t)))/kappa;
    }

    double Vasicek::bondSigma(double t, double T) {
        return (sig/sqrt(T-t))*sqrt((1.0-exp(-2.0*kappa*(T-t)))/(2.0*kappa))*((1.0-exp(-kappa*(T-t)))/kappa);
    }

    double Vasicek::bondPrice(double t, double T, double r) {
        //No need for explicit use of Vasicek
        return 0;
    }

    double Vasicek::shortRate(double t, double T, double bondP) {
        //No need for explicit use
        return 0;
    }

    double Vasicek::shortNext(double rNow, double step, double randomNorm) {
        //No need for explicit use
        return 0;
    }
}