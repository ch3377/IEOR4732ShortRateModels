//
// Created by 黄晨 on 2019-10-30.
//

#ifndef PYPIGS_VASICEK_H
#define PYPIGS_VASICEK_H

#include "IShortRates.h"

namespace rates
{
    class Vasicek : public IShortRates{
        double kappa, theta, sig;
    public:
        Vasicek(double kappa_, double theta_, double sig_):kappa(kappa_),theta(theta_),sig(sig_){};
        virtual double bondVol(double t, double T) override;
        virtual double bondSigma(double t, double T) override;
        virtual double bondPrice(double t, double T, double r) override;
        virtual double shortRate(double t, double T, double bondP) override;
        virtual double shortNext(double rNow,double step,double randomNorm) override;

    };

}






#endif //PYPIGS_VASICEK_H
