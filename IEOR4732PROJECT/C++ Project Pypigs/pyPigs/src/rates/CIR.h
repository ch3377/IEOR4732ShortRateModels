//
// Created by 黄晨 on 2019-11-13.
//

#ifndef PYPIGS_CIR_H
#define PYPIGS_CIR_H


#include "IShortRates.h"

namespace rates
{
    class CIR: public IShortRates{
        double k, theta, sigma;
        double alpha(double t, double T);
        double beta(double t, double T);
    public:
        CIR(double k_, double theta_, double sigma_):k(k_),theta(theta_),sigma(sigma_){};
        virtual double bondPrice(double t, double T, double r) override;
        virtual double shortRate(double t, double T, double bondP) override;
        virtual double shortNext(double rNow,double step,double randomNorm) override;

        virtual double bondVol(double t, double T) override;
        virtual double bondSigma(double t, double T) override;
    };

}

#endif //PYPIGS_CIR_H
