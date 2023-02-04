//
// Created by 黄晨 on 2019-11-23.
//

#ifndef PYPIGS_LONGSTAFF_H
#define PYPIGS_LONGSTAFF_H

#include "IShortRates.h"

namespace rates
{
    class Longstaff : public IShortRates{
        double r0,k,sigma;
        double alpha(double t, double T);
        double beta(double t, double T);
        double gamma(double t, double T);
    public:
        Longstaff(double r0_, double k_, double sigma_): r0(r0_),k(k_),sigma(sigma_){};
        virtual double bondPrice(double t, double T, double r) override;
        virtual double shortRate(double t, double T, double bondP) override;
        virtual double shortNext(double rNow,double step,double randomNorm) override;

        virtual double bondVol(double t, double T) override;
        virtual double bondSigma(double t, double T) override;
    };
}


#endif //PYPIGS_LONGSTAFF_H
