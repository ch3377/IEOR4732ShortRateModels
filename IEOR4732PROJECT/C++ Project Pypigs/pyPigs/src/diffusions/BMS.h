//
// Created by 黄晨 on 2019-09-11.
//

#ifndef PYPIGS_BMS_H
#define PYPIGS_BMS_H

#include "IDiff.h"
namespace diffusions
{
    class BMS : public IDiff{
        double S,r,q,sigma;
    public:
        BMS(double S_, double r_, double q_, double sigma_): S(S_),r(r_),q(q_),sigma(sigma_){};
        virtual pigComp characterFunc(const pigComp &v, const double &T) override;
        virtual double getSpot() override;
        virtual double getRiskFree() override;
    };
}
#endif //PYPIGS_BMS_H
