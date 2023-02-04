//
// Created by 黄晨 on 2019-09-19.
//

#ifndef PYPIGS_VG_H
#define PYPIGS_VG_H


#include "IDiff.h"
namespace diffusions
{
    class VG : public IDiff{
        double S,r,q,sigma,vmodel,theta;
    public:
        VG(double S_,double r_,double q_,double sigma_,double v_, double theta_): S(S_),r(r_),q(q_),sigma(sigma_),vmodel(v_),theta(theta_){};
        virtual pigComp characterFunc(const pigComp &v, const double &T) override;
        virtual double getSpot() override;
        virtual double getRiskFree() override;
    };
}


#endif //PYPIGS_VG_H
