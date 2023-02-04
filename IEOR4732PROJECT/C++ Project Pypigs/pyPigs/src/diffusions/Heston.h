//
// Created by 黄晨 on 2019-09-19.
//

#ifndef PYPIGS_HESTON_H
#define PYPIGS_HESTON_H

#include "IDiff.h"
namespace diffusions
{
    class Heston : public IDiff{
        double S,r,q,sigma,k,theta,lambda,rho,vmodel;
    public:
        Heston(double S_,double r_,double q_,double sigma_,double k_,double theta_,double lambda_,double rho_,double v_):
        S(S_),r(r_),q(q_),sigma(sigma_),k(k_),theta(theta_),lambda(lambda_),rho(rho_),vmodel(v_){};
        virtual pigComp characterFunc(const pigComp &v, const double &T) override;
        virtual double getSpot() override;
        virtual double getRiskFree() override;
    };
}

#endif //PYPIGS_HESTON_H
