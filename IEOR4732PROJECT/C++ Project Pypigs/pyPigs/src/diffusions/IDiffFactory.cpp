//
// Created by 黄晨 on 2019-09-11.
//

#include "IDiffFactory.h"
namespace diffusions
{
    IDiffContPtr IDiffFactory::createBMS(const double S, const double r, const double q, const double sigma) {
        return make_shared<BMS>(S,r,q,sigma);
    }

    IDiffContPtr
    IDiffFactory::createHeston(const double S, const double r, const double q, const double sigma, const double k,
                               const double theta, const double lambda, const double rho, const double v) {
        return make_shared<Heston>(S,r,q,sigma,k,theta,lambda,rho,v);
    }

    IDiffContPtr
    IDiffFactory::createVG(const double S, const double r, const double q, const double sigma, const double v,
                           const double theta) {
        return make_shared<VG>(S,r,q,sigma,v,theta);
    }
}