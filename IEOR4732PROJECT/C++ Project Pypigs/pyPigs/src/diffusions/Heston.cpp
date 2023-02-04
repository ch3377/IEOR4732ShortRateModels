//
// Created by 黄晨 on 2019-09-19.
//

#include "Heston.h"
namespace diffusions
{

    pigComp Heston::characterFunc(const pigComp &v, const double &T) {

        pigComp gam = sqrt(sigma*sigma*(v*v+pigI*v)+(k-pigI*rho*sigma*v)*(k-pigI*rho*sigma*v));
        pigComp partA = exp(pigI*v*log(S)+pigI*v*(r-q)*T+k*theta*T*(k-pigI*rho*sigma*v)/(sigma*sigma));
        pigComp partB = pow(cosh(gam*T/2.0)+sinh(gam*T/2.0)*(k-pigI*rho*sigma*v)/gam,2*k*theta/(sigma*sigma));
        pigComp partC = exp(-(v*v+pigI*v)*vmodel/(gam/tanh(gam*T/2.0))+k-pigI*rho*sigma*v);

        return partA*partC/partB;
    }

    double Heston::getSpot() {
        return S;
    }

    double Heston::getRiskFree() {
        return r;
    }


}
