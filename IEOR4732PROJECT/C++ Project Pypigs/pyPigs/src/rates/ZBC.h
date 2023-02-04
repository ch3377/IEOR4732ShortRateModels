//
// Created by 黄晨 on 2019-10-30.
//

#ifndef PYPIGS_ZBC_H
#define PYPIGS_ZBC_H

#include "Vasicek.h"
#include "IShortRatesFactory.h"
#include "LNLMM.h"
#include <cmath>
#include <math.h>
#include <functional>
namespace rates
{
    double ZBCVasicek(double t, double T, double s, double K, double BtT, double Bts, double kappa, double theta, double sigma,bool call=true);
    double ZBCCIR(double t, double T, double s, double K, double BtT, double Bts, double k, double theta, double sigma,bool call=true);
    double ZBCLFM(double t, double T, double s, double K, double BtT, double Bts, double lambda,bool call=true);
    double ZBCLongstaff(double t, double T, double s, double K, double BtT, double Bts, double r0, double k, double sigma, bool call=true);
}



#endif //PYPIGS_ZBC_H
