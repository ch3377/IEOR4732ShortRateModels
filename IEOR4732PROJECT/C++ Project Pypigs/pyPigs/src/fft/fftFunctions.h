//
// Created by 黄晨 on 2019-09-11.
//

#ifndef PYPIGS_FFTFUNCTIONS_H
#define PYPIGS_FFTFUNCTIONS_H

#include <IDiffFactory.h>
#include <map>
#include <vector>

using namespace diffusions;
namespace fft
{
// FFT for pricing European Call given BMS model
    map<double, double>
    myfftBMS(const int N, double S, double ita, double alpha, double r, double q, double sigma, double T, double lower,
             double higher);

    map<double, double>
    myfftHeston(const int N,double S, double ita, double alpha, double r, double q, double sigma, double T, double k, double theta,
            double lamb, double rho, double v, double lower, double higher);

    map<double, double>
    myfftVG(const int N, double S, double ita, double alpha, double r, double q, double sigma, double T, double v, double theta,
            double lower, double higher);

    map<double, double>
    myfrfftBMS(const int N, double S, double ita, double lambd, double alpha, double r, double q, double sigma, double T, double lower,
             double higher);

    map<double, double>
    myfrfftHeston(const int N,double S, double ita, double lambd, double alpha, double r, double q, double sigma, double T, double k, double theta,
                double lamb, double rho, double v, double lower, double higher);

    map<double, double>
    myfrfftVG(const int N, double S, double ita, double lambd, double alpha, double r, double q, double sigma, double T, double v, double theta,
            double lower, double higher);
}
#endif //PYPIGS_FFTFUNCTIONS_H
