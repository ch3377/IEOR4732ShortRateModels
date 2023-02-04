//
// Created by 黄晨 on 2019-10-30.
//

#ifndef PYPIGS_ISHORTRATES_H
#define PYPIGS_ISHORTRATES_H

#include <memory>


namespace rates
{
    using namespace std;
    class IShortRates;
    using IShortRatesPtr = shared_ptr<IShortRates>;

    class IShortRates{
    public:
        virtual ~IShortRates()=default;
        virtual double bondVol(double t, double T) = 0;//for explicit formulas
        virtual double bondSigma(double t, double T) = 0;//for explicit formulas
        virtual double bondPrice(double t, double T, double r) = 0;//for simulation
        virtual double shortRate(double t, double T, double bondP) = 0;//for simulation
        virtual double shortNext(double rNow,double step,double randomNorm) = 0;//for simulation


    };

}

#endif //PYPIGS_ISHORTRATES_H
