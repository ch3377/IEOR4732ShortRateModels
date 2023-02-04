//
// Created by 黄晨 on 2019-09-11.
//

#ifndef PYPIGS_IDIFF_H
#define PYPIGS_IDIFF_H

#include <memory>
#include <complex>
#include <math.h>

// diffusions base class and some complex settings, pigI = i
namespace diffusions
{
    using namespace std;

    class IDiff;

    using IDiffContPtr = shared_ptr<IDiff>;
    using pigComp = complex<double>;
    const pigComp pigI = pigComp(0,1);

    class IDiff{
    public:
        virtual ~IDiff()=default;
        virtual pigComp characterFunc(const pigComp &v, const double &T) = 0;
        virtual double getSpot() = 0;
        virtual double getRiskFree() = 0;

    };
}

#endif //PYPIGS_IDIFF_H
