//
// Created by 黄晨 on 2019-09-11.
//

#ifndef PYPIGS_IDIFFFACTORY_H
#define PYPIGS_IDIFFFACTORY_H

#include "BMS.h"
#include "Heston.h"
#include "VG.h"

namespace diffusions
{
    class IDiffFactory {
    public:
        static IDiffContPtr createBMS(const double S, const double r, const double q, const double sigma);
        static IDiffContPtr createHeston(const double S, const double r, const double q, const double sigma,
                const double k, const double theta, const double lambda, const double rho, const double v);
        static IDiffContPtr createVG(const double S, const double r, const double q, const double sigma,const double v,const double theta);

    };
}

#endif //PYPIGS_IDIFFFACTORY_H
