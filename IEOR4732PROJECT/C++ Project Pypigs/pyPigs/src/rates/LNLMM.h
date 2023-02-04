//
// Created by 黄晨 on 2019-11-20.
//

#ifndef PYPIGS_LNLMM_H
#define PYPIGS_LNLMM_H

/*
 * LNLMM is short for log normal libor market models
 * only consider gamma which is sqrt((1/(T-t))*int_t^T{lambda(s,T)})
 * assuming lambda(s,T) is the volatility function for maturity T
 */
namespace rates
{
double gammaSimple(double t, double T,double lambda);//lambda(t,T) is constant
}


#endif //PYPIGS_LNLMM_H
