//
// Created by 黄晨 on 2019-10-30.
//

#include "ZBC.h"
#include <random>
#include <vector>
#include <numeric>
namespace rates
{
    namespace
    {
        double dOne(double s, double t, double r, double sig, double K){
            return (1.0/(sig*sqrt(t)))*log(s/(K*exp(-r*t)))+sqrt(t)*sig/2;
        }

        double dZero(double s, double t, double r, double sig, double K) {
            return dOne(s,t,r,sig,K)-sqrt(t)*sig;
        }

        double normalCDF(double x){
            return erfc(-x/sqrt(2))/2;
        }

        double ZBCSigPricer(double t, double T, double s, double K, double BtT, double Bts, const function<double(double, double)> &sig,bool call=true){
            return call? Bts*normalCDF(dOne(Bts/BtT,T-t,0,sig(t,T),K))-K*BtT*normalCDF(dZero(Bts/BtT,T-t,0,sig(t,T),K))
                        :-Bts*normalCDF(-dOne(Bts/BtT,T-t,0,sig(t,T),K))+K*BtT*normalCDF(-dZero(Bts/BtT,T-t,0,sig(t,T),K));
        }

        double ZBCSimulatePricer(double t,double T,double s,double K, double BtT, double Bts,int steps, int N, const IShortRatesPtr& ptr,bool call = true){
            vector<double> norms;
            vector<double> values;
            values.reserve(N);
            norms.reserve(steps*N);
            default_random_engine generator;
            normal_distribution<double> dist(0.0,1.0);
            for(auto i=0;i<steps*N;++i)
                norms.push_back(dist(generator));
            double stepSize = (T-t)/steps;

            double rZero = ptr->shortRate(t,s,Bts);
            for(auto i=0;i<N;++i){
                double rLocal = rZero;
                for(auto j=i*steps;j<(i+1)*steps;++j){
                    rLocal=ptr->shortNext(rLocal,stepSize,norms[j]);
                }
                double bondT = ptr->bondPrice(T,s,rLocal);
                double value=0;
                if (call && bondT>K)
                    value = bondT-K;
                else if((!call) && bondT<K)
                    value=K-bondT;
                values.push_back(value);
            }

            return BtT*accumulate(values.begin(),values.end(),0.0)/N;
        }

        double ZBCLFMPricer(double t, double T, double s, double K, double BtT, double Bts, const function<double(double,double)> &gamma,bool call) {
            double forward = Bts/BtT;
            double lambda = gamma(t,T);
            if (call)
                return Bts*(1-K)*normalCDF(dOne(forward*(1-K),T-t,0,lambda,K*(1-forward)))-K*(BtT-Bts)*normalCDF(dZero(forward*(1-K),T-t,0,lambda,K*(1-forward)));
            else
                return -Bts*(1-K)*normalCDF(-dOne(forward*(1-K),T-t,0,lambda,K*(1-forward)))+K*(BtT-Bts)*normalCDF(-dZero(forward*(1-K),T-t,0,lambda,K*(1-forward)));
        }
    }

    double ZBCVasicek(double t, double T, double s, double K, double BtT, double Bts, double kappa, double theta,
                      double sigma,bool call) {
        auto vasicek = IShortRatesFactory::createVasicek(kappa,theta,sigma);
        auto sigFunc = [&](double t1, double t2){return vasicek->bondSigma(t1,t2);};
        return ZBCSigPricer(t,T,s,K,BtT,Bts,sigFunc,call);
    }

    double
    ZBCCIR(double t, double T, double s, double K, double BtT, double Bts, double k, double theta, double sigma,bool call) {

        auto cir = IShortRatesFactory::createCIR(k,theta,sigma);
        return ZBCSimulatePricer(t,T,s,K,BtT,Bts,32,1000,cir,call);
    }

    double ZBCLFM(double t, double T, double s, double K, double BtT, double Bts, double lambda,bool call) {
        auto gamma = [&](double t, double T){return gammaSimple(t,T,lambda);};
        return ZBCLFMPricer(t,T,s,K,BtT,Bts,gamma,call);
    }

    double
    ZBCLongstaff(double t, double T, double s, double K, double BtT, double Bts, double r0, double k, double sigma,
                 bool call) {
        auto longstaff = IShortRatesFactory::createLongstaff(r0,k,sigma);
        return ZBCSimulatePricer(t,T,s,K,BtT,Bts,32,1000,longstaff,call);
    }
}

