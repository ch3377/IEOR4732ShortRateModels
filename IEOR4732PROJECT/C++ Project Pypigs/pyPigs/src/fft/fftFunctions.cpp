//
// Created by 黄晨 on 2019-09-11.
//

#include "fftFunctions.h"
#include <fftw3.h>
namespace fft
{
    namespace
    {
        map<double,double> _fft(const int N, double ita, double alpha,double T, double lower, double higher, IDiffContPtr &model)
        {
            map<double, double> result;
            double lambd = 2 * M_PI / (N * ita);
            double beita = log(model->getSpot()) - lambd * N / 2;

            vector<pigComp> xVec;
            xVec.reserve(N);

            for (auto i = 1U; i <= N; ++i) {
                auto v = (i - 1) * ita;
                auto subx0 = i == 1 ? ita / 2 : ita;
                auto subx1 = exp(-model->getRiskFree() * T) / ((alpha + pigI * v) * (alpha + pigI * v + 1.0));
                auto subx2 = exp(-pigI * beita * v) * model->characterFunc(v - (alpha + 1) * pigI, T);
                xVec.push_back(subx0 * subx1 * subx2);
            }

            fftw_complex x[N];
            fftw_complex y[N];
            for (auto i = 0U; i < N; ++i) {
                x[i][0] = real(xVec[i]);
                x[i][1] = imag(xVec[i]);
            }

            fftw_plan plan = fftw_plan_dft_1d(N, x, y, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);
            fftw_cleanup();

            for (auto i = 0U; i < N; ++i) {
                auto k = beita + i * lambd;
                auto K = exp(k);
                if (K < higher and K > lower)
                    result.insert(make_pair(K, y[i][0] * exp(-alpha * k) / M_PI));
            }
            return result;
        }

        map<double,double> _frFft(const int N, double ita,double lambd, double alpha,double T, double lower, double higher, IDiffContPtr &model)
        {
            map<double, double> result;
            double gam = ita*lambd/(2*M_PI);
            double beita = log(model->getSpot()) - lambd * N / 2;

            vector<pigComp> xVec,yVec,zVec,epVec;
            xVec.reserve(N);
            yVec.reserve(2*N);
            zVec.reserve(2*N);
            epVec.reserve(2*N);

            for (auto i = 1U; i <= N; ++i) {
                auto v = (i - 1) * ita;
                auto subx0 = i == 1 ? ita / 2 : ita;
                auto subx1 = exp(-model->getRiskFree() * T) / ((alpha + pigI * v) * (alpha + pigI * v + 1.0));
                auto subx2 = exp(-pigI * beita * v) * model->characterFunc(v - (alpha + 1) * pigI, T);
                xVec.push_back(subx0 * subx1 * subx2);
            }

            for(auto i = 1U;i<=2*N;++i){
                auto partY = i<=N? xVec[i-1]*exp(-pigI*M_PI*double((i-1)*(i-1))*gam) : 0;
                auto partZ = i<=N? exp(pigI*M_PI*double((i-1)*(i-1))*gam) : exp(pigI*M_PI*double((2*N-i)*(2*N-i))*gam);
                yVec.push_back(partY);
                zVec.push_back(partZ);
            }

            fftw_complex x1[2*N];
            fftw_complex y1[2*N];
            fftw_complex x2[2*N];
            fftw_complex y2[2*N];
            fftw_complex x3[2*N];
            fftw_complex y3[2*N];
            for (auto i = 0U; i < 2*N; ++i) {
                x1[i][0] = real(yVec[i]);
                x1[i][1] = imag(yVec[i]);
                x2[i][0] = real(zVec[i]);
                x2[i][1] = imag(zVec[i]);
            }

            fftw_plan plan1 = fftw_plan_dft_1d(2*N, x1, y1, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_plan plan2 = fftw_plan_dft_1d(2*N, x2, y2, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(plan1);
            fftw_execute(plan2);

            for(auto i = 1U;i<=2*N;++i){
                auto partEp = pigComp(y1[i-1][0],y1[i-1][1])*pigComp(y2[i-1][0],y2[i-1][1]);
                epVec.push_back(partEp);
            }

            for (auto i = 0U; i < 2*N; ++i){
                x3[i][0] = real(epVec[i]);
                x3[i][1] = imag(epVec[i]);
            }

            fftw_plan plan3 = fftw_plan_dft_1d(2*N, x3, y3, FFTW_BACKWARD, FFTW_ESTIMATE);
            fftw_execute(plan3);
            fftw_destroy_plan(plan1);
            fftw_destroy_plan(plan2);
            fftw_destroy_plan(plan3);
            fftw_cleanup();

            for (auto i = 0U; i < N; ++i) {
                auto k = beita + i * lambd;
                auto k_ = (beita+i)*lambd;
                auto K = exp(k);
                if (K < higher and K > lower)
                {
                    auto myReal = real(pigComp(y3[i][0],y3[i][1])*exp(-pigI*M_PI*gam*double(i*i)));
                    result.insert(make_pair(K, myReal * exp(-alpha * k_) / M_PI));
                }

            }
            return result;
        }
    }

    map<double, double>
    myfftBMS(const int N, double S, double ita, double alpha, double r, double q, double sigma, double T, double lower,
             double higher){

        IDiffContPtr testBMS = IDiffFactory::createBMS(S, r, q, sigma);
        return _fft(N,ita,alpha,T,lower,higher,testBMS);
    }

    map<double, double>
    myfftHeston(const int N, double S, double ita, double alpha, double r, double q, double sigma, double T, double k,
                double theta, double lamb, double rho, double v, double lower, double higher) {
        IDiffContPtr testHeston = IDiffFactory::createHeston(S,r,q,sigma,k,theta,lamb,rho,v);
        return _fft(N,ita,alpha,T,lower,higher,testHeston);
    }

    map<double, double>
    myfftVG(const int N, double S, double ita, double alpha, double r, double q, double sigma, double T, double v,
            double theta,double lower, double higher) {
        IDiffContPtr testVG = IDiffFactory::createVG(S,r,q,sigma,v,theta);
        return _fft(N,ita,alpha,T,lower,higher,testVG);
    }

    map<double, double>
    myfrfftBMS(const int N, double S, double ita, double lambd, double alpha, double r, double q, double sigma, double T, double lower,
             double higher){

        IDiffContPtr testBMS = IDiffFactory::createBMS(S, r, q, sigma);
        return _frFft(N,ita,lambd,alpha,T,lower,higher,testBMS);
    }

    map<double, double>
    myfrfftHeston(const int N, double S, double ita, double lambd, double alpha, double r, double q, double sigma, double T, double k,
                double theta, double lamb, double rho, double v, double lower, double higher) {
        IDiffContPtr testHeston = IDiffFactory::createHeston(S,r,q,sigma,k,theta,lamb,rho,v);
        return _frFft(N,ita,lambd,alpha,T,lower,higher,testHeston);
    }

    map<double, double>
    myfrfftVG(const int N, double S, double ita, double lambd,double alpha, double r, double q, double sigma, double T, double v,
            double theta,double lower, double higher) {
        IDiffContPtr testVG = IDiffFactory::createVG(S,r,q,sigma,v,theta);
        return _frFft(N,ita,lambd,alpha,T,lower,higher,testVG);
    }


}