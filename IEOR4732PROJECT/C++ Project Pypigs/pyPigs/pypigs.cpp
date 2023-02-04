#include <iostream>
#include <IDiffFactory.h>
//#include <fftFunctions.h>
#include <ZBC.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace diffusions;
//using namespace fft;
using namespace rates;
namespace py = pybind11;
PYBIND11_DECLARE_HOLDER_TYPE(T,std::shared_ptr<T>);

// pypigs functions interface
//map<double, double> fftBMS(const int N, double S, double ita, double alpha, double r, double q, double sigma, double T,double lower, double higher)
//{
//    return myfftBMS(N,S,ita,alpha,r,q,sigma,T,lower,higher);
//}
//
//map<double, double> fftHeston(const int N,double S, double ita, double alpha, double r, double q, double sigma, double T, double k, double theta,
//                             double lamb, double rho, double v, double lower, double higher)
//{
//    return myfftHeston(N,S,ita,alpha,r,q,sigma,T,k,theta,lamb,rho,v,lower,higher);
//}
//
//map<double, double> fftVG(const int N, double S, double ita, double alpha,double r, double q,double sigma, double T, double v, double theta,
//        double lower, double higher)
//{
//    return myfftVG(N,S,ita,alpha,r,q,sigma,T,v,theta,lower,higher);
//}
//
//map<double, double> frfftBMS(const int N, double S, double ita, double lambd, double alpha, double r, double q, double sigma, double T,double lower, double higher)
//{
//    return myfrfftBMS(N,S,ita,lambd,alpha,r,q,sigma,T,lower,higher);
//}
//
//map<double, double> frfftHeston(const int N,double S, double ita, double lambd, double alpha, double r, double q, double sigma, double T, double k, double theta,
//                              double lamb, double rho, double v, double lower, double higher)
//{
//    return myfrfftHeston(N,S,ita,lambd,alpha,r,q,sigma,T,k,theta,lamb,rho,v,lower,higher);
//}
//
//map<double, double> frfftVG(const int N, double S, double ita,double lambd, double alpha,double r, double q,double sigma, double T, double v, double theta,
//                          double lower, double higher)
//{
//    return myfrfftVG(N,S,ita,lambd,alpha,r,q,sigma,T,v,theta,lower,higher);
//}

double ZBCVasicekPricer(double t, double T, double s, double K, double BtT, double Bts, double kappa, double theta, double sigma){
    return ZBCVasicek(t,T,s,K,BtT,Bts,kappa,theta,sigma);
}

double ZBPVasicekPricer(double t, double T, double s, double K, double BtT, double Bts, double kappa, double theta, double sigma){
    return ZBCVasicek(t,T,s,K,BtT,Bts,kappa,theta,sigma,false);
}
double ZBCCIRPricer(double t, double T, double s, double K, double BtT, double Bts, double k, double theta, double sigma){
    return ZBCCIR(t,T,s,K,BtT,Bts,k,theta,sigma);
}

double ZBPCIRPricer(double t, double T, double s, double K, double BtT, double Bts, double k, double theta, double sigma){
    return ZBCCIR(t,T,s,K,BtT,Bts,k,theta,sigma,false);
}

double ZBCLFMPricer(double t, double T, double s, double K, double BtT, double Bts, double lambda){
    return ZBCLFM(t,T,s,K,BtT,Bts,lambda);
}

double ZBPLFMPricer(double t, double T, double s, double K, double BtT, double Bts, double lambda){
    return ZBCLFM(t,T,s,K,BtT,Bts,lambda,false);
}

double ZBCLongstaffPricer(double t, double T, double s, double K, double BtT, double Bts, double r0, double k, double sigma){
    return ZBCLongstaff(t,T,s,K,BtT,Bts,r0,k,sigma);
}

double ZBPLongstaffPricer(double t, double T, double s, double K, double BtT, double Bts, double r0, double k, double sigma){
    return ZBCLongstaff(t,T,s,K,BtT,Bts,r0,k,sigma,false);
}

double capletVasicek(double t,double T,double s,double N, double K,double BtT, double Bts, double kappa, double theta, double sigma){
    return N*(1+(s-T)*K)*ZBPVasicekPricer(t,T,s,1.0/(1+(s-T)*K),BtT,Bts,kappa,theta,sigma);
}

double floorletVasicek(double t,double T,double s,double N, double K,double BtT, double Bts, double kappa, double theta, double sigma){
    return N*(1+(s-T)*K)*ZBCVasicekPricer(t,T,s,1.0/(1+(s-T)*K),BtT,Bts,kappa,theta,sigma);
}

double capletCIR(double t,double T,double s,double N, double K,double BtT, double Bts, double k, double theta, double sigma){
    return N*(1+(s-T)*K)*ZBPCIRPricer(t,T,s,1.0/(1+(s-T)*K),BtT,Bts,k,theta,sigma);
}

double floorletCIR(double t,double T,double s,double N, double K,double BtT, double Bts, double k, double theta, double sigma){
    return N*(1+(s-T)*K)*ZBCCIRPricer(t,T,s,1.0/(1+(s-T)*K),BtT,Bts,k,theta,sigma);
}

double capletLFM(double t,double T,double s,double N, double K,double BtT, double Bts, double lambda){
    return N*(1+(s-T)*K)*ZBPLFMPricer(t,T,s,1.0/(1+(s-T)*K),BtT,Bts,lambda);
}

double floorletLFM(double t,double T,double s,double N, double K,double BtT, double Bts, double lambda){
    return N*(1+(s-T)*K)*ZBCLFMPricer(t,T,s,1.0/(1+(s-T)*K),BtT,Bts,lambda);
}

double capletLongstaff(double t,double T,double s,double N,double K,double BtT,double Bts,double r0, double k, double sigma){
    return N*(1+(s-T)*K)*ZBPLongstaffPricer(t,T,s,1.0/(1+(s-T)*K),BtT,Bts,r0,k,sigma);
}

double floorletLongstaff(double t,double T,double s,double N,double K,double BtT,double Bts,double r0, double k, double sigma){
    return N*(1+(s-T)*K)*ZBCLongstaffPricer(t,T,s,1.0/(1+(s-T)*K),BtT,Bts,r0,k,sigma);
}

//set module pypigs
PYBIND11_MODULE(pypigs, m){

    m.doc() = "Baby pricing library, class based w.r.t. IEOR4732 IEOR4710 \n Student: Chen Huang\n Email: ch3377@columbia.edu";

//    m.def("fftBMS", &fftBMS,
//          "Generate European Call option prices with fft methods given BMS model",
//          py::arg("N"),py::arg("S"),py::arg("ita"),py::arg("alpha"),py::arg("r"),py::arg("q"),py::arg("sigma"),py::arg("T"),
//          py::arg("lower"),py::arg("higher"));
//
//    m.def("fftHeston", &fftHeston,
//            "Generate European Call option prices with fft methods given Heston model",
//          py::arg("N"),py::arg("S"),py::arg("ita"),py::arg("alpha"),py::arg("r"),py::arg("q"),py::arg("sigma"),py::arg("T"),
//          py::arg("k"),py::arg("theta"),py::arg("lambda"),py::arg("rho"),py::arg("v0"),
//          py::arg("lower"),py::arg("higher"));
//
//    m.def("fftVG",&fftVG,
//            "Generate European Call option prices with fft methods given VG model",
//          py::arg("N"),py::arg("S"),py::arg("ita"),py::arg("alpha"),py::arg("r"),py::arg("q"),py::arg("sigma"),py::arg("T"),
//          py::arg("v"),py::arg("theta"),
//          py::arg("lower"),py::arg("higher"));
//
//    m.def("frFftBMS", &frfftBMS,
//          "Generate European Call option prices with frFft methods given BMS model",
//          py::arg("N"),py::arg("S"),py::arg("ita"),py::arg("lambda"),py::arg("alpha"),py::arg("r"),py::arg("q"),py::arg("sigma"),py::arg("T"),
//          py::arg("lower"),py::arg("higher"));
//
//    m.def("frFftHeston", &frfftHeston,
//          "Generate European Call option prices with fft methods given Heston model",
//          py::arg("N"),py::arg("S"),py::arg("ita"),py::arg("lambda"),py::arg("alpha"),py::arg("r"),py::arg("q"),py::arg("sigma"),py::arg("T"),
//          py::arg("k"),py::arg("theta"),py::arg("hLambda"),py::arg("rho"),py::arg("v0"),
//          py::arg("lower"),py::arg("higher"));
//
//    m.def("frFftVG",&frfftVG,
//          "Generate European Call option prices with fft methods given VG model",
//          py::arg("N"),py::arg("S"),py::arg("ita"),py::arg("lambda"),py::arg("alpha"),py::arg("r"),py::arg("q"),py::arg("sigma"),py::arg("T"),
//          py::arg("v"),py::arg("theta"),
//          py::arg("lower"),py::arg("higher"));

    m.def("capletVasicek",&capletVasicek,
            "Calculate caplet price under Vasicek model",
            py::arg("t"),py::arg("T"),py::arg("s"),py::arg("N"),py::arg("K"),py::arg("BtT"),py::arg("Bts"),py::arg("kappa"),py::arg("theta"),py::arg("sigma"));

    m.def("floorletVasicek",&floorletVasicek,
          "Calculate floorlet price under Vasicek model",
          py::arg("t"),py::arg("T"),py::arg("s"),py::arg("N"),py::arg("K"),py::arg("BtT"),py::arg("Bts"),py::arg("kappa"),py::arg("theta"),py::arg("sigma"));

    m.def("capletCIR",&capletCIR,
          "Calculate caplet price under CIR model",
          py::arg("t"),py::arg("T"),py::arg("s"),py::arg("N"),py::arg("K"),py::arg("BtT"),py::arg("Bts"),py::arg("k"),py::arg("theta"),py::arg("sigma"));

    m.def("floorletCIR",&floorletCIR,
          "Calculate floorlet price under CIR model",
          py::arg("t"),py::arg("T"),py::arg("s"),py::arg("N"),py::arg("K"),py::arg("BtT"),py::arg("Bts"),py::arg("k"),py::arg("theta"),py::arg("sigma"));

    m.def("capletLFM",&capletLFM,
          "Calculate caplet price under simple LFM model",
          py::arg("t"),py::arg("T"),py::arg("s"),py::arg("N"),py::arg("K"),py::arg("BtT"),py::arg("Bts"),py::arg("lambda"));

    m.def("floorletLFM",&floorletLFM,
          "Calculate floorlet price under simple LFM model",
          py::arg("t"),py::arg("T"),py::arg("s"),py::arg("N"),py::arg("K"),py::arg("BtT"),py::arg("Bts"),py::arg("lambda"));

    m.def("capletLongstaff",&capletLongstaff,
          "Calculate caplet price under Longstaff model",
          py::arg("t"),py::arg("T"),py::arg("s"),py::arg("N"),py::arg("K"),py::arg("BtT"),py::arg("Bts"),py::arg("r0"),py::arg("k"),py::arg("sigma"));

    m.def("floorletLongstaff",&floorletLongstaff,
          "Calculate floorlet price under Longstaff model",
          py::arg("t"),py::arg("T"),py::arg("s"),py::arg("N"),py::arg("K"),py::arg("BtT"),py::arg("Bts"),py::arg("r0"),py::arg("k"),py::arg("sigma"));

}