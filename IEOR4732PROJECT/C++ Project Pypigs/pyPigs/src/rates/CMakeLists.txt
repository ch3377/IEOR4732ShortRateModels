project(rates)
set(rates_headers IShortRates.h Vasicek.h IShortRatesFactory.h ZBC.h CIR.h LNLMM.h Longstaff.h)
set(rates_sources Vasicek.cpp IShortRatesFactory.cpp ZBC.cpp CIR.cpp LNLMM.cpp Longstaff.cpp)
add_library(rates STATIC ${rates_headers} ${rates_sources})