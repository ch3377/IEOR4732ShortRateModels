project(diffusions)

set(diff_headers IDiff.h  BMS.h  IDiffFactory.h Heston.h VG.h)
set(diff_sources BMS.cpp IDiffFactory.cpp Heston.cpp VG.cpp)
add_library(diffusions STATIC ${diff_headers} ${diff_sources})