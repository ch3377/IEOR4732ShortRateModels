project(fft)

include_directories(/usr/local/include)
include_directories(${pigs_SOURCE_DIR}/diffusions)

set(fft_headers fftFunctions.h)
set(fft_sources fftFunctions.cpp)
add_library(fft STATIC ${fft_headers} ${fft_sources})
target_link_libraries(fft /usr/local/lib/libfftw3.a)