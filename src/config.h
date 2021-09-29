#pragma once

#include <math.h>
#include <Eigen.h>

// define system's float precision
typedef float float_prec;

namespace SYSTEM{
    constexpr uint16_t sampleRate = 50; // [Hz]
    constexpr float_prec dt = static_cast<float_prec>(1.0f / sampleRate);
    constexpr unsigned long dtMicros = static_cast<unsigned long>(dt * 1E6);
}

namespace CONSTANTS{
    constexpr float_prec pi = 3.141592;
    constexpr float_prec r2d = 180.0 / pi;
    constexpr float_prec d2r = pi / 180.0;
}

