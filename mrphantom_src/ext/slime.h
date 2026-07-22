#pragma once

#include <vector>
#include <cstdint>

extern bool genPhant
(
    int64_t nAx, int64_t nZ, int64_t nY, int64_t nX,
    double ampRes, double ampCar,
    std::vector<uint8_t>* voSlime
);