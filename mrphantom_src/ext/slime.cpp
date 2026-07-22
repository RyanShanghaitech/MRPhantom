#include <cstdio>
#include <vector>
#include <cstdint>
#include <cmath>
#include <stdexcept>

#include <omp.h>

// Tissue
enum Tissue : uint8_t
{
    Air = 0,
    Fat = 1,
    Myo = 2,
    Blood = 3,
    Liver = 4,
};

bool isInsideEllipsoid
(
    double x, double y, double z,
    double cx, double cy, double cz,
    double rX, double rY, double rZ
)
{
    if (rX <= 0e0 || rY <= 0e0 || rZ <= 0e0) return false;

    const double dx = x - cx;
    const double dy = y - cy;
    const double dz = z - cz;

    const double v =
        (dx*dx) / (rX*rX) +
        (dy*dy) / (rY*rY) +
        (dz*dz) / (rZ*rZ);

    return v <= 1e0;
}

bool genPhant
(
    int64_t nAx, int64_t nZ, int64_t nY, int64_t nX,
    double ampRes, double ampCar,
    std::vector<uint8_t>* voSlime
)
{
    size_t nPixSum = nZ*nY*nX;
    voSlime->resize(nPixSum);
    voSlime->assign(nPixSum, Tissue::Air);

    // shape parameter; r: radius, c: center
    const double dFatOt_rY = nY*400e-3 + nY*ampRes;
    const double dFatOt_rX = nX*400e-3 - 5e-1*nX*ampRes;
    const double dFatOt_rZ = nZ*480e-3;

    const double dFatIn_rY = nY*380e-3 + nY*ampRes;
    const double dFatIn_rX = nX*380e-3 - 5e-1*nX*ampRes;
    const double dFatIn_rZ = nZ*450e-3;

    const double dMyoOt_rY = nY*100e-3 + nY*ampCar;
    const double dMyoOt_rX = nX*120e-3 + nX*ampCar;
    const double dMyoOt_rZ = nZ*100e-3 + nZ*ampCar;;

    const double dMyoIn_rY = nY*60e-3 + nY*(2*ampCar);
    const double dMyoIn_rX = nX*60e-3 + nX*(2*ampCar);
    const double dMyoIn_rZ = nZ*60e-3 + nZ*(2*ampCar);;

    // Centers are at (0,0,0) in your centered coordinate system,
    const double dFat_cx = 0e0, dFat_cy = 0e0, dFat_cz = 0e0;
    const double dMyoOt_cx = 0e0, dMyoOt_cy = 0e0, dMyoOt_cz = 0e0;
    const double dMyoIn_cx = -nX*20e-3, dMyoIn_cy = 0e0, dMyoIn_cz = 0e0;

    /* generate phantom given by `ampRes` and `ampCar` here */
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < nPixSum; ++i)
    {
        // derive coordinates
        const int64_t x = i%nX - nX/2;
        const int64_t y = (i/nX)%nY - nY/2;
        const int64_t z = (nAx==3) ? (i/(nY*nX) - nZ/2) : 0;

        // decide what tissue current pixel is
        if (!isInsideEllipsoid(x,y,z, dFat_cx,dFat_cy,dFat_cz, dFatOt_rX,dFatOt_rY,dFatOt_rZ))
        {
            (*voSlime)[(size_t)i] = (uint8_t)Tissue::Air;
            continue;
        }

        if (!isInsideEllipsoid(x,y,z, dFat_cx,dFat_cy,dFat_cz, dFatIn_rX,dFatIn_rY,dFatIn_rZ))
        {
            (*voSlime)[(size_t)i] = (uint8_t)Tissue::Fat;
            continue;
        }
        
        // vessel balls
        #define V_HIT(cx,cy,cz,div) isInsideEllipsoid(x,y,z, ((cx)*nX), ((cy)*nY), ((cz)*nZ), (dFatIn_rX/(div)), (dFatIn_rY/(div)), (dFatIn_rZ/(div)))

        // 48 “random” vessels (cx,cy,cz in nX,nY,nZ fractions; div in ~[18..44])
        #define VLIST \
            /* 5 on XY plane (z = 0) — slightly farther from heart */ \
            V_HIT( 0.18,  0.09,  0.00, 12) || \
            V_HIT(-0.19,  0.07,  0.00, 16) || \
            V_HIT( 0.11, -0.17,  0.00, 19) || \
            V_HIT(-0.13, -0.16,  0.00, 21) || \
            V_HIT( 0.04,  0.20,  0.00, 22) || \
            /* 5 on YZ plane (x = 0) — adjusted to avoid closeness with XZ plane (y=0) */ \
            V_HIT( 0.00,  0.14,  0.13, 12) || \
            V_HIT( 0.00, -0.15,  0.12, 16) || \
            V_HIT( 0.00,  0.09, -0.17, 19) || \
            V_HIT( 0.00, -0.12, -0.16, 21) || \
            V_HIT( 0.00,  0.23,  0.06, 22) || \
            /* 5 on XZ plane (y = 0) — adjusted corresponding “near” one */ \
            V_HIT( 0.16,  0.00,  0.12, 12) || \
            V_HIT(-0.17,  0.00,  0.11, 16) || \
            V_HIT( 0.13,  0.00, -0.16, 19) || \
            V_HIT(-0.12,  0.00, -0.17, 21) || \
            V_HIT( 0.23,  0.00, -0.06, 22) || \
            /* 8 octants (near heart, slightly pushed outward) */ \
            V_HIT( 0.14,  0.13,  0.12, 16) || /* + + + */ \
            V_HIT(-0.15,  0.13,  0.12, 16) || /* - + + */ \
            V_HIT( 0.14, -0.14,  0.12, 16) || /* + - + */ \
            V_HIT(-0.15, -0.14,  0.12, 16) || /* - - + */ \
            V_HIT( 0.14,  0.13, -0.13, 16) || /* + + - */ \
            V_HIT(-0.15,  0.13, -0.13, 16) || /* - + - */ \
            V_HIT( 0.14, -0.14, -0.13, 16) || /* + - - */ \
            V_HIT(-0.15, -0.14, -0.13, 16)    /* - - - */

        if (VLIST)
        {
            (*voSlime)[(size_t)i] = (uint8_t)Tissue::Fat;
            continue;
        }
        // vessel balls (end)

        if (!isInsideEllipsoid(x,y,z, dMyoOt_cx,dMyoOt_cy,dMyoOt_cz, dMyoOt_rX,dMyoOt_rY,dMyoOt_rZ))
        {
            (*voSlime)[(size_t)i] = (uint8_t)Tissue::Liver;
            continue;
        }

        if (!isInsideEllipsoid(x,y,z, dMyoIn_cx,dMyoIn_cy,dMyoIn_cz, dMyoIn_rX,dMyoIn_rY,dMyoIn_rZ))
        {
            (*voSlime)[(size_t)i] = (uint8_t)Tissue::Myo;
            continue;
        }

        (*voSlime)[(size_t)i] = (uint8_t)Tissue::Blood;
    }

    return true;
}