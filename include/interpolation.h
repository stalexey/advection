#pragma once

#include <grid.h>
#include <utils.h>

enum class InterpolationType { Linear, CatmullRom, MonotoneCubic };

template <class T>
T interpolate(const GridData<T>& gridData, const T x,
              const InterpolationType type)
{
    const Grid<T>& grid = gridData.grid();
    int baseIndex;
    T alpha;
    grid.gridSpace(x, baseIndex, alpha);

    T result;
    switch (type) {
    case InterpolationType::Linear: {
        const T f0 = gridData.periodic(baseIndex);
        const T f1 = gridData.periodic(baseIndex + 1);
        result = f1 * alpha + f0 * (1 - alpha);
    } break;
    case InterpolationType::CatmullRom: {
        // Reference implementation
        // https://www.paulinternet.nl/?page=bicubic

        const T p0 = gridData.periodic(baseIndex - 1);
        const T p1 = gridData.periodic(baseIndex);
        const T p2 = gridData.periodic(baseIndex + 1);
        const T p3 = gridData.periodic(baseIndex + 1);

        result = p1 + 0.5 * alpha *
                          (p2 - p0 +
                           alpha * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3 +
                                    alpha * (3.0 * (p1 - p2) + p3 - p0)));

    } break;
    case InterpolationType::MonotoneCubic: {
        const T f0 = gridData.periodic(baseIndex - 1);
        const T f1 = gridData.periodic(baseIndex);
        const T f2 = gridData.periodic(baseIndex + 1);
        const T f3 = gridData.periodic(baseIndex + 1);

        const T c1 = alpha;
        const T c2 = c1 * alpha;
        const T c3 = c2 * alpha;

        T dkm, dkp;
        T delta = f2 - f1;
        if (delta < 0) {
            dkm = std::min((f2 - f0) * 0.5f, T(0.0f));
            dkp = std::min((f3 - f1) * 0.5f, T(0.0f));
        } else if (delta > 0) {
            dkm = std::max((f2 - f0) * 0.5f, T(0.0f));
            dkp = std::max((f3 - f1) * 0.5f, T(0.0f));
        } else {
            dkm = dkp = 0.0f;
        }

        const T a0 = f1;
        const T a1 = dkm;
        const T a2 = 3.0f * delta - 2.0f * dkm - dkp;
        const T a3 = dkm + dkp - 2.0f * delta;

        result = a3 * c3 + a2 * c2 + a1 * c1 + a0;

        if (delta < 0) {
            if (result > f1)
                result = f1;
            else if (result < f2)
                result = f2;
        } else if (delta > 0) {
            if (result > f2)
                result = f2;
            else if (result < f1)
                result = f1;
        }

    } break;
    default:
        ASSERT(false);
    }

    return result;
}
