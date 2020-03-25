#pragma once

#include <grid.h>
#include <utils.h>

enum class InterpolationType {
    Linear,
    CatmullRom,
    MonotonicCubicFedkiw,
    MonotonicCubicFritschCarlson
};

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
        // C^1 continuity, may exhibit high frequency oscillations

        const T f0 = gridData.periodic(baseIndex - 1);
        const T f1 = gridData.periodic(baseIndex);
        const T f2 = gridData.periodic(baseIndex + 1);
        const T f3 = gridData.periodic(baseIndex + 2);

        result = f1 + 0.5 * alpha *
                          (f2 - f0 +
                           alpha * (2.0 * f0 - 5.0 * f1 + 4.0 * f2 - f3 +
                                    alpha * (3.0 * (f1 - f2) + f3 - f0)));
    } break;
    case InterpolationType::MonotonicCubicFedkiw: {
        // From "Visual Simulation of Smoke", Fedkiw et al, 2001
        // C^0 continuity, monotonic

        const T f0 = gridData.periodic(baseIndex - 1);
        const T f1 = gridData.periodic(baseIndex);
        const T f2 = gridData.periodic(baseIndex + 1);
        const T f3 = gridData.periodic(baseIndex + 2);

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
    case InterpolationType::MonotonicCubicFritschCarlson: {
        // From https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
        // C^0 continuity, monotonic, does not seem to require clamping

        const T f0 = gridData.periodic(baseIndex - 1);
        const T f1 = gridData.periodic(baseIndex);
        const T f2 = gridData.periodic(baseIndex + 1);
        const T f3 = gridData.periodic(baseIndex + 2);

        const T c1 = alpha;
        const T c2 = c1 * alpha;
        const T c3 = c2 * alpha;

        T dkm, dkp;
        T delta = f2 - f1;
        if (delta == 0) {
            dkm = dkp = 0.0f;
        } else {
            dkm = (f2 - f0) * 0.5f;
            dkp = (f3 - f1) * 0.5f;
            if (delta * dkm < 0 || delta * dkp < 0) {
                const T dk2 = dkm * dkm + dkp * dkp;
                const T delta2times9 = delta * delta * 9;
                if (dk2 > delta2times9) {
                    const T scale = std::sqrt(delta2times9 / dk2);
                    dkm *= scale;
                    dkp *= scale;
                }
            }
        }

        const T a0 = f1;
        const T a1 = dkm;
        const T a2 = 3.0f * delta - 2.0f * dkm - dkp;
        const T a3 = dkm + dkp - 2.0f * delta;

        result = a3 * c3 + a2 * c2 + a1 * c1 + a0;
    } break;
    default:
        ASSERT(false);
    }

    return result;
}
