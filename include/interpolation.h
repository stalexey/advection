#pragma once

#include <grid.h>
#include <utils.h>

#include <string>

enum class InterpolationType
{
    Linear,
    CatmullRom,
    MonotonicCubicFedkiw,
    MonotonicCubicFritschCarlson,
    Last // dummy to indicate the end of list
};

std::string
interpolationTypeToString(const InterpolationType& type)
{
    switch (type) {
    case InterpolationType::Linear:
        return "Linear";
    case InterpolationType::CatmullRom:
        return "CatmullRom";
    case InterpolationType::MonotonicCubicFedkiw:
        return "MonotonicCubicFedkiw";
    case InterpolationType::MonotonicCubicFritschCarlson:
        return "MonotonicCubicFritschCarlson";
    default:
        ASSERT(false);
    }
}

template <class T>
T
interpolate(const GridData<T>& gridData, const T x, const InterpolationType type)
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
        // NOTE: same as what is below without clamping or monotonicity fixes

        const T f0 = gridData.periodic(baseIndex - 1);
        const T f1 = gridData.periodic(baseIndex);
        const T f2 = gridData.periodic(baseIndex + 1);
        const T f3 = gridData.periodic(baseIndex + 2);

        result =
            f1 + 0.5 * alpha *
                     (f2 - f0 +
                      alpha * (2.0 * f0 - 5.0 * f1 + 4.0 * f2 - f3 + alpha * (3.0 * (f1 - f2) + f3 - f0)));
    } break;
    case InterpolationType::MonotonicCubicFedkiw: {
        // From "Visual Simulation of Smoke", Fedkiw et al, 2001
        // C^0 continuity, monotonic, may not look smooth for monotonic regions

        const T f0 = gridData.periodic(baseIndex - 1);
        const T f1 = gridData.periodic(baseIndex);
        const T f2 = gridData.periodic(baseIndex + 1);
        const T f3 = gridData.periodic(baseIndex + 2);

        const T c1 = alpha;
        const T c2 = c1 * alpha;
        const T c3 = c2 * alpha;

        // clamping for non-monotonic region
        T dM, dP;
        T delta = f2 - f1;
        if (delta < 0) {
            dM = std::min((f2 - f0) * 0.5f, T(0.0f));
            dP = std::min((f3 - f1) * 0.5f, T(0.0f));
        } else if (delta > 0) {
            dM = std::max((f2 - f0) * 0.5f, T(0.0f));
            dP = std::max((f3 - f1) * 0.5f, T(0.0f));
        } else {
            dM = dP = 0.0f;
        }

        const T a0 = f1;
        const T a1 = dM;
        const T a2 = 3.0f * delta - 2.0f * dM - dP;
        const T a3 = dM + dP - 2.0f * delta;

        result = a3 * c3 + a2 * c2 + a1 * c1 + a0;

        // just clamp to be between end-points to avoid non-monotonicity
        // not the most ideal solution, but it works
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
        // C^0 continuity, monotonic, looks smooth, does not require clamping to
        // end-points

        const T f0 = gridData.periodic(baseIndex - 1);
        const T f1 = gridData.periodic(baseIndex);
        const T f2 = gridData.periodic(baseIndex + 1);
        const T f3 = gridData.periodic(baseIndex + 2);

        const T c1 = alpha;
        const T c2 = c1 * alpha;
        const T c3 = c2 * alpha;

        // clamping for non-monotonic regions
        T dM, dP;
        T delta = f2 - f1;
        if (delta < 0) {
            dM = std::min((f2 - f0) * 0.5f, T(0.0f));
            dP = std::min((f3 - f1) * 0.5f, T(0.0f));
        } else if (delta > 0) {
            dM = std::max((f2 - f0) * 0.5f, T(0.0f));
            dP = std::max((f3 - f1) * 0.5f, T(0.0f));
        } else {
            dM = dP = 0.0f;
        }

        // clamping for monotonic regions
        const T d2 = dM * dM + dP * dP;
        const T delta2times9 = delta * delta * 9;
        if (d2 > delta2times9) {
            const T scale = std::sqrt(delta2times9 / d2);
            dM *= scale;
            dP *= scale;
        }

        const T a0 = f1;
        const T a1 = dM;
        const T a2 = 3.0f * delta - 2.0f * dM - dP;
        const T a3 = dM + dP - 2.0f * delta;

        result = a3 * c3 + a2 * c2 + a1 * c1 + a0;
    } break;
    default:
        ASSERT(false);
    }

    return result;
}
