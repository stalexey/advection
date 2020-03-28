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
interpolationTypeToString(const InterpolationType& interpolationType)
{
    switch (interpolationType) {
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
cubicInterpolant(const T x, const T fM, const T fP, const T dM, const T dP)
{
    // powers
    const T x1 = x;
    const T x2 = x1 * x;
    const T x3 = x2 * x;

    const T delta = fP - fM;
    const T c0 = fM;
    const T c1 = dM;
    const T c2 = 3.0f * delta - 2.0f * dM - dP;
    const T c3 = dM + dP - 2.0f * delta;

    return c3 * x3 + c2 * x2 + c1 * x1 + c0;
}

template <class T>
T
interpolate(const GridData<T>& gridData, const T x, const InterpolationType interpolationType)
{
    const Grid<T>& grid = gridData.grid();
    int baseIndex;
    T alpha;
    grid.gridSpace(x, baseIndex, alpha);

    T result;
    switch (interpolationType) {
    case InterpolationType::Linear: {
        const T fM = gridData.periodic(baseIndex);
        const T fP = gridData.periodic(baseIndex + 1);
        result = fP * alpha + fM * (1 - alpha);
    } break;
    case InterpolationType::CatmullRom: {
        // Reference implementation
        // https://www.paulinternet.nl/?page=bicubic
        // C^1 continuity, may exhibit high frequency oscillations
        // NOTE: same as what is below without clamping or monotonicity fixes

        // values
        const T fM1 = gridData.periodic(baseIndex - 1);
        const T fM0 = gridData.periodic(baseIndex);
        const T fP0 = gridData.periodic(baseIndex + 1);
        const T fP1 = gridData.periodic(baseIndex + 2);

        // derivatives
        const T dM = (fP0 - fM1) * 0.5f;
        const T dP = (fP1 - fM0) * 0.5f;

        result = cubicInterpolant(alpha, fM0, fP0, dM, dP);
    } break;
    case InterpolationType::MonotonicCubicFedkiw: {
        // From "Visual Simulation of Smoke", Fedkiw et al, 2001
        // C^0 continuity, monotonic, may not look smooth for monotonic regions

        // values
        const T fM1 = gridData.periodic(baseIndex - 1);
        const T fM0 = gridData.periodic(baseIndex);
        const T fP0 = gridData.periodic(baseIndex + 1);
        const T fP1 = gridData.periodic(baseIndex + 2);

        // derivatives
        T dM = (fP0 - fM1) * 0.5f;
        T dP = (fP1 - fM0) * 0.5f;

        // clamping for non-monotonic region
        const T delta = fP0 - fM0;
        if (delta < 0) {
            dM = std::min(dM, T(0.0f));
            dP = std::min(dP, T(0.0f));
        } else if (delta > 0) {
            dM = std::max(dM, T(0.0f));
            dP = std::max(dP, T(0.0f));
        } else {
            dM = dP = 0.0f;
        }

        result = cubicInterpolant(alpha, fM0, fP0, dM, dP);

        // just clamp to be between end-points to avoid non-monotonicity
        // not the most ideal solution, but it works
        result = std::min(result, std::max(fM0, fP0));
        result = std::max(result, std::min(fM0, fP0));
    } break;
    case InterpolationType::MonotonicCubicFritschCarlson: {
        // From https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
        // C^0 continuity, monotonic, looks smooth, does not require clamping to
        // end-points

        const T fM1 = gridData.periodic(baseIndex - 1);
        const T fM0 = gridData.periodic(baseIndex);
        const T fP0 = gridData.periodic(baseIndex + 1);
        const T fP1 = gridData.periodic(baseIndex + 2);

        // derivatives
        T dM = (fP0 - fM1) * 0.5f;
        T dP = (fP1 - fM0) * 0.5f;

        // clamping for non-monotonic region
        const T delta = fP0 - fM0;
        if (delta < 0) {
            dM = std::min(dM, T(0.0f));
            dP = std::min(dP, T(0.0f));
        } else if (delta > 0) {
            dM = std::max(dM, T(0.0f));
            dP = std::max(dP, T(0.0f));
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

        result = cubicInterpolant(alpha, fM0, fP0, dM, dP);
    } break;
    default:
        ASSERT(false);
    }

    return result;
}
