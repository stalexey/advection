#pragma once

#include <interpolation.h>

#include <iostream>

enum class AdvectionType
{
    LaxWendroffCDF,
    Last // dummy to indicate the end of list
};

std::string
advectionTypeToString(const AdvectionType& advectionType)
{
    switch (advectionType) {
    case AdvectionType::LaxWendroffCDF:
        return "LaxWendroffCDF";
    default:
        ASSERT(false);
    }
}

enum class SteppingType
{
    SemiLagrangian,
    MacCormack,
    Last // dummy to indicate the end of list
};

std::string
steppingTypeToString(const SteppingType& steppingType)
{
    switch (steppingType) {
    case SteppingType::SemiLagrangian:
        return "SemiLagrangian";
    case SteppingType::MacCormack:
        return "MacCormack";
    default:
        ASSERT(false);
    }
}

template <class T>
void
advect(
    GridData<T>& gridData,
    const T dt,
    const T velocity,
    const SteppingType steppingType,
    const InterpolationType interpolationType)
{
    const GridData<T> gridDataPrev(gridData);
    const Grid<T>& grid = gridData.grid();

#pragma omp parallel for
    for (int i = 0; i < grid.samples(); ++i) {
        const T x = grid.position(i) - dt * velocity;
        gridData[i] = interpolate(gridDataPrev, x, interpolationType);
    }

    switch (steppingType) {
    case SteppingType::SemiLagrangian:
        break; // do nothing
    case SteppingType::MacCormack: {
        GridData<T> gridDataNext(gridData);
#pragma omp parallel for
        for (int i = 0; i < grid.samples(); ++i) {
            const T x = grid.position(i) + dt * velocity; // advect the other way
            gridDataNext[i] = interpolate(gridData, x, interpolationType);
        }
#pragma omp parallel for
        for (int i = 0; i < grid.samples(); ++i) {
            gridData[i] += (gridDataPrev[i] - gridDataNext[i]) / 2; // apply correction
        }
    } break;
    default:
        ASSERT(false);
    }
}

template <class T>
void
advect(GridData<T>& gridData, const T dt, const T velocity, const AdvectionType advectionType)
{
    const GridData<T> gridDataPrev(gridData);
    const Grid<T>& grid = gridData.grid();

#pragma omp parallel for
    for (int i = 0; i < grid.samples(); ++i) {
        const T x = grid.position(i) - dt * velocity;
        int baseIndex;
        T alpha;
        grid.gridSpace(x, baseIndex, alpha);

        T result;
        switch (advectionType) {
        case AdvectionType::LaxWendroffCDF: {
            // values
            const T fM1 = gridDataPrev.periodic(baseIndex - 1);
            const T fM0 = gridDataPrev.periodic(baseIndex);
            const T fP0 = gridDataPrev.periodic(baseIndex + 1);
            const T fP1 = gridDataPrev.periodic(baseIndex + 2);

            result = fM0 + alpha / 2 * (fP0 - fM1);

            const T ddM = fP0 - 2 * fM0 + fM1;
            result += ddM * alpha * alpha / 2;
        } break;
        default:
            ASSERT(false);
        }
        gridData[i] = result;
    }
}
