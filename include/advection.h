#pragma once

#include <interpolation.h>

#include <iostream>

enum class AdvectionType
{
    LaxWendroffCDS,
    Last // dummy to indicate the end of list
};

std::string
advectionTypeToString(const AdvectionType& advectionType)
{
    switch (advectionType) {
    case AdvectionType::LaxWendroffCDS:
        return "LaxWendroffCDS";
    default:
        ASSERT(false);
    }
}

enum class SteppingType
{
    SemiLagrangian,
    MacCormack,
    BFECC,
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
    case SteppingType::BFECC:
        return "BFECC";
    default:
        ASSERT(false);
    }
}

template <class T>
void
advectSemiLagrangian(
    const GridData<T>& gridDataIn,
    GridData<T>& gridDataOut,
    const InterpolationType interpolationType,
    const T dx)
{
    ASSERT(gridDataIn.grid() == gridDataOut.grid());
    const Grid<T> grid = gridDataOut.grid();
    const int N = grid.samples();
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        const T x = grid.position(i) - dx;
        gridDataOut[i] = interpolate(gridDataIn, x, interpolationType);
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
    switch (steppingType) {
    case SteppingType::SemiLagrangian: {
        // advect forward
        const GridData<T> gridDataPrev(gridData);
        advectSemiLagrangian(gridDataPrev, gridData, interpolationType, dt * velocity);
    } break;
    case SteppingType::MacCormack: {
        // advect forward
        const GridData<T> gridDataPrev(gridData);
        advectSemiLagrangian(gridDataPrev, gridData, interpolationType, dt * velocity);

        // advect backward
        GridData<T> gridDataNext(gridData);
        advectSemiLagrangian(gridData, gridDataNext, interpolationType, -dt * velocity);

        // correct
        const int N = gridData.grid().samples();
#pragma omp parallel for
        for (int i = 0; i < N; ++i) { gridData[i] += (gridDataPrev[i] - gridDataNext[i]) * 0.5f; }
    } break;
    case SteppingType::BFECC: {
        // advect forward
        GridData<T> gridDataPrev(gridData);
        advectSemiLagrangian(gridDataPrev, gridData, interpolationType, dt * velocity);

        // advect backward
        GridData<T> gridDataNext(gridData);
        advectSemiLagrangian(gridData, gridDataNext, interpolationType, -dt * velocity);

        // correct
        const int N = gridData.grid().samples();
#pragma omp parallel for
        for (int i = 0; i < N; ++i) { gridDataPrev[i] = (3.0f * gridDataPrev[i] - gridDataNext[i]) * 0.5f; }

        // readvect
        advectSemiLagrangian(gridDataPrev, gridData, interpolationType, dt * velocity);
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

    const int N = gridData.grid().samples();
    T alpha = dt * velocity / gridData.grid().dx();
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {

        T result;
        switch (advectionType) {
        case AdvectionType::LaxWendroffCDS: {
            // values
            const T fM = gridDataPrev.periodic(i - 1);
            const T fC = gridDataPrev.periodic(i);
            const T fP = gridDataPrev.periodic(i + 1);

            result = fC - alpha * 0.5 * (fP - fM) + alpha * alpha * 0.5f * (fP - 2.0f * fC + fM);
        } break;
        default:
            ASSERT(false);
        }
        gridData[i] = result;
    }
}
