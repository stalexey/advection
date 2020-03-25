#pragma once

#include <interpolation.h>
#include <iostream>

template <class T>
void advect(GridData<T>& gridData, const T dx, const InterpolationType type)
{
    const GridData<T> gridDataPrev(gridData);
    const Grid<T>& grid = gridData.grid();

#pragma omp parallel for
    for (int i = 0; i < grid.samples(); ++i) {
        const T x = grid.position(i) - dx;
        gridData[i] = interpolate(gridDataPrev, x, type);
    }
}
