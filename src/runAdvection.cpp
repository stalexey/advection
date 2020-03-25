#include <advection.h>
#include <fstream>
#include <iostream>
#include <string>

using T = double;

void saveData(const GridData<T>& gridData, std::string fileName)
{
    std::ofstream file;
    file.open(fileName.append(".dat"));
    for (int i = 0; i < gridData.size(); i++) {
        file << gridData.position(i) << " " << gridData[i] << std::endl;
    }
    file.close();
}

void initializeStep(GridData<T>& gridData)
{
    for (int i = 0; i < gridData.size(); i++) {
        const T x = gridData.position(i);
        if (x > 2.5f && x < 7.5f)
            gridData[i] = 1.0f;
        else
            gridData[i] = 0.0f;
    }
}

void run(InterpolationType type, std::string fileName)
{
    std::cout << "Running advection for '" << fileName << "'" << std::endl;

    const T domainSize = 10.0f;
    const int samples = 200;
    const Grid<T> grid(domainSize, samples);
    GridData<T> gridData(grid);
    initializeStep(gridData);

    const T velocity = 5.0f;
    const T cfl = 6.3f;
    const T maxDt = grid.dx() * cfl / velocity;

    const T finalTime = 1.0f;
    const int substeps = std::ceil(finalTime / maxDt);
    const T dt = finalTime / substeps;
    const T dx = dt * velocity;

    for (int i = 0; i < substeps; ++i) {
        advect(gridData, dx, type);
    }

    saveData(gridData, fileName);
}

int main()
{
    run(InterpolationType::Linear, "interpolationLinear");
    run(InterpolationType::CatmullRom, "interpolationCatmullRom");
    run(InterpolationType::MonotonicCubicFedkiw, "interpolationMonotonicCubicFedkiw");
    run(InterpolationType::MonotonicCubicFritschCarlson, "interpolationMonotonicCubicFritschCarlson");

    return 0;
}