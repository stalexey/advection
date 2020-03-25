#include <fstream>
#include <interpolation.h>
#include <iostream>
#include <string>

using T = double;

void saveData(const GridData<T>& gridData, std::string fileName,
              const InterpolationType type, const int samples)
{
    std::ofstream file;
    file.open(fileName.append(".dat"));
    const T dx = gridData.grid().domainSize() / samples;
    for (int i = 0; i < samples; i++) {
        const T x = (i + 0.5f) * dx;
        const T y = interpolate(gridData, x, type);
        file << x << " " << y << std::endl;
    }
    file.close();
}

void initialize(GridData<T>& gridData)
{
    gridData[0] = 3;
    gridData[1] = 2.9;
    gridData[2] = 2.5;
    gridData[3] = 1;
    gridData[4] = 0.9;
    gridData[5] = 0.8;
    gridData[6] = 0.5;
    gridData[7] = 0.2;
    gridData[8] = 0.1;
    gridData[9] = 0.1;
}

void run(const InterpolationType type, const std::string& fileName)
{
    std::cout << "Running interpolation for '" << fileName << "'" << std::endl;

    const T domainSize = 10.0f;
    const int samples = 10;
    const Grid<T> grid(domainSize, samples);
    GridData<T> gridData(grid);
    initialize(gridData);
    saveData(gridData, fileName, type, 1000);
}

int main()
{
    run(InterpolationType::Linear, "interpolationLinear");
    run(InterpolationType::CatmullRom, "interpolationCatmullRom");
    run(InterpolationType::MonotonicCubicFedkiw,
        "interpolationMonotonicCubicFedkiw");
    run(InterpolationType::MonotonicCubicFritschCarlson,
        "interpolationMonotonicCubicFritschCarlson");

    return 0;
}
