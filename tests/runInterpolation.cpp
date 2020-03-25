#include <fstream>
#include <interpolation.h>
#include <iostream>
#include <string>

using T = double;

void saveData(const GridData<T>& gridData, const InterpolationType type,
              const int samples, std::string fileName)
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

void run(const InterpolationType type)
{
    std::string fileName =
        std::string("interpolation") + interpolationTypeToString(type);
    std::cout << "Running '" << fileName << "'" << std::endl;

    const T domainSize = 10.0f;
    const int samples = 10;
    const Grid<T> grid(domainSize, samples);
    GridData<T> gridData(grid);
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

    saveData(gridData, type, 1000, fileName);
}

int main()
{
    for (int i = 0; i < static_cast<int>(InterpolationType::Last); ++i)
        run(static_cast<InterpolationType>(i));

    return 0;
}
