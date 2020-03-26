#include <advection.h>

#include <fstream>
#include <iostream>
#include <string>

using T = double;

void
saveData(const GridData<T>& gridData, std::string fileName)
{
    std::ofstream file;
    file.open(fileName.append(".dat"));
    for (int i = 0; i < gridData.size(); i++) {
        file << gridData.position(i) << " " << gridData[i] << std::endl;
    }
    file.close();
}

void
initializeStep(GridData<T>& gridData)
{
    for (int i = 0; i < gridData.size(); i++) {
        const T x = gridData.position(i);
        if (x > 2.5f && x < 7.5f)
            gridData[i] = 1.0f;
        else
            gridData[i] = 0.0f;
    }
}

void
run(const InterpolationType type)
{
    std::string fileName = std::string("advection") + interpolationTypeToString(type);
    std::cout << "Running '" << fileName << "'" << std::endl;

    const T domainSize = 10.0f;
    const int samples = 500;
    const Grid<T> grid(domainSize, samples);
    GridData<T> gridData(grid);
    initializeStep(gridData);

    const T velocity = 5.0f;
    const T cfl = 0.25f;
    const T maxDt = grid.dx() * cfl / velocity;

    const T finalTime = 4.0f;
    const int substeps = std::ceil(finalTime / maxDt);
    const T dt = finalTime / substeps;
    const T dx = dt * velocity;

    for (int i = 0; i < substeps; ++i) { advect(gridData, dx, type); }

    saveData(gridData, fileName);
}

int
main()
{
    for (int i = 0; i < static_cast<int>(InterpolationType::Last); ++i)
        run(static_cast<InterpolationType>(i));

    return 0;
}
