#include <advection.h>

#include <fstream>
#include <iostream>
#include <string>

using T = double;

const T DOMAIN_SIZE = 10.0f;
const int SAMPLES = 500;
const T VELOCITY = 5.0f;
const T CFL = 0.25f;
const T FINAL_TIME = 4.0f;

const T DX = DOMAIN_SIZE / SAMPLES;
const T MAX_DT = DX * CFL / std::abs(VELOCITY);
const int SUBSTEPS = std::ceil(FINAL_TIME / MAX_DT);
const T DT = FINAL_TIME / SUBSTEPS;

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
run(const SteppingType steppingType, const InterpolationType interpolationType)
{
    std::string fileName = std::string("advection_") + steppingTypeToString(steppingType) + std::string("_") +
                           interpolationTypeToString(interpolationType);
    std::cout << "Running '" << fileName << "'" << std::endl;

    const Grid<T> grid(DOMAIN_SIZE, SAMPLES);
    GridData<T> gridData(grid);
    initializeStep(gridData);

    for (int i = 0; i < SUBSTEPS; ++i) {
        advect(gridData, DT, VELOCITY, steppingType, interpolationType);
        break;
    }

    saveData(gridData, fileName);
}

void
run(const AdvectionType advectionType)
{
    std::string fileName = std::string("advection_") + advectionTypeToString(advectionType);
    std::cout << "Running '" << fileName << "'" << std::endl;

    const Grid<T> grid(DOMAIN_SIZE, SAMPLES);
    GridData<T> gridData(grid);
    initializeStep(gridData);

    for (int i = 0; i < SUBSTEPS; ++i) {
        advect(gridData, DT, VELOCITY, advectionType);
        break;
    }

    saveData(gridData, fileName);
}

void
runReference()
{
    const Grid<T> grid(DOMAIN_SIZE, SAMPLES);
    GridData<T> gridData(grid);
    initializeStep(gridData);

    saveData(gridData, "advection_Reference");
}

int
main()
{
    for (int s = 0; s < static_cast<int>(SteppingType::Last); ++s)
        for (int i = 0; i < static_cast<int>(InterpolationType::Last); ++i)
            run(static_cast<SteppingType>(s), static_cast<InterpolationType>(i));

    for (int a = 0; a < static_cast<int>(AdvectionType::Last); ++a) run(static_cast<AdvectionType>(a));

    runReference();

    return 0;
}
