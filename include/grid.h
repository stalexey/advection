#pragma once

#include <cmath>
#include <vector>

template <class T>
class Grid
{
public:
    Grid(const T domainSize, const int samples)
        : m_domainSize(domainSize)
        , m_samples(samples)
        , m_dx(domainSize / samples)
    {
    }
    Grid(const Grid<T>& grid) = default;

    T domainSize() const { return m_domainSize; }
    int samples() const { return m_samples; }
    T dx() const { return m_dx; }

    T position(const int index) const { return (index + 0.5f) * m_dx; }
    void gridSpace(const T x, int& baseIndex, T& alpha) const
    {
        const T gridSpaceX = x / m_dx - 0.5f;
        baseIndex = std::floor(gridSpaceX);
        alpha = gridSpaceX - baseIndex;
    }

    bool operator==(const Grid<T> other) const
    {
        return m_domainSize == other.domainSize() && m_samples == other.samples() && m_dx == other.dx();
    }

protected:
    T m_domainSize;
    int m_samples;
    T m_dx;
};

template <class T>
class GridData
{
public:
    GridData(const Grid<T>& grid)
        : m_grid(grid)
        , m_data(grid.samples(), T())
    {
    }
    GridData(const GridData<T>& gridData) = default;

    const Grid<T>& grid() const { return m_grid; }
    const std::vector<T>& data() const { return m_data; }

    const T& operator[](const int index) const { return m_data[index]; }
    T& operator[](const int index) { return m_data[index]; }

    int size() const { return m_grid.samples(); }
    T position(const int index) const { return m_grid.position(index); }

    const T& periodic(const int index) const
    {
        const int size = m_grid.samples();
        const int periodicIndex = (index % size + size) % size;
        return m_data[periodicIndex];
    }

private:
    Grid<T> m_grid;
    std::vector<T> m_data;
};
