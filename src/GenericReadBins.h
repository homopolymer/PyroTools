#ifndef GENERICREADBINS_H
#define GENERICREADBINS_H

#include <tuple>
#include <vector>
#include <string>
#include <limits>

namespace GenericSequenceTools{

class GenericReadBins
{
    // ctor & dtor
    public:
        GenericReadBins(int binMin, int binMax, int binStep);
        ~GenericReadBins() {}

    public:
        int binIndex(int x);
        GenericReadBins& operator=(const GenericReadBins& toCopy);

    public:
        int m_binNum;
        std::vector<std::string> m_binLabels;
        std::vector<std::tuple<int,int>> m_binIntervals;

    public:
        int m_binMin;
        int m_binMax;
        int m_binStep;
};

}  // namespace

#endif // READBINS_H
