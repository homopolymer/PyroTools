#include "GenericReadBins.h"
using namespace GenericSequenceTools;

#include <sstream>


void GenericReadBins_Clear_Stringstream(std::stringstream& ss)
{
    ss.str(std::string());
    ss.clear();
}

GenericReadBins::GenericReadBins(int binMin, int binMax, int binStep)
    : m_binMin(binMin)
    , m_binMax(binMax)
    , m_binStep(binStep)
{
    // knots of bins
    std::vector<int> binKnots;
    for (int i=m_binMin; i<=m_binMax; i+=m_binStep)
        binKnots.push_back(i);
    binKnots.push_back(std::numeric_limits<int>::max());

    // number of bins
    m_binNum = binKnots.size()-1;

    // the lower and upper boundaries of bins
    for (int i=0; i<m_binNum; i++)
        m_binIntervals.push_back(std::tuple<int,int>(binKnots[i], binKnots[i+1]));

    // the labels of bins
    for (int i=0; i<m_binNum; i++)
    {
        std::stringstream strBuf;
        strBuf << "[";
        strBuf << binKnots[i];
        strBuf << ",";
        if (binKnots[i+1]==std::numeric_limits<int>::max())
            strBuf << "inf";
        else
            strBuf << binKnots[i+1];
        strBuf << ")";
        m_binLabels.push_back(strBuf.str());

        GenericReadBins_Clear_Stringstream(strBuf);
    }

}


GenericReadBins& GenericReadBins::operator =(const GenericReadBins& toCopy)
{
    this->m_binMin  = toCopy.m_binMin;
    this->m_binMax  = toCopy.m_binMax;
    this->m_binStep = toCopy.m_binStep;

    this->m_binNum  = toCopy.m_binNum;
    this->m_binLabels.assign(toCopy.m_binLabels.begin(), toCopy.m_binLabels.end());
    this->m_binIntervals.assign(toCopy.m_binIntervals.begin(), toCopy.m_binIntervals.end());

    return *this;
}


int GenericReadBins::binIndex(int x)
{
    int idx = 0;

    std::vector<std::tuple<int,int>>::iterator iter;
    for (iter=m_binIntervals.begin(); iter!=m_binIntervals.end(); ++iter, ++idx)
    {
        if (x>=std::get<0>(*iter) && x<std::get<1>(*iter))
            break;
    }

    return idx;
}

