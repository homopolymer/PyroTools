#ifndef GENERICVARIANTCALL_H
#define GENERICVARIANTCALL_H

#include "utils/bamtools_fasta.h"
#include "api/BamReader.h"
#include "api/BamAux.h"
using namespace BamTools;

#include "GenericVariant.h"
#include "GenericProbabilisticAlignment.h"
using namespace GenericSequenceTools;

#include <set>
#include <string>
#include <vector>
using namespace std;


namespace GenericSequenceTools
{

#define SNPCALL_PRIOR_FULL "full"
#define SNPCALL_PRIOR_FLAT "flat"

struct VariantCallSetting
{
    string m_priorType;
    double m_prior;

    int    m_ploidy;
    double m_variantQualityFilter;

    int    m_topK;

    int    m_flankingSize;

    int    m_graphPruneLevel;

    int    m_band;

    set<string> m_sampleList;
};

class GenericVariantCall
{
    public:
        GenericVariantCall();
        virtual ~GenericVariantCall() {}

    public:
        virtual int call(Fasta& fastaObj, BamReader& bamObj, BamRegion& roi, GenericProbabilisticAlignment& probAligner, VariantCallSetting& snpCallSettings, vector<GenericVariant>& variantSet) = 0;
};

}   // namespace

#endif // GENERICVARIANTCALL_H
