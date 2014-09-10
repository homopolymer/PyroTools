#ifndef GENERICVCFTOOLS_H
#define GENERICVCFTOOLS_H

#include "api/BamAux.h"
using namespace BamTools;

#include "GenericVariant.h"
#include "GenericVariantCall.h"
using namespace GenericSequenceTools;

namespace GenericSequenceTools
{
class GenericVcfTools
{
    public:
        GenericVcfTools();

    public:
        static void write(RefVector& chromosomes, string& vcfFileName, VariantCallSetting& variantCallSettings, vector<GenericVariant>& variantsToReport);
};
}   // namespace

#endif // GENERICVCFTOOLS_H
