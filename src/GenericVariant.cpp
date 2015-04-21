#include "GenericVariant.h"
using namespace GenericSequenceTools;

GenericVariant::GenericVariant()
{
    m_variantType = VARIANT_UNDEFINE;
    m_index = 0;
}

void GenericVariant::AlleleInsert(Allele &allele){
    m_alleleIndex[allele.m_allele] = m_index;
    m_alleles.emplace_back(allele);
    m_index++;
}

int GenericVariant::AlleleIndex(string &allele){
    if (m_alleleIndex.count(allele)==0)
        return -1;
    return m_alleleIndex[allele];
}
