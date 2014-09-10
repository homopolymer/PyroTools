#ifndef GENERICVARIANT_H
#define GENERICVARIANT_H

#include <vector>
#include <list>
#include <string>
#include <tuple>
using namespace std;

namespace GenericSequenceTools
{

#define VARIANT_UNDEFINE "UNDEFINE"
#define VARIANT_SNP      "SNP"
#define VARIANT_INSERT   "INSERT"
#define VARIANT_DELETE   "DELETE"
#define VARIANT_MNP      "MNP"


class GenericAllele
{
    public:
        GenericAllele() {}

    public:
        // allele info
        int    m_chrID;
        int    m_chrPosition;
        string m_allele;
        string m_reference;

        // data info
        int    m_alleleDepth;
        double m_alleleMapAvgQual;
        int    m_alleleStrandPos;
        int    m_alleleStrandNeg;

        // global info
        int    m_globalDepth;
        double m_globalMapAvgQual;
        int    m_globalStrandPos;
        int    m_globalStrandNeg;

        // sequencing data
        list<tuple<char,int,int>> m_bamData;
};

typedef GenericAllele         Allele;
typedef vector<GenericAllele> AlleleSet;

class GenericVariant
{
    public:
        GenericVariant();

    public:
        // variant info
        int         m_chrID;
        int         m_chrPosition;       // 0-based
        AlleleSet   m_alleles;
        string      m_variantType;
        long double m_probScoreRef;
        long double m_probScoreVar;
        long double m_quality;
        string      m_reference;

        vector<int> m_haploidType;
};

}   // namespace

#endif // GENERICVARIANT_H
