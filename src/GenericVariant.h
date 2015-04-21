#ifndef GENERICVARIANT_H
#define GENERICVARIANT_H

#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <tuple>
#include <map>
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
        GenericAllele(int chrId, int chrPos, string ref, string alt)
            : m_chrID(chrId)
            , m_chrPosition(chrPos)
            , m_allele(alt)
            , m_reference(ref)
        {}

    public:
        // allele info
        int    m_chrID;
        int    m_chrPosition;
        string m_allele;
        string m_reference;

        // data info
        int    m_alleleDepth;
        double m_alleleMapAvgQual;
        int    m_alleleStrandPos;   // positive strand
        int    m_alleleStrandNeg;   // negative strand

        // global info
        int    m_globalDepth;
        double m_globalMapAvgQual;
        int    m_globalStrandPos;
        int    m_globalStrandNeg;

        // sequencing data
        // base-pair resolution information
        // 1:allele
        // 2:base quality
        // 3:map quality
        // 4:strand
        list<tuple<char,int,int,int>> m_bamData;

        static string key(int chrId, int chrPos, char ref, char alt)
        {
            stringstream ss;
            ss << chrId << ":" << chrPos << "_" << ref << ">" << alt;
            return ss.str();
        }
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

        // for allele indexing
        int m_index;
        map<string,int> m_alleleIndex;
        void AlleleInsert(Allele &allele);
        int AlleleIndex(string &allele);

        vector<int> m_haploidType;  // for diploidy

        // sample information
        typedef vector<int> PhasedAlleles;    // order sensitive
        map<string,PhasedAlleles> m_samples;

        // the key value of the variant
        string key(){
            stringstream ss;
            ss << m_chrID << ":" << m_chrPosition;
            return ss.str();
        }
};

}   // namespace

#endif // GENERICVARIANT_H
