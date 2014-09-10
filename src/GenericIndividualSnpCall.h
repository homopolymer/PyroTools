#ifndef GENERICINDIVIDUALSNPCALL_H
#define GENERICINDIVIDUALSNPCALL_H

#include "GenericVariantCall.h"
#include "GenericProbabilisticAlignment.h"
using namespace GenericSequenceTools;

namespace GenericSequenceTools
{

class GenericIndividualSnpCall : public GenericVariantCall
{
    public:
        GenericIndividualSnpCall(int downSample, int minReadLength, int minMapQuality, double maxMismatchFrac, int minSnpRead, int minSnpFrac, int verbosity);

    public:
        int call(Fasta &fastaObj, BamReader &bamObj, BamRegion &roi, GenericProbabilisticAlignment &probAligner, VariantCallSetting& snpCallSettings, vector<GenericVariant> &variantSet);

    public:
        // call SNP by simply scanning the pileup
        void simpleSnpCall(string &fastaObj, BamReader &bamObj, int chrID, int leftPosition, int rightPosition, vector<Allele> &variantCandidates, map<int,list<tuple<char,int,int,double>>> &bamData);

        // call SNP by simple Bayesian method at one site
        void simpleBayesianSnpCall(Fasta& fastaObj, BamReader& bamObj, int chrID, int leftPosition, int rightPosition, list<Allele>& alleles, list<tuple<char,int,int,double>>& data, VariantCallSetting& snpCallSettings, vector<GenericVariant>& variantResult);

        // call SNP by PyroHMMsnp
        void PyroHMMsnp(Fasta &fastaObj, BamReader &bamObj, int chrID, int leftPosition, int rightPosition, GenericProbabilisticAlignment &probAligner, list<Allele> &allelesInBlock, VariantCallSetting& snpCallSettings, vector<GenericVariant> &variantResults);

    public:
        // read filtration
        int    m_minReadLength;
        int    m_minMapQuality;
        int    m_minBaseQuality;
        double m_maxMismatchFrac;

        // sampling
        int    m_downSample;

        // snp filtration
        int    m_minSnpRead;    // minimum number of reads supporting SNP
        double m_minSnpFrac;    // minimum fraction of reads supporting SNP

        int    m_verbosity;

};


class IndividualSnpCallTool : public GenericAbstractTool
{
    public:
        IndividualSnpCallTool();

    public:
        int Run(int argc, char *argv[]);
        int Help();

//    public:
//        GenericIndividualSnpCall IndividualSnpCaller;

    public:
        BamReader m_bamReader;
        BamRegion m_bamRegion;
        Fasta     m_fasta;
        BamWriter m_bamWriter;

    public:
        bool   HasInput;
        bool   HasFasta;
        bool   HasRegion;
        bool   HasConfig;
        bool   HasOutput;

        string InputFile;
        string FastaFile;
        string Region;
        string ConfigFile;
        string OutputFile;

    public:
        bool   HasDownSample;
        string DownSample;
        int    m_downsample;

    public:
        bool   HasMinReadLen;
        string MinReadLen;
        int    m_minReadLen;

        bool   HasMinMapQual;
        string MinMapQual;
        int    m_minMapQual;

        bool   HasMaxMismatchFrac;
        string MaxMismatchFrac;
        double m_maxMismatchFrac;

    public:
        bool   HasMinSnpRead;
        string MinSnpRead;
        int    m_minSnpRead;

        bool   HasMinSnpFrac;
        string MinSnpFrac;
        double m_minSnpFrac;

    public:
        bool   HasVariantQualityFilter;
        string VariantQualityFilter;
        double m_variantQualityFilter;

    public:
        bool   HasPloidy;
        string Ploidy;
        int    m_ploidy;

    public:
        bool   HasPriorType;
        string PriorType;

        bool   HasPrior;
        string Prior;
        double m_prior;

    public:
        bool   HasTopK;
        string TopK;
        int    m_topK;

    public:
        bool   HasFlankingSize;
        string FlankingSize;
        int    m_flankingSize;

    public:
        bool   HasBand;
        string Band;
        double m_band;

    public:
        bool   HasGraphPruneLevel;
        string GraphPruneLevel;
        int    m_graphPruneLevel;

    public:
        bool   HasSample;
        bool   HasSampleList;
        vector<string> SampleList;

    public:
        bool   HasVerbosity;
        string Verbosity;
        int    m_verbosity;
};

}

#endif // GENERICINDIVIDUALSNPCALL_H
