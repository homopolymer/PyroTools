#ifndef GENERICPROBABILISTICALIGNMENT_H
#define GENERICPROBABILISTICALIGNMENT_H

#include <string>
#include <map>
#include <tuple>
#include <set>
#include <vector>
#include <iostream>
using namespace std;

#include "GenericTool.h"
#include "GenericReadBins.h"
#include "GenericSequenceGlobal.h"
using namespace GenericSequenceTools;

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "utils/bamtools_fasta.h"
using namespace BamTools;

#define USE_MARKOV_FEATURE false
#define USE_NLOPT_LBFGS    true

namespace GenericSequenceTools
{

static string ProbalnStateName[] = {"BEGIN", "MATCH", "MISMATCH", "INSERT", "DELETE", "END"};

typedef int ProbalnState;
#define PROBALN_STATE_SPACE_SIZE 7
#define PROBALN_BEGIN            0
#define PROBALN_MATCH            1
#define PROBALN_MISMATCH         2
#define PROBALN_INSERT           3
#define PROBALN_DELETE           4
#define PROBALN_END              5
#define PROBALN_UNDEFINE         -1

typedef tuple<int,int>     ProbalnTransitionFeature;        // tuple<PREVIOUS_STATE, CURRENT_STATE>
typedef tuple<int,int,int> ProbalnReadHomopolymerFeature;   // tuple<CURRENT_STATE, HOMOPOLYMER_LENGTH>
typedef tuple<int,int,int> ProbalnGenomeHomopolymerFeature; // tuple<CURRENT_STATE, HOMOPOLYMER_LENGTH>
typedef tuple<int,int>     ProbalnNeighborContextFeature;   // tuple<CURRENT_STATE, NEIGHBOR_CONTEXT>
typedef tuple<int,int,int> ProbalnEmissionFeature;          // tuple<CURRENT_STATE, GENOME_ALLELE, READ_ALLELE>

typedef map<ProbalnTransitionFeature, long double>        ProbalnTransitionFeatureValue;
typedef map<ProbalnReadHomopolymerFeature, long double>   ProbalnReadHomopolymerFeatureValue;
typedef map<ProbalnGenomeHomopolymerFeature, long double> ProbalnGenomeHomopolymerFeatureValue;
typedef map<ProbalnNeighborContextFeature, long double>   ProbalnNeighborContextFeatureValue;
typedef map<ProbalnEmissionFeature, long double>          ProbalnEmissionFeatureValue;

struct TrainSetting
{
    int    m_maxIter;
    double m_stopVal;
    double m_band;
    int    m_verbosity;
};

class GenericProbabilisticAlignment
{
    public:
        GenericProbabilisticAlignment();

    public:
        void initTransitionFeature();
        void initReadHomopolymerFeature();
        void initGenomeHomopolymerFeature();
        void initNeighborContextFeature();
        void initEmissionFeature();

        void setTransitionFeature(ProbalnState previousState, ProbalnState currentState);
        void setReadHomopolymerFeature(ProbalnState currentState, int homopolymerLength, int delta);
        void setGenomeHomopolymerFeature(ProbalnState currentState, int homopolymerLength, int delta);
        void setNeighborContextFeature(ProbalnState currentState, int homopolymerLength);
        void setEmissionFeature(ProbalnState currentState, SequenceBase genomeAllele, SequenceBase readAllele);

        void printTransitionFeatureCount(VectorDouble& count);
        void printReadHomopolymerFeatureCount(VectorDouble& count);
        void printGenomeHomopolymerFeatureCount(VectorDouble& count);
        void printNeighborContextFeatureCount(VectorDouble& count);
        void printEmissionFeatureCount(VectorDouble& count);

    public:
        void train(vector<string>& alignReadSeqs, vector<string>& alignGenomeSeqs, TrainSetting& settings);
        void initModelSetting(vector<string>& alignReadSeqs, vector<string>& alignGenomeSeqs);
        void initTransitionFeatureWeight(VectorDouble& counts);
        void initReadHomopolymerFeatureWeight(VectorDouble& counts, VectorDouble& hlen);
        void initGenomeHomopolymerFeatureWeight(VectorDouble& counts, VectorDouble& hlen);
        void initNeighborContextFeatureWeight(VectorDouble& counts, VectorDouble& hlen);
        void initEmissionFeatureWeight(VectorDouble& counts);


    public:
        void print(ostream& output);
        friend ostream& operator << (ostream& output, GenericProbabilisticAlignment& mdl)
        {
            output << "Parameter setting of probabilistic alignment" << endl;
            mdl.print(output);
            return output;
        }

        void read(istream& input);
        friend istream& operator >> (istream& input, GenericProbabilisticAlignment& mdl)
        {
            mdl.read(input);
            return input;
        }

    public:


        long double ViterbiComputation(string& readSeq, string& genomeSeq, double band = 1.0,

                                       string& alignReadSeq   = *static_cast<string*>(0),
                                       string& alignGenomeSeq = *static_cast<string*>(0),

                                       long double& ViterbiScore = *static_cast<long double*>(0),
                                       long double& ForwardScore = *static_cast<long double*>(0),

                                       vector<long double>& TransitionFeatureCounts        = *static_cast<vector<long double>*>(0),
                                       vector<long double>& ReadHomopolymerFeatureCounts   = *static_cast<vector<long double>*>(0),
                                       vector<long double>& GenomeHomopolymerFeatureCounts = *static_cast<vector<long double>*>(0),
                                       vector<long double>& NeighborContextFeatureCounts   = *static_cast<vector<long double>*>(0),
                                       vector<long double>& EmissionFeatureCounts          = *static_cast<vector<long double>*>(0),

                                       vector<long double>& TransitionFeatureExpectations        = *static_cast<vector<long double>*>(0),
                                       vector<long double>& ReadHomopolymerFeatureExpectations   = *static_cast<vector<long double>*>(0),
                                       vector<long double>& GenomeHomopolymerFeatureExpectations = *static_cast<vector<long double>*>(0),
                                       vector<long double>& NeighborContextFeatureExpectations   = *static_cast<vector<long double>*>(0),
                                       vector<long double>& EmissionFeatureExpectations          = *static_cast<vector<long double>*>(0)
                );



        long double featureWeightSum(ProbalnTransitionFeature&        TransitionFeature,
                                     ProbalnReadHomopolymerFeature&   ReadHomopolymerFeature,
                                     ProbalnGenomeHomopolymerFeature& GenomeHomopolymerFeature,
                                     ProbalnNeighborContextFeature&   NeighborContextFeature,
                                     ProbalnEmissionFeature&          EmissionFeature);


    public:
        void getAlignmentFeatureCount(string& alignReadSeq, string& alignGenomeSeq,
                                      vector<long double>& transitionFeatureCount,
                                      vector<long double>& readHomopolymerFeatureCount,
                                      vector<long double>& genomeHomopolymerFeatureCount,
                                      vector<long double>& neighborContextFeatureCount,
                                      vector<long double>& emissionFeatureCount);


    public:
        // number of features
        int                                 NumberTransitionFeatures;
        int                                 NumberReadHomopolymerFeatures;
        int                                 NumberGenomeHomopolymerFeatures;
        int                                 NumberNeighborContextFeatures;
        int                                 NumberEmissionFeatures;

        // set of features
        vector<ProbalnTransitionFeature>                   TransitionFeatureSet;
        vector<ProbalnReadHomopolymerFeature>              ReadHomopolymerFeatureSet;
        vector<ProbalnGenomeHomopolymerFeature>            GenomeHomopolymerFeatureSet;
        vector<ProbalnNeighborContextFeature>              NeighborContextFeatureSet;
        vector<ProbalnEmissionFeature>                     EmissionFeatureSet;

        // feature index
        map<ProbalnTransitionFeature, int>                 TransitionFeatureIndex;
        map<ProbalnReadHomopolymerFeature, int>            ReadHomopolymerFeatureIndex;
        map<ProbalnGenomeHomopolymerFeature, int>          GenomeHomopolymerFeatureIndex;
        map<ProbalnNeighborContextFeature, int>            NeighborContextFeatureIndex;
        map<ProbalnEmissionFeature, int>                   EmissionFeatureIndex;

        // feature weight
        map<ProbalnTransitionFeature, long double>         TransitionFeatureWeight;
        map<ProbalnReadHomopolymerFeature, long double>    ReadHomopolymerFeatureWeight;
        map<ProbalnGenomeHomopolymerFeature, long double>  GenomeHomopolymerFeatureWeight;
        map<ProbalnNeighborContextFeature, long double>    NeighborContextFeatureWeight;
        map<ProbalnEmissionFeature, long double>           EmissionFeatureWeight;

        // feature label
        map<ProbalnTransitionFeature, string>              TransitionFeatureLabel;
        map<ProbalnReadHomopolymerFeature, string>         ReadHomopolymerFeatureLabel;
        map<ProbalnGenomeHomopolymerFeature, string>       GenomeHomopolymerFeatureLabel;
        map<ProbalnNeighborContextFeature, string>         NeighborContextFeatureLabel;
        map<ProbalnEmissionFeature, string>                EmissionFeatureLabel;
};

class ProbabilisticAlignmentTool : public GenericAbstractTool
{
    public:
        ProbabilisticAlignmentTool();

    public:
        int Run(int argc, char *argv[]);
        int Help();

    public:
        int Execute();
        int Train();
        int Realignment();
        int RealignmentAmplicon();
        int Alignment();

    public:
        GenericProbabilisticAlignment ProbAligner;

        bool HasInput;
        bool HasFasta;
        bool HasRegion;
        bool HasConfig;
        bool HasOutput;
        bool HasNumber;
        bool HasRealignMode;
        bool HasAlignMode;
        bool HasTrainMode;
        bool HasCompress;

        string InputFile;
        string FastaFile;
        string RegionStr;
        string ConfigFile;
        string OutputFile;
        string NumberStr;

    public:
        int m_numTrainData;

        BamReader m_bamReader;
        BamRegion m_bamRegion;
        RefVector m_bamGenomes;
        Fasta     m_fasta;
        BamWriter m_bamWriter;

    public:
        bool HasMinLength;
        bool HasMinMapQual;
        bool HasMaxMismatches;

        string MinLengthStr;
        string MinMapQualStr;
        string MaxMismatchesStr;

        int m_minLength;
        int m_minMapQual;
        double m_maxMismatches;

    public:
        bool HasBand;
        string BandStr;
        double m_band;

    public:
        bool   HasMaxIter;
        bool   HasStopVal;

        string MaxIterStr;
        string StopValStr;

        int    m_maxIterNum;
        double m_stopVal;

    public:
        bool   HasVerbose;
        string VerboseStr;
        int    m_verboseLevel;

    public:
        bool   HasUpdateGraph;
        bool   HasTopKPath;
        string TopKPathStr;
        int    m_topK;
        bool   HasHaplotypeLikeRank;

    public:
        bool   HasSeqA;
        string m_seqA;
        bool   HasSeqB;
        string m_seqB;

    public:
        bool   HasAmplicon;

    public:
        bool   HasGraphEdgePruneLevel;
        string GraphEdgePruneLevel;
        int    m_graphEdgePruneLevel;

    public:
        bool   HasDownSample;
        string DownSample;
        int    m_downSample;

    public:
        bool   HasMinSnpRead;
        string MinSnpRead;
        int    m_minSnpRead;

        bool   HasMinInDelRead;
        string MinInDelRead;
        int    m_minInDelRead;

    public:
        bool   HasProcLargeHomopolymer;

    public:
        bool   HasRepeatLength;
        string RepeatLength;
        int    m_repeatLength;

    public:
        bool   HasKeepDuplicate;

    public:
        bool   HasMNP;

    public:
        bool   HasFlankSize;
        string FlankSize;
        int    m_flankSize;
};

}

#endif // GENERICPROBABILISTICALIGNMENT_H
