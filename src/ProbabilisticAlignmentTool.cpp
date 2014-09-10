#include "GenericReadBins.h"
#include "GenericProbabilisticAlignment.h"
#include "GenericBamAlignmentTools.h"
#include "GenericFastaTools.h"
#include "GenericGraphTools.h"
using namespace GenericSequenceTools;

#include "api/BamAlignment.h"
#include "utils/bamtools_options.h"
#include "utils/bamtools_utilities.h"
using namespace BamTools;

#include <ctime>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
using namespace std;

/**
 * @brief ProbabilisticAlignmentTool::ProbabilisticAlignmentTool
 * the constructor of the re-alignme tool is to set up the program options, and
 * initial the internal variables.
 */
ProbabilisticAlignmentTool::ProbabilisticAlignmentTool()
    : HasInput(false)
    , HasFasta(false)
    , HasRegion(false)
    , HasConfig(false)
    , HasOutput(false)
    , HasNumber(false)
    , HasRealignMode(false)
    , HasAlignMode(false)
    , HasTrainMode(false)
    , HasCompress(false)
    , m_numTrainData(500)
    , HasMinLength(false)
    , HasMinMapQual(false)
    , HasMaxMismatches(false)
    , m_minLength(100)
    , m_minMapQual(10)
    , m_maxMismatches(0.03)
    , HasBand(false)
    , m_band(0.1)
    , HasMaxIter(false)
    , HasStopVal(false)
    , m_maxIterNum(30)
    , m_stopVal(1e-1)
    , HasVerbose(false)
    , m_verboseLevel(0)
    , HasUpdateGraph(false)
    , HasTopKPath(false)
    , m_topK(2)
    , HasHaplotypeLikeRank(false)
    , HasSeqA(false)
    , HasSeqB(false)
    , HasAmplicon(false)
    , HasGraphEdgePruneLevel(false)
    , m_graphEdgePruneLevel(1)
    , HasDownSample(false)
    , m_downSample(1e6)
    , HasMinSnpRead(false)
    , m_minSnpRead(2)
    , HasMinInDelRead(false)
    , m_minInDelRead(2)
    , HasProcLargeHomopolymer(false)
    , HasRepeatLength(false)
    , m_repeatLength(6)
    , HasKeepDuplicate(false)
    , HasMNP(false)
    , HasFlankSize(false)
    , m_flankSize(5)
{
    // set program details
    Options::SetProgramInfo("GenericSequenceTools paln",
                            "probabilistic re-alignment tool",
                            "[-region <REGION>] [-config <CONFIG>] [-train] [-compress] -bam <INPUT> -ref <FASTA> [-out <OUTPUT>]");

    // set option group
    OptionGroup* SPG_Opts = Options::CreateOptionGroup("Subprogram Options");
    Options::AddOption("-align", "alignment of sequence B against sequence A", HasAlignMode, SPG_Opts);
    Options::AddOption("-realign", "re-alignment over BAM file", HasRealignMode, SPG_Opts);
    Options::AddOption("-train","estimate the alignment parameters", HasTrainMode, SPG_Opts);

    // set option group
    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input and Output Options");
    Options::AddValueOption("-bam", "FILE", "the input BAM file", "", HasInput, InputFile, IO_Opts);
    Options::AddValueOption("-ref", "FILE", "the genome file", "", HasFasta, FastaFile, IO_Opts);
    Options::AddValueOption("-a","STR","sequence A","",HasSeqA,m_seqA,IO_Opts);
    Options::AddValueOption("-b","STR","sequence B","",HasSeqB,m_seqB,IO_Opts);
    Options::AddValueOption("-out", "FILE", "the output BAM file", "", HasOutput, OutputFile, IO_Opts);
    Options::AddOption("-c", "compressed output", HasCompress, IO_Opts);
    Options::AddValueOption("-config", "FILE", "model configure file", "", HasConfig, ConfigFile, IO_Opts);
    Options::AddValueOption("-down-sample","INT","sample data down to the specified depth [null]","",HasDownSample,DownSample,IO_Opts);
    Options::AddOption("-amplicon", "data is AmpliconSeq", HasAmplicon, IO_Opts);

    // set option group
    OptionGroup* Rg_Opts = Options::CreateOptionGroup("Genome Region Options");
    Options::AddValueOption("-region", "STR", "set region of interesting, e.g. chr:p0..chr:p1 for a segment, chr:0 for a chromosome", "", HasRegion, RegionStr, Rg_Opts);

    // set option group
    OptionGroup* Flt_Opts = Options::CreateOptionGroup("Filtration Options");
    Options::AddOption("-keep-duplicate", "keep duplicates for PacBio", HasKeepDuplicate, Flt_Opts);
    Options::AddOption("-include-mnp", "Also realign in MNP sites", HasMNP, Flt_Opts);
    Options::AddValueOption("-min-repeat-len", "INT", "the minimum number of repetitive elements [6]", "", HasRepeatLength, RepeatLength, Flt_Opts);
    Options::AddValueOption("-min-read-len", "INT", "the minimum read length [100]", "", HasMinLength, MinLengthStr, Flt_Opts);
    Options::AddValueOption("-min-map-qual", "INT", "the minimum map quality score [10]", "", HasMinMapQual, MinMapQualStr, Flt_Opts);
    Options::AddValueOption("-max-mismatches", "FLT", "the maximum ratio of mismatches in read [0.03]", "", HasMaxMismatches, MaxMismatchesStr, Flt_Opts);
    Options::AddValueOption("-min-snp-read", "INT", "the minimum number of reads supporting SNPs [2]", "", HasMinSnpRead, MinSnpRead, Flt_Opts);
    Options::AddValueOption("-min-indel-read", "INT", "the minimum number of reads supporting InDels [2]", "", HasMinInDelRead, MinInDelRead, Flt_Opts);

    // set option group
    OptionGroup* Gh_Opts = Options::CreateOptionGroup("Realignment Options");
    Options::AddValueOption("-edge-prune-level", "INT", "prune away graph edges below the threshold [1]", "", HasGraphEdgePruneLevel, GraphEdgePruneLevel, Gh_Opts);
    Options::AddOption("-update-graph-by-read", "turn on graph updating by new computed alignment", HasUpdateGraph, Gh_Opts);
    Options::AddValueOption("-top-k", "INT", "the number of top ranked paths that will be consensus sequences [2]", "", HasTopKPath, TopKPathStr, Gh_Opts);
    Options::AddOption("-haplo-like-rank", "turn on the likelihood analysis of haplotype (consensus) sequences, ranking the haplotypes by the likelihood scores in the place of path weights", HasHaplotypeLikeRank, Gh_Opts);
    Options::AddOption("-large-homopolymer", "for 454 and Ion Torrent", HasProcLargeHomopolymer, Gh_Opts);
    Options::AddValueOption("-flank-size", "INT", "size of flanking fragment at both sides [5]", "", HasFlankSize, FlankSize, Gh_Opts);

    // set option group
    OptionGroup* Mdl_Opts = Options::CreateOptionGroup("Training Options");
    Options::AddValueOption("-band", "FLT", "the band of alignment", "", HasBand, BandStr, Mdl_Opts);
    Options::AddValueOption("-num", "INT", "the size of training data [500]", "", HasNumber, NumberStr, Mdl_Opts);
    Options::AddValueOption("-max-iter", "INT", "maximum iterations of iterative training [30]", "", HasMaxIter, MaxIterStr, Mdl_Opts);
    Options::AddValueOption("-stop-val", "FLT", "minimum change value of objective function evaluation [0.1]", "", HasStopVal, StopValStr, Mdl_Opts);

    // set option group
    OptionGroup* Vbs_Opts = Options::CreateOptionGroup("Verbosity Options");
    Options::AddValueOption("-verbosity", "INT", "verbosity level", "", HasVerbose, VerboseStr, Vbs_Opts);
}

int ProbabilisticAlignmentTool::Help()
{
    Options::DisplayHelp();
    return 0;
}

int ProbabilisticAlignmentTool::Run(int argc, char *argv[])
{
    // just 'GenericSequenceTools paln'
    if (argc==2)
        return Help();

    // parse the command line options
    Options::Parse(argc, argv, 1);

    if (HasRealignMode || HasTrainMode)
    {
        // -------------------------------
        // BAM file input

        // BAM file
        if (!HasInput){
            cerr << "GenericSequenceTools paln ERROR: BAM file not provided... Aborting" << endl;
            return false;
        }

        // open input files
        if (!m_bamReader.Open(InputFile)){
            cerr << "GenericSequenceTools paln ERROR: could not open input BAM file... Aborting" << endl;
            return false;
        }

        // if input is not stdin and a region is provided, look for index files
        if (HasInput && HasRegion){
            if (!m_bamReader.LocateIndex()){
                cerr << "GenericSequenceTools paln ERROR: could not locate index file... Aborting" << endl;
                return false;
            }
        }

        // retrieve reference data
        m_bamGenomes = m_bamReader.GetReferenceData();

        // set region if specified
        if (HasRegion){
            if (BamTools::Utilities::ParseRegionString(RegionStr, m_bamReader, m_bamRegion)){
                if (m_bamReader.HasIndex()){
                    if (!m_bamReader.SetRegion(m_bamRegion)){
                        cerr << "GenericSequenceTools paln ERROR: set region failed... Aborting" << endl;
                        m_bamReader.Close();
                        return false;
                    }
                }

            }else{
                cerr << "GenericSequenceTools paln ERROR: could not parse REGION: " << RegionStr << endl;
                m_bamReader.Close();
                return false;
            }
        }

        //---------------------------------------------
        // Fasta file input

        // Fasta file must be provided
        if (FastaFile.empty()){
            cerr << "GenericSequenceTools paln ERROR: Fasta file not provided... Aborting" << endl;
            return false;
        }

        // check Fasta index
        string FastaIndexFilename = "";
        if (BamTools::Utilities::FileExists(FastaFile+".fai"))
            FastaIndexFilename = FastaFile + ".fai";

        // open Fasta file
        if (!m_fasta.Open(FastaFile, FastaIndexFilename)){
            cerr << "GenericSequenceTools paln ERROR: " << FastaFile << " "
                 << "can not be open... Aborting" << endl;
            return false;
        }


        // ---------------------------------------------
        // number of training data
        if (HasNumber)
            m_numTrainData = stoi(NumberStr);

        //---------------------------------------
        // filtration
        if (HasRepeatLength)
            m_repeatLength = stoi(RepeatLength);

        if (HasMinLength)
            m_minLength = stoi(MinLengthStr);

        if (HasMinMapQual)
            m_minMapQual = stoi(MinMapQualStr);

        if (HasMaxMismatches)
            m_maxMismatches = stod(MaxMismatchesStr);

        if (HasMinSnpRead)
            m_minSnpRead = stoi(MinSnpRead);

        if (HasMinInDelRead)
            m_minInDelRead = stoi(MinInDelRead);

        // ---------------------------------------
        // graph pruning
        if (HasGraphEdgePruneLevel)
            m_graphEdgePruneLevel = stoi(GraphEdgePruneLevel);

        // ---------------------------------------
        // top-k
        if (HasTopKPath)
            m_topK = stoi(TopKPathStr);

        // ---------------------------------------
        if (HasFlankSize)
            m_flankSize = stoi(FlankSize);

        // ---------------------------------------
        // downsampling
        if (HasDownSample)
            m_downSample = stoi(DownSample);
    }

    // ---------------------------------------
    // band
    if (HasBand)
        m_band = stod(BandStr);

    // ---------------------------------------
    // verbosity
    if (HasVerbose)
        m_verboseLevel = stoi(VerboseStr);



    // executation
    if (Execute())
        return 1;

    return 0;
}

int ProbabilisticAlignmentTool::Execute()
{
    if (HasTrainMode)
    {
        if (!Train())
            return 0;
    }else if (HasRealignMode)
    {
        if (HasAmplicon)
        {
            if (!RealignmentAmplicon())
                return 0;
        }else
        {
            if (!Realignment())
                return 0;
        }
    }else
    {
        if (!Alignment())
            return 0;
    }
    return 1;
}


bool readIsFiltered(BamAlignment& alignObj, bool keepDuplicate, int& minReadLen=*static_cast<int*>(0), int& minMapQual=*static_cast<int*>(0), double& maxMismatchFrac=*static_cast<double*>(0))
{
    // skip if it is bad alignment
    if (!GenericBamAlignmentTools::goodAlignment(alignObj, keepDuplicate))
    {
        return true;
    }

    // skip if there is no MD tag
    if (!alignObj.HasTag("MD"))
    {
        return true;
    }

    // skip if it is too short
    if (&minReadLen!=0)
    {
        int len = GenericBamAlignmentTools::getBamAlignmentReadLength(alignObj);
        if (len<minReadLen)
        {
            return true;
        }
    }

    // skip if map quality is low scored
    if (&minMapQual!=0)
    {
        if (alignObj.MapQuality<minMapQual)
        {
            return true;
        }
    }

    // skip if too many mismatches
    if (&maxMismatchFrac!=0)
    {
        int len = GenericBamAlignmentTools::getBamAlignmentReadLength(alignObj);
        int nMs = GenericBamAlignmentTools::numBamAlignmentMismatches(alignObj);
        if (nMs>len*maxMismatchFrac)
        {
            return true;
        }
    }

    return false;
}


/**
 * @brief ProbabilisticAlignmentTool::Train
 * @return
 * It is designed to train the conditional random field by using the BAM alignments in a given region.
 */
int ProbabilisticAlignmentTool::Train()
{
    // make sure output configure file is provided
    if (!HasConfig)
    {
        cerr << "GenericSequenceTools paln ERROR: no output configure file provided... Aborting" << endl;
        return 0;
    }

    // -----------------------------------------------------
    // select training data

    int rgTotalAlns;

    BamAlignment al;

    vector<string> alignReadSeqs;
    vector<string> alignGenomeSeqs;


    // random dicer
    unsigned UniformSeed = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 UniformGenerator(UniformSeed);
    uniform_int_distribution<int> UniformDist;

    // total number of alignments
    rgTotalAlns = 0;
    while (m_bamReader.GetNextAlignment(al))
    {
        rgTotalAlns ++;
    }

    // rewind bam reader
    m_bamReader.Rewind();

    // reset bam region
    if(HasRegion)
        m_bamReader.SetRegion(m_bamRegion);

    // iterate
    UniformDist = uniform_int_distribution<int>(0, rgTotalAlns-1);
    auto UniformDicer = bind(UniformDist, UniformGenerator);

    int rgThresh = (m_numTrainData>rgTotalAlns*0.1) ? m_numTrainData : rgTotalAlns*0.1;

    while (m_bamReader.GetNextAlignment(al)) {

        // skip if read fail quality assessment
        if (readIsFiltered(al, m_minLength, m_minMapQual))
            continue;

        int len = GenericBamAlignmentTools::getBamAlignmentReadLength(al);

        if (UniformDicer()<rgThresh)
        {
            string alignRead, alignGenome;
            GenericBamAlignmentTools::getAlignmentSequences(al, alignRead, alignGenome);

            // skip if genome has ambiguous base
            int numAmbiguous = 0;
            for (string::iterator iter=alignGenome.begin(); iter!=alignGenome.end(); ++iter)
            {
                if ((*iter)==Amb)
                    numAmbiguous++;
            }
            if (numAmbiguous>0)
                continue;

            // save the item
            int numBlock;
            int sizeBlock;

            if (len>100)
            {
                numBlock  = ceil(len/100);
                sizeBlock = len/numBlock;

                for (int i=0; i<alignRead.length(); i+=sizeBlock)
                {
                    if (i+sizeBlock>alignRead.length())
                        continue;

                    string alignReadBlock   = alignRead.substr(i, sizeBlock);
                    string alignGenomeBlock = alignGenome.substr(i, sizeBlock);

                    // too many mismatches
                    int t_numMatch    = GenericBamAlignmentTools::numMatch(alignReadBlock, alignGenomeBlock);
                    int t_numMismatch = GenericBamAlignmentTools::numMismatch(alignReadBlock, alignGenomeBlock);
                    if (t_numMismatch>(t_numMatch+t_numMismatch)*m_maxMismatches)
                        continue;

                    // skip if head part is not equal
                    if (alignReadBlock.substr(0, 5)!=alignGenomeBlock.substr(0,5))
                        continue;

                    // skip if there is N
                    if (GenericFastaTools::numAmbiguousBase(alignReadBlock)>0)
                        continue;
                    if (GenericFastaTools::numAmbiguousBase(alignGenomeBlock)>0)
                        continue;

                    // compute the length difference between read and genome
                    int lr=0, lg=0;
                    for (string::iterator iter=alignReadBlock.begin(); iter!=alignReadBlock.end(); ++iter)
                    {
                        if ((*iter)!=Spa)
                            lr++;
                    }
                    for (string::iterator iter=alignGenomeBlock.begin(); iter!=alignGenomeBlock.end(); ++iter)
                    {
                        if ((*iter)!=Spa)
                            lg++;
                    }
                    int delta = abs(lr-lg);

                    // skip if length difference is too much
                    // because that it may be a structural variant
                    if (delta>m_band*0.5*(lr+lg) || delta>5)
                        continue;

                    // skip if gap is large
                    if (GenericBamAlignmentTools::maxGapLength(alignReadBlock, alignGenomeBlock)>5)
                        continue;

                    // ease '-' at the end
                    for (int t=alignReadBlock.length()-1; t>=0; --t)
                    {
                        if (alignReadBlock[t]!=Spa && alignGenomeBlock[t]!=Spa)
                            break;

                        alignReadBlock.erase(t, 1);
                        alignGenomeBlock.erase(t, 1);
                    }

                    // save
                    alignReadSeqs.push_back(alignReadBlock);
                    alignGenomeSeqs.push_back(alignGenomeBlock);
                }

            }else
            {

                // skip if too many mismatch
                int t_numMatch    = GenericBamAlignmentTools::numMatch(alignRead, alignGenome);
                int t_numMismatch = GenericBamAlignmentTools::numMismatch(alignRead, alignGenome);
                if (t_numMismatch>(t_numMatch+t_numMismatch)*m_maxMismatches)
                    continue;

                // skip if head part is not equal
                if (alignRead.substr(0,5)!=alignGenome.substr(0,5))
                    continue;

                // skip if there is N
                if (GenericFastaTools::numAmbiguousBase(alignRead)>0)
                    continue;
                if (GenericFastaTools::numAmbiguousBase(alignGenome)>0)
                    continue;

                // compute the length difference between read and genome
                int lr=0, lg=0;
                for (string::iterator iter=alignRead.begin(); iter!=alignRead.end(); ++iter)
                {
                    if ((*iter)!=Spa)
                        lr++;
                }
                for (string::iterator iter=alignGenome.begin(); iter!=alignGenome.end(); ++iter)
                {
                    if ((*iter)!=Spa)
                        lg++;
                }
                int delta = abs(lr-lg);

                // skip if length difference is too much
                // because that it may be a structural variant
                if (delta>m_band*0.5*(lr+lg) || delta>5)
                    continue;

                // skip if gap is large
                if (GenericBamAlignmentTools::maxGapLength(alignRead, alignGenome)>5)
                    continue;

                // ease '-' at the end
                for (int t=alignRead.length()-1; t>=0; --t)
                {
                    if (alignRead[t]!=Spa && alignGenome[t]!=Spa)
                        break;

                    alignRead.erase(t, 1);
                    alignGenome.erase(t, 1);
                }

                // save
                alignReadSeqs.push_back(alignRead);
                alignGenomeSeqs.push_back(alignGenome);
            }
        }

        if (alignReadSeqs.size()>=m_numTrainData)
            break;
    }


    // erase items more than the specified train data size
    if (alignReadSeqs.size()>m_numTrainData)
    {
        vector<int> itemIndex;
        for (int i=0; i<alignReadSeqs.size(); ++i)
        {
            itemIndex.push_back(i);
        }

        // random shuffle
        shuffle(itemIndex.begin(), itemIndex.end(), UniformGenerator);

        // temp vector
        vector<string> tempAlignReadSeqs, tempAlignGenomeSeqs;

        for (int i=0; i<m_numTrainData; ++i)
        {
            tempAlignReadSeqs.push_back(alignReadSeqs[itemIndex[i]]);
            tempAlignGenomeSeqs.push_back(alignGenomeSeqs[itemIndex[i]]);
        }

        // save
        alignReadSeqs.clear();
        alignGenomeSeqs.clear();

        alignReadSeqs   = tempAlignReadSeqs;
        alignGenomeSeqs = tempAlignGenomeSeqs;
    }

    // -----------------------------------------------------
    // training
    if (HasMaxIter)
        m_maxIterNum = stoi(MaxIterStr);

    if (HasStopVal)
        m_stopVal = stod(StopValStr);


    if (m_verboseLevel>=0)
    {
        cout << "Datasize:"  << m_numTrainData                           << "   ";
        cout << "MaxIter:"   << m_maxIterNum                             << "   ";
        cout << "StopValue:" << m_stopVal                                << "   ";
        cout << "Method:"    << ((USE_NLOPT_LBFGS) ? "L-BFGS" : "SLSQP") << "   ";
        cout << endl;
    }

    TrainSetting settings;
    settings.m_band      = m_band;
    settings.m_maxIter   = m_maxIterNum;
    settings.m_stopVal   = m_stopVal;
    settings.m_verbosity = m_verboseLevel;

    ProbAligner.train(alignReadSeqs, alignGenomeSeqs, settings);

    // -----------------------------------------------------
    // output

    ofstream output(ConfigFile.c_str());
    output << ProbAligner << endl;
    output.close();

    // -----------------------------------------------------
    // close BAM and Fasta reader
    m_bamReader.Close();
    m_fasta.Close();

    return 1;
}


//---------------------------------------------------------------------------
// Alignment of B to A

int ProbabilisticAlignmentTool::Alignment()
{
    // make sure the config file is provide
    if (!HasConfig)
    {
        cerr << "GenericSequenceTools paln ERROR: probabilistic model config file not provided... Aborting" << endl;
        return 0;
    }
    ifstream inputConfig;
    inputConfig.open(ConfigFile);
    inputConfig >> ProbAligner;
    inputConfig.close();

    // compute alignment
    string t_alnA, t_alnB;
    double t_probScore = ProbAligner.ViterbiComputation(m_seqB, m_seqA, m_band, t_alnB, t_alnA);

    cout << "> ProbScore=" << t_probScore << endl;
    cout << "SeqB: " << t_alnB << endl;
    cout << "SeqA: " << t_alnA << endl;

    return 1;
}


//---------------------------------------------------------------------------
// declaration

struct Realigner_Sequence_t
{

    void init()
    {
        t_length = m_sequence.length();

        GenericFastaTools::markPrecedeHomopolymer(m_sequence, m_precedeHomopolymerLength);
        GenericFastaTools::markLeftNeighborHomopolymer(m_sequence, m_leftContextLength);
        GenericFastaTools::markRightNeighborHomopolymer(m_sequence, m_rightContextLength);

        update();
    }
    void update()
    {
        int readPointer=-1;
        for (int i=0; i<m_sequence.length(); ++i)
        {
            if (m_sequence[i]!=Spa)
                readPointer++;
            m_alnIdxToSeqIdx.push_back(readPointer);
        }
    }

    char& operator[](int i)
    {
        return m_sequence[i];
    }

    string          m_ID;
    string          m_sequence;
    vector<int>     m_alnIdxToSeqIdx;
    vector<int>     m_precedeHomopolymerLength;
    vector<int>     m_leftContextLength;
    vector<int>     m_rightContextLength;
    int             m_startPositionShift;
    int             m_endPositionShift;
    map<int,string> m_inserts;

    string          t_sequence;
    Cigar           t_cigar;
    BamMD           t_md;
    int             t_numMismatch;
    int             t_numInDel;
    int             t_length;
    int             t_type;
    string          t_alnLocalRead;
    string          t_alnLocalGenome;
    int             t_mapQualScore;
    bool            t_realign;


};

struct Realigner_Profile_t
{
    int m_length;
    int m_dimension;
    vector<Realigner_Sequence_t> m_profileAlignSequences;
};

struct Realigner_Sequence_Sort_Rule
{
    bool operator()(const Realigner_Sequence_t& a, const Realigner_Sequence_t& b)
    {
        bool a_before_b = false;

        int mag = GenericBamAlignmentTools::maxGapLength((Cigar&)a.t_cigar);
        int mbg = GenericBamAlignmentTools::maxGapLength((Cigar&)b.t_cigar);

        if (a.t_type==READ_CROSS_REGION && b.t_type==READ_CROSS_REGION)
        {
            if (a.t_numInDel<b.t_numInDel)
            {
                return true;
            }

            if (a.t_numInDel==b.t_numInDel)
            {
                if (mag<mbg)
                    return false;

                if (mag==mbg)
                {
                    if (a.t_numMismatch<b.t_numMismatch)
                        return true;
                }
            }
        }
        if (a.t_type==READ_CROSS_REGION && b.t_type!=READ_CROSS_REGION)
        {
            return false;
        }
        if (a.t_type!=READ_CROSS_REGION && b.t_type==READ_CROSS_REGION)
        {
            return true;
        }
        if (a.t_type!=READ_CROSS_REGION && b.t_type!=READ_CROSS_REGION)
        {

            if (a.t_numInDel<b.t_numInDel)
            {
                return true;
            }

            if (a.t_numInDel==b.t_numInDel)
            {
                if (mag<mbg)
                    return false;

                if (mag==mbg)
                {
                    if (a.t_numMismatch<b.t_numMismatch)
                        return true;
                }
            }
        }

        return a_before_b;
    }
};


// re-alignment interface
int windowRealignment(GenericProbabilisticAlignment& ProbAligner, BamReader& BamReader, Fasta& fasta,
                      int leftRefID, int leftRefPos, int rightRefID, int rightRefPos,
                      unordered_map<string,BamAlignment>& alignObjects, bool hasDownsample, int downsample,
                      int graphPruneLevel, int topK, bool haplotypeLikeRank, bool updateGraphByNewAlignment, bool procLargeHomopolymer,
                      bool keepDuplicate, int minReadLen, int minMapQual, double maxMismatchFrac);

// find locus to be re-aligned
int windowRealignSites(BamReader& bamReader, Fasta& fasta,
                       int leftRefID, int leftRefPos, int rightRefID, int rightRefPos,
                       vector<int>& realignStartPos, vector<int>& realignEndPos, int flankSize, bool keepDuplicate,
                       int minRepeatLength, int minReadLen, int minMapQual, double maxMismatchFrac, int minSnpRead, int minInDelRead, bool extendFlankSize=true);

// sort and merge realigning loci
int windowRealignSitesMerge(vector<int>& realignStartPos, vector<int>& realignEndPos, vector<int>& mergeRealignStartPos, vector<int>& mergeRealignEndPos);

// SNPs in window
int windowSnpSites(BamReader& bamReader,
                   int leftRefID, int leftRefPos, int rightRefID, int rightRefPos,
                   map<int,int>& snpSites, bool keepDuplicate,
                   int minReadLen, int minMapQual, double maxMismatchFrac);
void windowSnpSitesFilter(map<int,int>& snpSites, map<int,int>& passSnpSites, int minSnpRead);

// Deletes in window
int windowDeleteSites(BamReader& bamReader,
                      int leftRefID, int leftRefPos, int rightRefID, int rightRefPos,
                      map<int,int>& deleteSites, bool keepDuplicate,
                      int minReadLen, int minMapQual, double maxMismatchFrac);
void windowDeleteSitesFilter(map<int,int>& deleteSites, map<int,int>& passDeleteSites, int minIndelRead);

// SNPs in tandem repeat
int windowSnpInTandemRepeat(int tandemRepeatStartPos, int tandemRepeatEndPos, map<int,int>& snpInWindow, map<int,int>::iterator& iter);

// search MNP in a region
void windowMnpSites(BamReader& bamReader, int leftRefID, int leftRefPos, int rightRefID, int rightRefPos,
                    vector<int>& mnpStartPos, vector<int>& mnpEndPos, int flankSize,
                    bool keepDuplicate, int minReadLen, int minMapQual, double maxMismatchFrac, int minSnpRead);

// class function implementation

int ProbabilisticAlignmentTool::Realignment()
{
    // make sure the output file is provided
    if (!HasOutput)
    {
        cerr << "GenericSequenceTools paln ERROR: BAM output file not provided... Aborting" << endl;
        return 0;
    }

    // make sure the config file is provide
    if (!HasConfig)
    {
        cerr << "GenericSequenceTools paln ERROR: probabilistic model config file not provided... Aborting" << endl;
        return 0;
    }
    ifstream inputConfig;
    inputConfig.open(ConfigFile);
    inputConfig >> ProbAligner;
    inputConfig.close();

    // open file for bam writing
    m_bamWriter.Open(OutputFile, m_bamReader.GetHeaderText(), m_bamGenomes);

    if (HasCompress)
        m_bamWriter.SetCompressionMode(BamWriter::Compressed);

    // region size is 1000000 bp
    int t_RegionSize = 1e6;
    int t_LeftRefID, t_RightRefID;
    int t_LeftRefPos, t_RightRefPos;

    unordered_map<string,BamAlignment> realignObjects;
    int numRealignRegion = 0;
    int numRegion = 0;
    // iterate through small region
    if (HasRegion)
    {
        t_LeftRefID  = m_bamRegion.LeftRefID;
        t_RightRefID = m_bamRegion.RightRefID;

        int t_lastEndPos=-1;
        // iterate over regions
        for (t_LeftRefPos=m_bamRegion.LeftPosition; t_LeftRefPos<m_bamRegion.RightPosition; t_LeftRefPos+=t_RegionSize)
        {

            clock_t tBegin = clock();

            if (t_LeftRefPos+t_RegionSize<m_bamRegion.RightPosition)
                t_RightRefPos = t_LeftRefPos+t_RegionSize;
            else
                t_RightRefPos = m_bamRegion.RightPosition;

            // collect possible realignment sites in window
            vector<int> realnStartPos, realnEndPos;
            windowRealignSites(m_bamReader, m_fasta, t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos, realnStartPos, realnEndPos, m_flankSize, HasKeepDuplicate, m_repeatLength, m_minLength, m_minMapQual, m_maxMismatches, m_minSnpRead, m_minInDelRead);

            if (HasMNP)
            {
                vector<int> mnpStartPos, mnpEndPos;
                windowMnpSites(m_bamReader, t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos, mnpStartPos, mnpEndPos, m_flankSize, HasKeepDuplicate, m_minLength, m_minMapQual, m_maxMismatches, m_minSnpRead);
                for (int i=0; i<mnpStartPos.size(); i++)
                {
                    realnStartPos.push_back(mnpStartPos[i]);
                    realnEndPos.push_back(mnpEndPos[i]);
                }
            }

            vector<int> mergeRealnStartPos, mergeRealnEndPos;
            windowRealignSitesMerge(realnStartPos, realnEndPos, mergeRealnStartPos, mergeRealnEndPos);
            realnStartPos = mergeRealnStartPos;
            realnEndPos   = mergeRealnEndPos;

            // realignment in each site
            for (int i=0; i<realnStartPos.size(); i++)
            {
                int t_startPos, t_endPos;
                t_startPos = realnStartPos[i];
                t_endPos   = realnEndPos[i];


//                if (t_startPos<t_lastEndPos)
//                    t_startPos = t_lastEndPos+1;

                if (m_verboseLevel>=2)
                {
                    cout << m_bamGenomes[t_LeftRefID].RefName << "\t"
                         << t_startPos+1 << "\t"
                         << t_endPos+1+1 << endl;
                }

                int RealignSuccess = windowRealignment(ProbAligner, m_bamReader, m_fasta, t_LeftRefID, t_startPos, t_LeftRefID, t_endPos+1, realignObjects, HasDownSample, m_downSample, m_graphEdgePruneLevel, m_topK, HasHaplotypeLikeRank, HasUpdateGraph, HasProcLargeHomopolymer, HasKeepDuplicate, m_minLength, m_minMapQual, m_maxMismatches);

                numRealignRegion += RealignSuccess;

                numRegion ++;

                t_lastEndPos = t_endPos+1;
            }

            // reader rewind
            m_bamReader.Rewind();
            m_bamReader.SetRegion(t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos);

            BamAlignment al;
            while (m_bamReader.GetNextAlignment(al))
            {
                unordered_map<string,BamAlignment>::iterator iter = realignObjects.find(GenericBamAlignmentTools::getBamAlignmentID(al));
                if (iter==realignObjects.end())
                {
                    m_bamWriter.SaveAlignment(al);
                }
                else
                {
                    m_bamWriter.SaveAlignment(iter->second);
                }
            }

            clock_t tEnd = clock();

            if (m_verboseLevel>=1)
            {
                cout << "time elapsed " << double(tEnd-tBegin)/CLOCKS_PER_SEC/60.0 << " minutes"
                     << ", process " << numRealignRegion << " of " << numRegion << " regions"
                     << ", and "
                     << "re-align " << realignObjects.size() << " reads." << endl;
            }
        }
    }else
    {
        // iterate over genomes
        for (int i=0; i<m_bamGenomes.size(); ++i)
        {

            t_LeftRefID  = i;
            t_RightRefID = i;

            int t_lastEndPos = -1;
            // iterate over regions
            for (t_LeftRefPos=0; t_LeftRefPos<m_bamGenomes[i].RefLength; t_LeftRefPos+=t_RegionSize)
            {
                clock_t tBegin = clock();

                if (t_LeftRefPos+t_RegionSize<m_bamGenomes[i].RefLength)
                    t_RightRefPos = t_LeftRefPos+t_RegionSize;
                else
                    t_RightRefPos = m_bamGenomes[i].RefLength;

                // collect possible realignment sites in window
                vector<int> realnStartPos, realnEndPos;
                windowRealignSites(m_bamReader, m_fasta, t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos, realnStartPos, realnEndPos, m_flankSize, HasKeepDuplicate, m_repeatLength, m_minLength, m_minMapQual, m_maxMismatches, m_minSnpRead, m_minInDelRead);

                if (HasMNP)
                {
                    vector<int> mnpStartPos, mnpEndPos;
                    windowMnpSites(m_bamReader, t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos, mnpStartPos, mnpEndPos, m_flankSize, HasKeepDuplicate, m_minLength, m_minMapQual, m_maxMismatches, m_minSnpRead);
                    for (int i=0; i<mnpStartPos.size(); i++)
                    {
                        realnStartPos.push_back(mnpStartPos[i]);
                        realnEndPos.push_back(mnpEndPos[i]);
                    }
                }

                vector<int> mergeRealnStartPos, mergeRealnEndPos;
                windowRealignSitesMerge(realnStartPos, realnEndPos, mergeRealnStartPos, mergeRealnEndPos);
                realnStartPos = mergeRealnStartPos;
                realnEndPos   = mergeRealnEndPos;

                // realignment in each site
                for (int i=0; i<realnStartPos.size(); i++)
                {
                    int t_startPos, t_endPos;
                    t_startPos = realnStartPos[i];
                    t_endPos   = realnEndPos[i];

//                    if (t_startPos<t_lastEndPos)
//                        t_startPos = t_lastEndPos+1;

                    if (m_verboseLevel>=2)
                    {
                        cout << m_bamGenomes[t_LeftRefID].RefName << "\t"
                             << t_startPos+1 << "\t"
                             << t_endPos+1+1 << endl;
                    }

                    numRealignRegion += windowRealignment(ProbAligner, m_bamReader, m_fasta, t_LeftRefID, t_startPos, t_LeftRefID, t_endPos+1, realignObjects, HasDownSample, m_downSample, m_graphEdgePruneLevel, m_topK, HasHaplotypeLikeRank, HasUpdateGraph, HasProcLargeHomopolymer, HasKeepDuplicate, m_minLength, m_minMapQual, m_maxMismatches);
                    numRegion ++;

                    t_lastEndPos = t_endPos+1;
                }

                // reader rewind
                m_bamReader.Rewind();
                m_bamReader.SetRegion(t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos);

                BamAlignment al;
                while (m_bamReader.GetNextAlignment(al))
                {
                    unordered_map<string,BamAlignment>::iterator iter = realignObjects.find(GenericBamAlignmentTools::getBamAlignmentID(al));
                    if (iter==realignObjects.end())
                    {
                        m_bamWriter.SaveAlignment(al);
                    }
                    else
                    {
                        m_bamWriter.SaveAlignment(iter->second);
                    }
                }

                clock_t tEnd = clock();

                if (m_verboseLevel>=1)
                {
                    cout << "time elapsed " << double(tEnd-tBegin)/CLOCKS_PER_SEC/60.0 << " minutes"
                         << ", process " << numRealignRegion << " of " << numRegion << " regions"
                         <<", and "
                         << "re-align " << realignObjects.size() << " reads." << endl;
                }
            }
        }
    }

    // processing information
    cout << "Totally re-align " << numRealignRegion << " of " << numRegion << " regions"
         <<", and "
         << "process " << realignObjects.size() << " reads." << endl;

    // close reader, wirter, and fasta
    m_bamReader.Close();
    m_fasta.Close();
    m_bamWriter.Close();


    return 1;
}






//--------------------------
// implementation

int windowRealignSites(BamReader &bamReader, Fasta &fasta,
                       int leftRefID, int leftRefPos, int rightRefID, int rightRefPos,
                       vector<int>& realignStartPos, vector<int>& realignEndPos, int flankSize,
                       bool keepDuplicate, int minRepeatLength, int minReadLen, int minMapQual,
                       double maxMismatchFrac, int minSnpRead, int minInDelRead, bool extendFlankSize)
{
    int t_nN = 0;
    int t_headShift = 0;
    string t_fastaSeq;
    fasta.GetSequence(leftRefID, leftRefPos, rightRefPos, t_fastaSeq);

    // remove head N
    for (string::iterator iter=t_fastaSeq.begin(); iter!=t_fastaSeq.end(); ++iter)
    {
        if ((*iter)==Amb)
        {
            t_nN++;
            t_headShift++;
            leftRefPos++;
        }
        else
            break;
    }

    // remove tail N
    for (string::reverse_iterator iter=t_fastaSeq.rbegin(); iter!=t_fastaSeq.rend(); ++iter)
    {
        if ((*iter)==Amb)
        {
            t_nN++;
            rightRefPos--;
        }
        else
            break;
    }

    // window is N
    if (leftRefPos>=rightRefPos)
        return 0;

    if (t_nN>0)
    {
        t_fastaSeq = t_fastaSeq.substr(t_headShift, rightRefPos-leftRefPos);
    }

    // precede homopolymer length
    vector<int> precedeHomopolymerLength;
    GenericFastaTools::markPrecedeHomopolymer(t_fastaSeq, precedeHomopolymerLength);
    // left neighbor homopolymer length
    vector<int> leftNeighborHomopolymerLength;
    GenericFastaTools::markLeftNeighborHomopolymer(t_fastaSeq, leftNeighborHomopolymerLength);
    // succede homopolymer length
    vector<int> succedeHomopolymerLength;
    GenericFastaTools::markSuccedeHomopolymer(t_fastaSeq, succedeHomopolymerLength);
    // right neighbor homopolymer length
    vector<int> rightNeighborHomopolymerLength;
    GenericFastaTools::markRightNeighborHomopolymer(t_fastaSeq, rightNeighborHomopolymerLength);

    // search SNP sites in window
    map<int,int> snpInWindow;
    map<int,int>::iterator iterSnpInWindow;
    windowSnpSites(bamReader, leftRefID, leftRefPos, rightRefID, rightRefPos, snpInWindow, keepDuplicate, minReadLen, minMapQual, maxMismatchFrac);

    map<int,int> passSnpInWindow;
    windowSnpSitesFilter(snpInWindow, passSnpInWindow, minSnpRead);

    // search Delete sites in window
    map<int,int> deleteInWindow;
    windowDeleteSites(bamReader, leftRefID, leftRefPos, rightRefID, rightRefPos, deleteInWindow, keepDuplicate, minReadLen, minMapQual, maxMismatchFrac);

    map<int,int> passDeleteInWindow;
    windowDeleteSitesFilter(deleteInWindow, passDeleteInWindow, minInDelRead);

    // search tandem repeat sites in window
    vector<int> tandemRepeatPos, tandemRepeatSize;
    GenericFastaTools::findTandemRepeatRegions(t_fastaSeq, minRepeatLength, tandemRepeatPos, tandemRepeatSize);

    // iterate from first repeat to last repeat
    iterSnpInWindow = snpInWindow.begin();
    for (int i=0; i<tandemRepeatPos.size(); ++i)
    {
        int t_leftPos, t_rightPos;

        t_leftPos  = leftRefPos+tandemRepeatPos[i]-flankSize;
        t_rightPos = t_leftPos+tandemRepeatSize[i]+flankSize;

        if (extendFlankSize)
        {
            // enlong left position
            do{
                t_leftPos--;
            }while(t_leftPos>leftRefPos && (passSnpInWindow.count(t_leftPos)>0 || passDeleteInWindow.count(t_leftPos)>0 || precedeHomopolymerLength[t_leftPos-leftRefPos]>0 || leftNeighborHomopolymerLength[t_leftPos-leftRefPos]>2));

            // enlong right position
            do{
                t_rightPos++;
            }while (t_rightPos<rightRefPos && (passSnpInWindow.count(t_rightPos)>0 || passDeleteInWindow.count(t_rightPos)>0 || succedeHomopolymerLength[t_rightPos-leftRefPos]>0 || rightNeighborHomopolymerLength[t_rightPos-leftRefPos]>2));
        }

        // count SNPs in repeat region
        if (windowSnpInTandemRepeat(t_leftPos, t_rightPos, snpInWindow, iterSnpInWindow)>0)
        {
            realignStartPos.push_back(t_leftPos);
            realignEndPos.push_back(t_rightPos);
        }
    }

    // return the number of regions need realignment
    return realignStartPos.size();
}


struct RealignSiteSortByPosition
{
    bool operator()(const tuple<int,int>& a, const tuple<int,int>& b)
    {
        return (get<0>(a) < get<0>(b));
    }
};

int windowRealignSitesMerge(vector<int> &realignStartPos, vector<int> &realignEndPos, vector<int> &mergeRealignStartPos, vector<int> &mergeRealignEndPos)
{
    vector<tuple<int,int>> realignSites;
    for (int i=0; i<realignStartPos.size(); i++)
    {
        realignSites.push_back(tuple<int,int>(realignStartPos[i], realignEndPos[i]));
    }

    // sort by position
    sort(realignSites.begin(), realignSites.end(), RealignSiteSortByPosition());

    // merge
    for (int i=0; i<realignSites.size(); i++)
    {
        if (i==0)
        {
            mergeRealignStartPos.push_back(get<0>(realignSites[i]));
            mergeRealignEndPos.push_back(get<1>(realignSites[i]));
            continue;
        }

        int j = mergeRealignStartPos.size()-1;
        int t0 = mergeRealignStartPos[j];
        int t1 = mergeRealignEndPos[j];
        int s0 = get<0>(realignSites[i]);
        int s1 = get<1>(realignSites[i]);

        int tl = t1-t0;

        if (t0<=s0 && t1>=s1)
            continue;

        if (mergeRealignEndPos[j]+2>=get<0>(realignSites[i]) && tl<30)
        {
            if (mergeRealignEndPos[j]<get<1>(realignSites[i]))
                mergeRealignEndPos[j] = get<1>(realignSites[i]);
        }else
        {
            mergeRealignStartPos.push_back(get<0>(realignSites[i]));
            mergeRealignEndPos.push_back(get<1>(realignSites[i]));
        }
    }

    return mergeRealignStartPos.size();
}

int windowSnpSites(BamReader& bamReader,
                   int leftRefID, int leftRefPos, int rightRefID, int rightRefPos,
                   map<int,int>& snpSites,
                   bool keepDuplicate, int minReadLen, int minMapQual, double maxMismatchFrac)
{
    // rewind
    bamReader.Rewind();

    // set region
    bamReader.SetRegion(leftRefID, leftRefPos, rightRefID, rightRefPos);

    // iteratively read alignment
    BamAlignment al;
    while(bamReader.GetNextAlignment(al))
    {
        // skip if it is not good alignment
        if (!GenericBamAlignmentTools::goodAlignment(al, keepDuplicate))
        {
            continue;
        }

        // skip if there is no MD tag
        if (!al.HasTag("MD"))
            continue;

        // skip if read is too short
        int len = GenericBamAlignmentTools::getBamAlignmentReadLength(al);
        if (len<minReadLen)
        {
            continue;
        }

        // skip if map quality is poor
        if (al.MapQuality<minMapQual)
        {
            continue;
        }

        // skip if mismatches is too much
        if (GenericBamAlignmentTools::numBamAlignmentMismatches(al)>maxMismatchFrac*len)
        {
            continue;
        }

        // retrieve SNPs in the read
        vector<long> snpInRead;
        GenericBamAlignmentTools::getBamAlignmentMismatches(al, snpInRead);

        if (snpInRead.empty())
            continue;

        // save SNPs
        for (vector<long>::iterator pos=snpInRead.begin(); pos!=snpInRead.end(); ++pos)
        {
            if ((*pos)<leftRefPos || (*pos)>=rightRefPos)
                continue;

            if (snpSites.find((*pos))==snpSites.end())
            {
                snpSites[(*pos)] = 1;
            }else
            {
                snpSites[(*pos)] += 1;
            }
        }
    }

    // rewind
    bamReader.Rewind();

    return snpSites.size();
}


int windowDeleteSites(BamReader& bamReader,
                      int leftRefID, int leftRefPos, int rightRefID, int rightRefPos,
                      map<int,int>& deleteSites,
                      bool keepDuplicate, int minReadLen, int minMapQual, double maxMismatchFrac)
{
    // rewind
    bamReader.Rewind();

    // set region
    bamReader.SetRegion(leftRefID, leftRefPos, rightRefID, rightRefPos);

    // iteratively read alignment
    BamAlignment al;
    while(bamReader.GetNextAlignment(al))
    {
        // skip if it is not good alignment
        if (!GenericBamAlignmentTools::goodAlignment(al, keepDuplicate))
            continue;

        // skip if there is no MD tag
        if (!al.HasTag("MD"))
            continue;

        // skip if read is too short
        int len = GenericBamAlignmentTools::getBamAlignmentReadLength(al);
        if (len<minReadLen)
            continue;

        // skip if map quality is poor
        if (al.MapQuality<minMapQual)
            continue;

        // skip if mismatches is too much
        if (GenericBamAlignmentTools::numBamAlignmentMismatches(al)>maxMismatchFrac*len)
            continue;

        // retrieve SNPs in the read
        vector<long> deleteInRead;
        vector<string> deleteStringInRead;
        GenericBamAlignmentTools::getBamAlignmentDeletes(al, deleteInRead, deleteStringInRead);

        if (deleteInRead.empty())
            continue;

        // save deletes
        int i=0;
        for (vector<long>::iterator pos=deleteInRead.begin(); pos!=deleteInRead.end(); ++pos,++i)
        {
            if ((*pos)<leftRefPos || (*pos)>=rightRefPos)
                continue;

            for (int j=0; j<deleteStringInRead[i].length(); j++)
            {
                if (deleteSites.find((*pos)+j)==deleteSites.end())
                {
                    deleteSites[(*pos)+j] = 1;
                }else
                {
                    deleteSites[(*pos)+j] += 1;
                }
            }
        }
    }

    // rewind
    bamReader.Rewind();

    return deleteSites.size();
}

void windowSnpSitesFilter(map<int,int>& snpSites, map<int,int>& passSnpSites, int minSnpRead)
{
    for (map<int,int>::iterator iter=snpSites.begin(); iter!=snpSites.end(); iter++)
    {
        if (iter->second>=minSnpRead)
            passSnpSites[iter->first] = iter->second;
    }
}

void windowDeleteSitesFilter(map<int,int>& deleteSites, map<int,int>& passDeleteSites, int minIndelRead)
{
    for (map<int,int>::iterator iter=deleteSites.begin(); iter!=deleteSites.end(); iter++)
    {
        if (iter->second>=minIndelRead)
            passDeleteSites[iter->first] = iter->second;
    }
}

void smallGapRemove(string& i_seq, Cigar& i_cigar, string& o_seq, Cigar& o_cigar)
{
    o_seq = "";
    int seqPointer = 0;
    for (Cigar::iterator iter=i_cigar.begin(); iter!=i_cigar.end(); ++iter)
    {
        if (GenericBamAlignmentTools::isMatch(iter->Type) || GenericBamAlignmentTools::isMismatch(iter->Type))
        {
            for (int i=0; i<iter->Length; ++i, ++seqPointer)
            {
                o_seq += i_seq[seqPointer];
            }
            o_cigar.push_back(*iter);
        }

        if (GenericBamAlignmentTools::isInsert(iter->Type))
        {
            if (iter->Length==1)
            {
                ++seqPointer;
            }else
            {
                for (int i=0; i<iter->Length; ++i, ++seqPointer)
                {
                    o_seq += i_seq[seqPointer];
                }
                o_cigar.push_back(*iter);
            }
        }

        if (GenericBamAlignmentTools::isDelete(iter->Type))
        {
            o_cigar.push_back(*iter);
        }
    }
}

struct Sort_Haplotype_By_LOR
{
    bool operator()(const tuple<double,int>& a, const tuple<double,int>& b)
    {
        return (get<0>(a)>get<0>(b));
    }
};

void rankHaplotypeByLikelihood(GenericProbabilisticAlignment& ProbAligner, string& fastaSeq, vector<string>& haploSeqs, vector<Realigner_Sequence_t>& readSeqs, vector<int>& topRankHaplotypes)
{
    double ll;
    // number of reads
    int n = readSeqs.size();

    // likelihood of reference
    double likeli0 = 0;
    vector<double> likeli0_datum(n, 0);
    for (int i=0; i<n; ++i)
    {
        ll = ProbAligner.ViterbiComputation(readSeqs[i].t_sequence, fastaSeq, 0.1);
        likeli0_datum[i] = ll;
        likeli0 += log(ll);
    }
    likeli0 = exp(likeli0);

    // number of alternative haplotypes
    int m = haploSeqs.size();
    // likelihood of haplotypes
    vector<double> likeli1(m, 0);
    vector<vector<double>> likeli1_datum(m, vector<double>(n, 0));
    for (int j=0; j<m; j++)
    {
        string haploSeq = haploSeqs[j];
        for (int i=0; i<n; i++)
        {
            ll = ProbAligner.ViterbiComputation(readSeqs[i].t_sequence, haploSeq);
            likeli1_datum[j][i] = ll;
            likeli1[j] += log(0.5*ll+0.5*likeli0_datum[i]);
        }
        likeli1[j] = exp(likeli1[j]);
    }

    // log odds ratio
    vector<double> lor(m, 0);
    for (int j=0; j<m; ++j)
    {
        lor[j] = log(likeli1[j]) - log(likeli0);
    }

    // sort haplotype by log odds ratio
    vector<tuple<double,int>> haploLors;
    for (int j=0; j<m; j++)
    {
        haploLors.push_back(tuple<double,int>(lor[j],j));
    }
    sort(haploLors.begin(), haploLors.end(), Sort_Haplotype_By_LOR());

    // picked haplotypes
    for (int j=0; j<m; j++)
    {
        cout << get<0>(haploLors[j]) << endl;

//        if (get<0>(haploLors[j])>0.3)
        {
            topRankHaplotypes.push_back(get<1>(haploLors[j]));
        }
    }
}

int windowRealignment(GenericProbabilisticAlignment &ProbAligner, BamReader &BamReader, Fasta &fasta,
                      int leftRefID, int leftRefPos, int rightRefID, int rightRefPos,
                      unordered_map<string, BamAlignment> &alignObjects, bool hasDownsample, int downsample,
                      int graphPruneLevel, int topK, bool haplotypeLikeRank, bool updateGraphByNewAlignment, bool procLargeHomopolymer,
                      bool keepDuplicate, int minReadLen, int minMapQual, double maxMismatchFrac)
{
     BamAlignment al;

    // number of alignments in the region
    int numAlignTotal = 0;
    if (hasDownsample)
    {
        // rewind BAM reader
        BamReader.Rewind();
        // set BAM region
        BamReader.SetRegion(leftRefID, leftRefPos, rightRefID, rightRefPos);
        while(BamReader.GetNextAlignment(al))
        {
            // skip if alignment is not in the region
            if (al.Position>=rightRefPos || al.GetEndPosition()<=leftRefPos)
                continue;
            // skip if read fail quality assessments
            if (readIsFiltered(al, keepDuplicate, minReadLen, minMapQual, maxMismatchFrac))
                continue;

            // counting
            numAlignTotal++;
        }
    }

    // random dicer
    mt19937_64 UniformGenerator;
    uniform_int_distribution<int> UniformDist;
    if (hasDownsample)
    {
        unsigned UniformSeed = chrono::system_clock::now().time_since_epoch().count();
        UniformGenerator = mt19937_64(UniformSeed);
        UniformDist = uniform_int_distribution<int>(0, numAlignTotal);
    }
    auto UniformDice = bind(UniformDist, UniformGenerator);

    // genome sequence
    string t_fastaSeq;
    fasta.GetSequence(leftRefID, leftRefPos, rightRefPos-1, t_fastaSeq);

    // set of sequences in local region to be re-aligned
    vector<Realigner_Sequence_t> t_seqInLocalRegion;

    // rewind BAM reader
    BamReader.Rewind();

    // set BAM region
    BamReader.SetRegion(leftRefID, leftRefPos, rightRefID, rightRefPos);

    int numMisAligns = 0;
    // retrieve alignment
    while(BamReader.GetNextAlignment(al))
    {
        BamAlignment alignObj;
        // skip if alignment is not in the region
        if (al.Position>=rightRefPos || al.GetEndPosition()<=leftRefPos)
            continue;

        // skip if read fail quality assessments
        if (readIsFiltered(al, keepDuplicate, minReadLen, minMapQual, maxMismatchFrac))
            continue;

        // skip if read is not sampled
        if (hasDownsample)
        {
            if (UniformDice()>downsample)
                continue;
        }

        // check whether it is in realigned set
        string t_ID = GenericBamAlignmentTools::getBamAlignmentID(al);
        unordered_map<string,BamAlignment>::iterator iter = alignObjects.find(t_ID);

        if (iter!=alignObjects.end())
        {
            alignObj = BamAlignment(iter->second);
        }else
        {
            alignObj = BamAlignment(al);
        }

        // extract local sequences
        string t_localRead, t_localGenome;
        Cigar  t_cigar;
        BamMD  t_md;
        int    t_numMismatch, t_numInDel;
        GenericBamAlignmentTools::getLocalAlignment(alignObj, leftRefPos, rightRefPos-leftRefPos,
                                                    t_localRead, t_localGenome, t_cigar, t_md,
                                                    t_numMismatch, t_numInDel);

        if (t_localRead.empty() || t_localGenome.empty())
            continue;

        // put into set
        Realigner_Sequence_t t_seq;
        t_seq.m_ID           = GenericBamAlignmentTools::getBamAlignmentID(alignObj);
        t_seq.m_sequence     = t_localRead;
        t_seq.t_sequence     = t_localRead;
        t_seq.t_cigar        = t_cigar;
        t_seq.t_md           = t_md;
        t_seq.t_numMismatch  = t_numMismatch;
        t_seq.t_numInDel     = t_numInDel;
        t_seq.t_mapQualScore = alignObj.MapQuality;
        t_seq.t_realign      = false;


        if (alignObj.Position>leftRefPos)
            t_seq.m_startPositionShift = alignObj.Position-leftRefPos;
        else
            t_seq.m_startPositionShift = 0;

        if (alignObj.GetEndPosition()<rightRefPos)
            t_seq.m_endPositionShift = rightRefPos-alignObj.GetEndPosition();
        else
            t_seq.m_endPositionShift = 0;

        t_seq.t_type = READ_CROSS_REGION;
        if (t_seq.m_startPositionShift>0)
            t_seq.t_type = READ_BEGIN_REGION;
        if (t_seq.m_endPositionShift>0)
            t_seq.t_type = READ_END_REGION;

        t_seq.init();

        t_seqInLocalRegion.push_back(t_seq);

        if (t_numMismatch>0)
            numMisAligns += 1;
    }

    if (numMisAligns<=1)
        return 0;

    // sort reads
    sort(t_seqInLocalRegion.begin(), t_seqInLocalRegion.end(), Realigner_Sequence_Sort_Rule());

    // construct the graph to check the alignment consistency
    GenericDagGraph alignGraph;
    vector<string>  alignGraphReads;
    vector<Cigar>   alignGraphReadCigars;
    vector<int>     alignGraphReadStarts;

    // set of aligned reads to construct the graph
    for (int i=0; i<t_seqInLocalRegion.size(); ++i)
    {
        alignGraphReads.push_back(t_seqInLocalRegion[i].t_sequence);
        alignGraphReadCigars.push_back(t_seqInLocalRegion[i].t_cigar);
        alignGraphReadStarts.push_back(t_seqInLocalRegion[i].m_startPositionShift);
    }
    // construct the alignment graph
    alignGraph.buildDagGraph(t_fastaSeq, alignGraphReads, alignGraphReadCigars, alignGraphReadStarts);
    if (graphPruneLevel>0)
    {
        alignGraph.edgePruning2(graphPruneLevel-1, graphPruneLevel);
//        alignGraph.edgePruning(graphPruneLevel);
    }

    // debug
//    cout << alignGraph << endl;

    // re-align if the graph is inconsistent
    vector<int> genomePrecedeHomopolymerLength;
    vector<int> genomeSuccedeHomopolymerLength;
    GenericFastaTools::markPrecedeHomopolymer(t_fastaSeq, genomePrecedeHomopolymerLength);
    GenericFastaTools::markSuccedeHomopolymer(t_fastaSeq, genomeSuccedeHomopolymerLength);

    if (!alignGraph.isConsistentGraph())
    {
        int maxGenomeHomopolymer = GenericFastaTools::maxHomopolymerLength(t_fastaSeq);

        GenericDagGraph haplotypeGraph;

        // pick top rank path(s) as alternative haplotype sequence(s)
        vector<string>       topRankAlignGraphPaths;
        vector<list<Vertex>> topRankAlignGraphPathVertexs;
        vector<double>       topRankAlignGraphPathWeights;
        alignGraph.topRankPathsExcludeGenome(10+topK, topRankAlignGraphPaths, topRankAlignGraphPathVertexs, topRankAlignGraphPathWeights);

        vector<int> haplotypeRank;
        if (haplotypeLikeRank)
        {
            rankHaplotypeByLikelihood(ProbAligner, t_fastaSeq, topRankAlignGraphPaths, t_seqInLocalRegion, haplotypeRank);
        }else
        {
            for (int i=0; i<topRankAlignGraphPaths.size(); i++)
                haplotypeRank.push_back(i);
        }

        vector<string>       haplotypeSequences;
        vector<Cigar>        haplotypeCigars;
        vector<int>          haplotypeStarts;
        for (int j=0; j<topK; ++j)
        {

            if (haplotypeRank.empty() || j>=haplotypeRank.size())
                continue;

            if (maxGenomeHomopolymer>8 && j>=1 && procLargeHomopolymer)
                break;

            int i = haplotypeRank[j];
            haplotypeSequences.push_back(topRankAlignGraphPaths[i]);

            Cigar h_cigar;
            alignGraph.pathCigar(topRankAlignGraphPathVertexs[i], h_cigar);
            haplotypeCigars.push_back(h_cigar);

            haplotypeStarts.push_back(0);
        }


        // construct the haplotype graph
        haplotypeGraph.buildDagGraph(t_fastaSeq, haplotypeSequences, haplotypeCigars, haplotypeStarts);

        // debug
//        cout << haplotypeGraph << endl;

        // re-align reads
        for (int i=0; i<t_seqInLocalRegion.size(); i++)
        {
            if (t_seqInLocalRegion[i].t_numMismatch==0 && t_seqInLocalRegion[i].t_numInDel==0)
                continue;

            string r_seq = t_seqInLocalRegion[i].t_sequence;

            int r_type = t_seqInLocalRegion[i].t_type;

            // realignment computation
            int d_start;
            string d_alnRead;
            string d_alnGenome;


            Cigar r_cigar;
            haplotypeGraph.alignReadToGraph(r_seq, ProbAligner, r_type, 1.0, procLargeHomopolymer, r_cigar, d_start, d_alnGenome, d_alnRead);

            int r_numMismatch = GenericBamAlignmentTools::numMismatch(d_alnRead, d_alnGenome);
            if (r_numMismatch<=t_seqInLocalRegion[i].t_numMismatch || r_numMismatch==0)
            {
                t_seqInLocalRegion[i].t_realign = true;
            }
//            // debug
//            cout << ">" << t_seqInLocalRegion[i].m_ID << endl;
//            cout << d_alnRead << endl;
//            cout << d_alnGenome << endl;

            // update the haplotype graph
//            int r_maxGap = GenericBamAlignmentTools::maxGapLength(r_cigar);
            int r_numMis = GenericBamAlignmentTools::numMismatch(d_alnRead, d_alnGenome);

            if (t_seqInLocalRegion[i].t_realign)
            {
                if (r_type==READ_CROSS_REGION && updateGraphByNewAlignment && t_seqInLocalRegion[i].t_mapQualScore>=20 && r_numMis<=5)
                {
                    int r_start = 0;
                    int r_count = 1;
                    haplotypeGraph.updateDagGraphByRead(r_seq, r_cigar, r_start, r_count, true);
                }

                // save the realignment information
                t_seqInLocalRegion[i].t_cigar = r_cigar;
                t_seqInLocalRegion[i].t_alnLocalRead   = d_alnRead;
                t_seqInLocalRegion[i].t_alnLocalGenome = d_alnGenome;
            }
        }


        // save alignments into BAM file
        map<string, int> t_seqInLocalRegionSet;
        for (int i=0; i<t_seqInLocalRegion.size(); ++i)
        {
            if (t_seqInLocalRegion[i].t_numMismatch==0 && t_seqInLocalRegion[i].t_numInDel==0)
                continue;

            t_seqInLocalRegionSet[t_seqInLocalRegion[i].m_ID] = i;
        }

        // rewind the reader
        BamReader.Rewind();

        // set the region
        BamReader.SetRegion(leftRefID, leftRefPos, rightRefID, rightRefPos);

        // retrieve alignment
        while(BamReader.GetNextAlignment(al))
        {
            // skip if alignment is not in the region
            if (al.Position>=rightRefPos || al.GetEndPosition()<=leftRefPos)
                continue;

            // skip if read fail quality assessment
            if (readIsFiltered(al, keepDuplicate, minReadLen, minMapQual, maxMismatchFrac))
                continue;

            map<string,int>::iterator iter = t_seqInLocalRegionSet.find(GenericBamAlignmentTools::getBamAlignmentID(al));
            // skip if read is not in the re-alignment set
            if (iter==t_seqInLocalRegionSet.end())
                continue;

            // skip if read is not re-aligned
            if (!t_seqInLocalRegion[iter->second].t_realign)
                continue;

            BamAlignment alignObj;
            // check whether it is realigned previously
            unordered_map<string,BamAlignment>::iterator iter2 = alignObjects.find(t_seqInLocalRegion[iter->second].m_ID);
            if (iter2!=alignObjects.end())
            {
                alignObj = BamAlignment(iter2->second);
            }else
            {
                alignObj = BamAlignment(al);
            }

            // change old alignment to new alignment
            BamAlignment n_alignObj;
            GenericBamAlignmentTools::changeLocalAlignment(alignObj, leftRefPos, rightRefPos-leftRefPos,
                                                           t_seqInLocalRegion[iter->second].t_alnLocalRead,
                                                           t_seqInLocalRegion[iter->second].t_alnLocalGenome,
                                                           n_alignObj);

            // update recording
            if (iter2==alignObjects.end())
            {
                alignObjects[t_seqInLocalRegion[iter->second].m_ID] = n_alignObj;
            }
            else
            {
                alignObjects[t_seqInLocalRegion[iter->second].m_ID] = n_alignObj;
            }
        }

    }else
    {
        return 0;
    }

    return 1;
}


int windowSnpInTandemRepeat(int tandemRepeatStartPos, int tandemRepeatEndPos, map<int, int> &snpInWindow, map<int,int>::iterator& iter)
{
    set<int> snpInTandemRepeat;

    for (; iter!=snpInWindow.end(); ++iter)
    {
        if (iter->first>=tandemRepeatStartPos && iter->first<tandemRepeatEndPos && iter->second<17)
            snpInTandemRepeat.insert(iter->first);

        if (iter->first>=tandemRepeatEndPos)
            break;
    }

    return (snpInTandemRepeat.size()-1);
}




//----------------------------------------------------------------
// Amplicon Sequencing

int windowRealignmentAmplicon(GenericProbabilisticAlignment& ProbAligner, BamReader& BamReader, Fasta& fasta,
                              int leftRefID, int leftRefPos, int rightRefID, int rightRefPos,
                              unordered_map<string,BamAlignment>& alignObjects,
                              int graphPruneLevel, int topK, bool haplotypeLikeRank, bool updateGraphByNewAlignment, bool procLargeHomopolymer,
                              bool keepDuplicate, int minReadLen=100, int minMapQual=10, double maxMismatchFrac=0.1);


int ProbabilisticAlignmentTool::RealignmentAmplicon()
{

    // make sure the output file is provided
    if (!HasOutput)
    {
        cerr << "GenericSequenceTools paln ERROR: BAM output file not provided... Aborting" << endl;
        return 0;
    }

    // make sure the config file is provide
    if (!HasConfig)
    {
        cerr << "GenericSequenceTools paln ERROR: probabilistic model config file not provided... Aborting" << endl;
        return 0;
    }
    ifstream inputConfig;
    inputConfig.open(ConfigFile);
    inputConfig >> ProbAligner;
    inputConfig.close();

    // open file for bam writing
    m_bamWriter.Open(OutputFile, m_bamReader.GetHeaderText(), m_bamGenomes);

    if (HasCompress)
        m_bamWriter.SetCompressionMode(BamWriter::Compressed);

    // region size is 1000000 bp
    int t_RegionSize = 1e6;
    int t_LeftRefID, t_RightRefID;
    int t_LeftRefPos, t_RightRefPos;

    unordered_map<string,BamAlignment> realignObjects;
    int numRealignRegion = 0;
    int numRegion = 0;
    // iterate through small region
    if (HasRegion)
    {
        t_LeftRefID  = m_bamRegion.LeftRefID;
        t_RightRefID = m_bamRegion.RightRefID;

        int t_lastEndPos=-1;
        // iterate over regions
        for (t_LeftRefPos=m_bamRegion.LeftPosition; t_LeftRefPos<m_bamRegion.RightPosition; t_LeftRefPos+=t_RegionSize)
        {
            clock_t tBegin = clock();

            if (t_LeftRefPos+t_RegionSize<m_bamRegion.RightPosition)
                t_RightRefPos = t_LeftRefPos+t_RegionSize;
            else
                t_RightRefPos = m_bamRegion.RightPosition;


            // collect possible realignment sites in window
            vector<int> realnStartPos, realnEndPos;
            windowRealignSites(m_bamReader, m_fasta, t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos, realnStartPos, realnEndPos, m_flankSize, HasKeepDuplicate, m_repeatLength, m_minLength, m_minMapQual, m_maxMismatches, m_minSnpRead, m_minInDelRead, false);

            if (HasMNP)
            {
                vector<int> mnpStartPos, mnpEndPos;
                windowMnpSites(m_bamReader, t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos, mnpStartPos, mnpEndPos, m_flankSize, HasKeepDuplicate, m_minLength, m_minMapQual, m_maxMismatches, m_minSnpRead);
                for (int i=0; i<mnpStartPos.size(); i++)
                {
                    realnStartPos.push_back(mnpStartPos[i]);
                    realnEndPos.push_back(mnpEndPos[i]);
                }
            }

            vector<int> mergeRealnStartPos, mergeRealnEndPos;
            windowRealignSitesMerge(realnStartPos, realnEndPos, mergeRealnStartPos, mergeRealnEndPos);
            realnStartPos = mergeRealnStartPos;
            realnEndPos   = mergeRealnEndPos;

            // realignment in each site
            for (int i=0; i<realnStartPos.size(); i++)
            {
                int t_startPos, t_endPos;
                t_startPos = realnStartPos[i];
                t_endPos   = realnEndPos[i];

//                if (t_startPos<t_lastEndPos)
//                    t_startPos = t_lastEndPos+1;

                if (m_verboseLevel>=1)
                {
                    cout << m_bamGenomes[t_LeftRefID].RefName << "\t"
                         << t_startPos+1 << "\t"
                         << t_endPos+1+1 << endl;
                }

                numRealignRegion += windowRealignmentAmplicon(ProbAligner, m_bamReader, m_fasta, t_LeftRefID, t_startPos, t_LeftRefID, t_endPos+1, realignObjects, m_graphEdgePruneLevel, m_topK, HasHaplotypeLikeRank, HasUpdateGraph, HasProcLargeHomopolymer, HasKeepDuplicate, m_minLength, m_minMapQual, m_maxMismatches);
                numRegion ++;

                t_lastEndPos = t_endPos+1;
            }

            // reader rewind
            m_bamReader.Rewind();
            m_bamReader.SetRegion(t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos);

            BamAlignment al;
            while (m_bamReader.GetNextAlignment(al))
            {
                unordered_map<string,BamAlignment>::iterator iter = realignObjects.find(GenericBamAlignmentTools::getBamAlignmentID(al));
                if (iter==realignObjects.end())
                {
                    m_bamWriter.SaveAlignment(al);
                }
                else
                {
                    m_bamWriter.SaveAlignment(iter->second);
                }
            }

            clock_t tEnd = clock();

            if (m_verboseLevel>=1)
            {
                cout << "time elapsed " << double(tEnd-tBegin)/CLOCKS_PER_SEC/60.0 << " minutes"
                     << ", process " << numRealignRegion << " of " << numRegion << " regions"
                     << ", and "
                     << "re-align " << realignObjects.size() << " reads." << endl;
            }
        }
    }else
    {
        // iterate over genomes
        for (int i=0; i<m_bamGenomes.size(); ++i)
        {

            t_LeftRefID  = i;
            t_RightRefID = i;

            int t_lastEndPos = -1;
            // iterate over regions
            for (t_LeftRefPos=0; t_LeftRefPos<m_bamGenomes[i].RefLength; t_LeftRefPos+=t_RegionSize)
            {
                clock_t tBegin = clock();

                if (t_LeftRefPos+t_RegionSize<m_bamGenomes[i].RefLength)
                    t_RightRefPos = t_LeftRefPos+t_RegionSize;
                else
                    t_RightRefPos = m_bamGenomes[i].RefLength;

                // collect possible realignment sites in window
                vector<int> realnStartPos, realnEndPos;
                windowRealignSites(m_bamReader, m_fasta, t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos, realnStartPos, realnEndPos, m_flankSize, HasKeepDuplicate, m_repeatLength, m_minLength, m_minMapQual, m_maxMismatches, m_minSnpRead, m_minInDelRead, false);

                if (HasMNP)
                {
                    vector<int> mnpStartPos, mnpEndPos;
                    windowMnpSites(m_bamReader, t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos, mnpStartPos, mnpEndPos, m_flankSize, HasKeepDuplicate, m_minLength, m_minMapQual, m_maxMismatches, m_minSnpRead);
                    for (int i=0; i<mnpStartPos.size(); i++)
                    {
                        realnStartPos.push_back(mnpStartPos[i]);
                        realnEndPos.push_back(mnpEndPos[i]);
                    }
                }

                vector<int> mergeRealnStartPos, mergeRealnEndPos;
                windowRealignSitesMerge(realnStartPos, realnEndPos, mergeRealnStartPos, mergeRealnEndPos);
                realnStartPos = mergeRealnStartPos;
                realnEndPos   = mergeRealnEndPos;

                // realignment in each site
                for (int i=0; i<realnStartPos.size(); i++)
                {
                    int t_startPos, t_endPos;
                    t_startPos = realnStartPos[i];
                    t_endPos   = realnEndPos[i];

//                    if (t_startPos<t_lastEndPos)
//                        t_startPos = t_lastEndPos+1;

                    if (m_verboseLevel>=1)
                    {
                        cout << m_bamGenomes[t_LeftRefID].RefName << "\t"
                             << t_startPos+1 << "\t"
                             << t_endPos+1+1 << endl;
                    }

                    numRealignRegion += windowRealignmentAmplicon(ProbAligner, m_bamReader, m_fasta, t_LeftRefID, t_startPos, t_LeftRefID, t_endPos+1, realignObjects, m_graphEdgePruneLevel, m_topK, HasHaplotypeLikeRank, HasUpdateGraph, HasProcLargeHomopolymer, HasKeepDuplicate, m_minLength, m_minMapQual, m_maxMismatches);
                    numRegion ++;

                    t_lastEndPos = t_endPos+1;
                }

                // reader rewind
                m_bamReader.Rewind();
                m_bamReader.SetRegion(t_LeftRefID, t_LeftRefPos, t_RightRefID, t_RightRefPos);

                BamAlignment al;
                while (m_bamReader.GetNextAlignment(al))
                {
                    unordered_map<string,BamAlignment>::iterator iter = realignObjects.find(GenericBamAlignmentTools::getBamAlignmentID(al));
                    if (iter==realignObjects.end())
                    {
                        m_bamWriter.SaveAlignment(al);
                    }
                    else
                    {
                        m_bamWriter.SaveAlignment(iter->second);
                    }
                }

                clock_t tEnd = clock();

                if (m_verboseLevel>=1)
                {
                    cout << "time elapsed " << double(tEnd-tBegin)/CLOCKS_PER_SEC/60.0 << " minutes"
                         << ", process " << numRealignRegion << " of " << numRegion << " regions"
                         << ", and "
                         << "re-align " << realignObjects.size() << " reads." << endl;
                }
            }
        }
    }

    // processing information
    cout << "Totally re-align " << numRealignRegion << " of " << numRegion << " regions"
         << ", and "
         << "process " << realignObjects.size() << " reads." << endl;

    // close reader, wirter, and fasta
    m_bamReader.Close();
    m_fasta.Close();
    m_bamWriter.Close();

    return 1;
}



int windowRealignmentAmplicon(GenericProbabilisticAlignment &ProbAligner, BamReader &BamReader, Fasta &fasta,
                              int leftRefID, int leftRefPos, int rightRefID, int rightRefPos,
                              unordered_map<string, BamAlignment> &alignObjects,
                              int graphPruneLevel, int topK, bool haplotypeLikeRank, bool updateGraphByNewAlignment, bool procLargeHomopolymer,
                              bool keepDuplicate, int minReadLen, int minMapQual, double maxMismatchFrac)
{

    // genome sequence
    string t_fastaSeq;
    fasta.GetSequence(leftRefID, leftRefPos, rightRefPos-1, t_fastaSeq);

    // set of sequences in local region to be re-aligned
    vector<Realigner_Sequence_t> t_seqInLocalRegion;

    // set of unique sequences
    map<tuple<int,string,int,int>, list<int>> t_uniqueSeqSet;

    // rewind BAM reader
    BamReader.Rewind();

    // set BAM region
    BamReader.SetRegion(leftRefID, leftRefPos, rightRefID, rightRefPos);

    // retrieve alignment
    BamAlignment al;
    while(BamReader.GetNextAlignment(al))
    {
        BamAlignment alignObj;
        // skip if alignment is not in the region
        if (al.Position>=rightRefPos || al.GetEndPosition()<=leftRefPos)
            continue;

        // skip if read fail quality assessments
        if (readIsFiltered(al, keepDuplicate, minReadLen, minMapQual, maxMismatchFrac))
            continue;

        // check whether it is in realigned set
        string t_ID = GenericBamAlignmentTools::getBamAlignmentID(al);
        unordered_map<string,BamAlignment>::iterator iter = alignObjects.find(t_ID);

        if (iter!=alignObjects.end())
        {
            alignObj = BamAlignment(iter->second);
        }else
        {
            alignObj = BamAlignment(al);
        }

        // extract local sequences
        string t_localRead, t_localGenome;
        Cigar  t_cigar;
        BamMD  t_md;
        int    t_numMismatch, t_numInDel;
        GenericBamAlignmentTools::getLocalAlignment(alignObj, leftRefPos, rightRefPos-leftRefPos,
                                                    t_localRead, t_localGenome, t_cigar, t_md,
                                                    t_numMismatch, t_numInDel);

        if (t_localRead.empty() || t_localGenome.empty())
            continue;

        // put into set
        Realigner_Sequence_t t_seq;
        t_seq.m_ID           = GenericBamAlignmentTools::getBamAlignmentID(alignObj);
        t_seq.m_sequence     = t_localRead;
        t_seq.t_sequence     = t_localRead;
        t_seq.t_cigar        = t_cigar;
        t_seq.t_md           = t_md;
        t_seq.t_numMismatch  = t_numMismatch;
        t_seq.t_numInDel     = t_numInDel;
        t_seq.t_mapQualScore = alignObj.MapQuality;
        t_seq.t_realign      = false;

        if (alignObj.Position>leftRefPos)
            t_seq.m_startPositionShift = alignObj.Position-leftRefPos;
        else
            t_seq.m_startPositionShift = 0;

        if (alignObj.GetEndPosition()<rightRefPos)
            t_seq.m_endPositionShift = rightRefPos-alignObj.GetEndPosition();
        else
            t_seq.m_endPositionShift = 0;

        t_seq.t_type = READ_CROSS_REGION;
        if (t_seq.m_startPositionShift>0)
            t_seq.t_type = READ_BEGIN_REGION;
        if (t_seq.m_endPositionShift>0)
            t_seq.t_type = READ_END_REGION;

        t_seq.init();

        t_seqInLocalRegion.push_back(t_seq);
    }

    // resort sequences into unique sequence
    for (int i=0; i<t_seqInLocalRegion.size(); ++i)
    {
        Realigner_Sequence_t t_seq = t_seqInLocalRegion[i];

        tuple<int,string,int,int> t_key(t_seq.t_type, t_seq.t_sequence, t_seq.t_numMismatch, t_seq.t_numInDel);
        map<tuple<int,string,int,int>,list<int>>::iterator t_uniqueSeqSetIter = t_uniqueSeqSet.find(t_key);
        if (t_uniqueSeqSetIter==t_uniqueSeqSet.end())
        {
            t_uniqueSeqSet[t_key] = list<int>(1, i);
        }else
        {
            t_uniqueSeqSet[t_key].emplace_back(i);
        }
    }

    // debug
//    cout << t_uniqueSeqSet.size() << endl;

    // construct the graph to check the alignment consistency
    GenericDagGraph alignGraph;
    vector<string>  alignGraphReads;
    vector<Cigar>   alignGraphReadCigars;
    vector<int>     alignGraphReadStarts;

    // set of aligned reads to construct the graph
    for (int i=0; i<t_seqInLocalRegion.size(); ++i)
    {
        alignGraphReads.push_back(t_seqInLocalRegion[i].t_sequence);
        alignGraphReadCigars.push_back(t_seqInLocalRegion[i].t_cigar);
        alignGraphReadStarts.push_back(t_seqInLocalRegion[i].m_startPositionShift);
    }
    // construct the alignment graph
    alignGraph.buildDagGraph(t_fastaSeq, alignGraphReads, alignGraphReadCigars, alignGraphReadStarts);
    if (graphPruneLevel>0)
    {
        alignGraph.edgePruning2(graphPruneLevel-1, graphPruneLevel);
//        alignGraph.edgePruning(graphPruneLevel);
    }

    // debug
//    cout << alignGraph << endl;

    // re-align if the graph is inconsistent
    if (!alignGraph.isConsistentGraph())
    {
        GenericDagGraph haplotypeGraph;

        // pick top rank path(s) as alternative haplotype sequence(s)
        vector<string>       topRankAlignGraphPaths;
        vector<list<Vertex>> topRankAlignGraphPathVertexs;
        vector<double>       topRankAlignGraphPathWeights;
        alignGraph.topRankPathsExcludeGenome(10+topK, topRankAlignGraphPaths, topRankAlignGraphPathVertexs, topRankAlignGraphPathWeights);

        vector<int> haplotypeRank;
        if (haplotypeLikeRank)
        {
            rankHaplotypeByLikelihood(ProbAligner, t_fastaSeq, topRankAlignGraphPaths, t_seqInLocalRegion, haplotypeRank);
        }else
        {
            for (int i=0; i<topRankAlignGraphPaths.size(); i++)
                haplotypeRank.push_back(i);
        }

        vector<string>       haplotypeSequences;
        vector<Cigar>        haplotypeCigars;
        vector<int>          haplotypeStarts;
        for (int j=0; j<topK; ++j)
        {
            if (haplotypeRank.empty() || j>=haplotypeRank.size())
                continue;

            int i = haplotypeRank[j];
            haplotypeSequences.push_back(topRankAlignGraphPaths[i]);

            Cigar h_cigar;
            alignGraph.pathCigar(topRankAlignGraphPathVertexs[i], h_cigar);
            haplotypeCigars.push_back(h_cigar);

            haplotypeStarts.push_back(0);
        }


        // construct the haplotype graph
        haplotypeGraph.buildDagGraph(t_fastaSeq, haplotypeSequences, haplotypeCigars, haplotypeStarts);

        // debug
//        cout << haplotypeGraph << endl;

        // re-align unique sequences
        map<tuple<int,string,int,int>,list<int>>::iterator t_uniqueSeqIter;
        for (t_uniqueSeqIter=t_uniqueSeqSet.begin(); t_uniqueSeqIter!=t_uniqueSeqSet.end(); ++t_uniqueSeqIter)
        {
            if (get<2>(t_uniqueSeqIter->first)==0 && get<3>(t_uniqueSeqIter->first)==0)
                continue;

            string r_seq  = get<1>(t_uniqueSeqIter->first);
            int    r_type = get<0>(t_uniqueSeqIter->first);


            // realignment computation
            int     d_start;
            string  d_alnQuery;
            string  d_alnGenome;

            Cigar r_cigar;
            haplotypeGraph.alignReadToGraph(r_seq, ProbAligner, r_type, 1.0, procLargeHomopolymer, r_cigar, d_start, d_alnGenome, d_alnQuery);

            int r_numMismatch = GenericBamAlignmentTools::numMismatch(d_alnQuery, d_alnGenome);
            if (r_numMismatch<=get<2>(t_uniqueSeqIter->first))
            {
                list<int>::iterator t_readIter = t_uniqueSeqIter->second.begin();
                for (; t_readIter!=t_uniqueSeqIter->second.end(); ++t_readIter)
                {
                    int i = *t_readIter;

                    t_seqInLocalRegion[i].t_realign        = true;
                    t_seqInLocalRegion[i].t_cigar          = r_cigar;
                    t_seqInLocalRegion[i].t_alnLocalRead   = d_alnQuery;
                    t_seqInLocalRegion[i].t_alnLocalGenome = d_alnGenome;
                }
            }
        }

        // save alignments into BAM file
        map<string, int> t_seqInLocalRegionSet;
        for (int i=0; i<t_seqInLocalRegion.size(); ++i)
        {
            if (t_seqInLocalRegion[i].t_numMismatch==0 && t_seqInLocalRegion[i].t_numInDel==0)
                continue;

            t_seqInLocalRegionSet[t_seqInLocalRegion[i].m_ID] = i;
        }

        // rewind the reader
        BamReader.Rewind();

        // set the region
        BamReader.SetRegion(leftRefID, leftRefPos, rightRefID, rightRefPos);

        // retrieve alignment
        while(BamReader.GetNextAlignment(al))
        {
            // skip if alignment is not in the region
            if (al.Position>=rightRefPos || al.GetEndPosition()<=leftRefPos)
                continue;

            // skip if read fail quality assessment
            if (readIsFiltered(al, keepDuplicate, minReadLen, minMapQual, maxMismatchFrac))
                continue;

            map<string,int>::iterator iter = t_seqInLocalRegionSet.find(GenericBamAlignmentTools::getBamAlignmentID(al));
            // skip if read is not in the re-alignment set
            if (iter==t_seqInLocalRegionSet.end())
                continue;

            // skip if read is not re-aligned
            if (!t_seqInLocalRegion[iter->second].t_realign)
                continue;

            BamAlignment alignObj;
            // check whether it is realigned previously
            unordered_map<string,BamAlignment>::iterator iter2 = alignObjects.find(t_seqInLocalRegion[iter->second].m_ID);
            if (iter2!=alignObjects.end())
            {
                alignObj = BamAlignment(iter2->second);
            }else
            {
                alignObj = BamAlignment(al);
            }

            // change old alignment to new alignment
            BamAlignment n_alignObj;
            GenericBamAlignmentTools::changeLocalAlignment(alignObj, leftRefPos, rightRefPos-leftRefPos,
                                                           t_seqInLocalRegion[iter->second].t_alnLocalRead,
                                                           t_seqInLocalRegion[iter->second].t_alnLocalGenome,
                                                           n_alignObj);

            // update recording
            if (iter2==alignObjects.end())
            {
                alignObjects[t_seqInLocalRegion[iter->second].m_ID] = n_alignObj;
            }
            else
            {
                alignObjects[t_seqInLocalRegion[iter->second].m_ID] = n_alignObj;
            }
        }

    }else
    {
        return 0;
    }

    return 1;
}

struct SortSnpByPosition
{
    bool operator()(tuple<long,char>& a ,tuple<long,char>& b)
    {
        return (get<0>(a) < get<0>(b));
    }
};

int distanceToMnpBlock(int t0, int t1, int s)
{
    if (s<t0)
        return -2;

    if (s<t1)
        return -1;

    return (s-t1);
}

void windowMnpSites(BamReader &bamReader, int leftRefID, int leftRefPos, int rightRefID, int rightRefPos, vector<int> &mnpStartPos, vector<int> &mnpEndPos, int flankSize, bool keepDuplicate, int minReadLen, int minMapQual, double maxMismatchFrac, int minSnpRead)
{
    // rewind
    bamReader.Rewind();
    // set region
    bamReader.SetRegion(leftRefID, leftRefPos, rightRefID, rightRefPos);

    // search SNP in this region
    map<tuple<long,char>,int> SnpSiteCounts;

    // iterate over alignments in this region
    BamAlignment al;
    while (bamReader.GetNextAlignment(al))
    {
        if (readIsFiltered(al, keepDuplicate, minReadLen, minMapQual, maxMismatchFrac))
            continue;

        // get mismatches in the alignment
        vector<long> localSnpPos;
        vector<char> localSnpBase;
        GenericBamAlignmentTools::getBamAlignmentMismatches(al, localSnpPos, localSnpBase);

        // save
        for (int i=0; i<localSnpPos.size(); i++)
        {
            tuple<long,char> key(localSnpPos[i],localSnpBase[i]);
            if (SnpSiteCounts.find(key)==SnpSiteCounts.end())
            {
                SnpSiteCounts[key] = 1;
            }else
            {
                SnpSiteCounts[key] += 1;
            }
        }
    }

    // filter
    vector<tuple<long,char>> SnpSites;
    for (map<tuple<long,char>,int>::iterator iter=SnpSiteCounts.begin(); iter!=SnpSiteCounts.end(); iter++)
    {
        if (get<0>(iter->first)<leftRefPos || get<0>(iter->first)>=rightRefPos)
            continue;

        if (iter->second<=minSnpRead)
            continue;

        SnpSites.push_back(iter->first);
    }

    // sort SNP
    sort(SnpSites.begin(), SnpSites.end(), SortSnpByPosition());



    // cluster snp to mnp
    vector<tuple<int,int,map<char,int>>> mnp;
    for (int i=0; i<SnpSites.size(); i++)
    {

        if (i==0)
        {
            // allele count
            map<char,int> ac;
            ac[get<1>(SnpSites[i])] = 1;
            // save
            mnp.push_back(tuple<int,int,map<char,int>>(get<0>(SnpSites[i]),get<0>(SnpSites[i])+1,ac));
            continue;
        }

        int j = mnp.size()-1;
        int t0 = get<0>(mnp[j]);
        int t1 = get<1>(mnp[j]);
        int s = get<0>(SnpSites[i]);

        int tl = t1-t0;

        if (s<leftRefPos || s>=rightRefPos)
            continue;

        if (distanceToMnpBlock(t0,t1,s)<0)
            continue;

        if (s-t1<2 && tl<17)
        {
            get<1>(mnp[j]) = s+1;
            if (get<2>(mnp[j]).find(get<1>(SnpSites[i]))==get<2>(mnp[j]).end())
            {
                get<2>(mnp[j])[get<1>(SnpSites[i])] = 1;
            }else
            {
                get<2>(mnp[j])[get<1>(SnpSites[i])] += 1;
            }
        }else
        {
            // allele count
            map<char,int> ac;
            ac[get<1>(SnpSites[i])] = 1;
            // save
            mnp.push_back(tuple<int,int,map<char,int>>(get<0>(SnpSites[i]),get<0>(SnpSites[i])+1,ac));
        }
    }

    // save
    for (int i=0; i<mnp.size(); i++)
    {
        int t0 = get<0>(mnp[i]);
        int t1 = get<1>(mnp[i]);

        if (t1>t0+1)
        {
            mnpStartPos.push_back(t0-flankSize);
            mnpEndPos.push_back(t1+flankSize);
        }
    }
}
