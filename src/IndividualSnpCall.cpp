#include "GenericIndividualSnpCall.h"
#include "GenericVcfTools.h"
using namespace GenericSequenceTools;

#include "utils/bamtools_options.h"
#include "utils/bamtools_utilities.h"
using namespace BamTools;


IndividualSnpCallTool::IndividualSnpCallTool()
    : HasInput(false)
    , HasFasta(false)
    , HasRegion(false)
    , HasConfig(false)
    , HasOutput(false)
    , HasDownSample(false)
    , m_downsample(300)
    , HasMinReadLen(false)
    , m_minReadLen(100)
    , HasMinMapQual(false)
    , m_minMapQual(10)
    , HasMaxMismatchFrac(false)
    , m_maxMismatchFrac(0.03)
    , HasMinSnpRead(false)
    , m_minSnpRead(2)
    , HasMinSnpFrac(false)
    , m_minSnpFrac(0.1)
    , HasVariantQualityFilter(false)
    , m_variantQualityFilter(30)
    , HasPloidy(false)
    , m_ploidy(2)
    , HasPriorType(false)
    , PriorType("full")
    , HasPrior(false)
    , m_prior(0.001)
    , HasTopK(false)
    , m_topK(3)
    , HasFlankingSize(false)
    , m_flankingSize(10)
    , HasBand(false)
    , m_band(1.)
    , HasGraphPruneLevel(false)
    , m_graphPruneLevel(1)
    , HasSample(false)
    , HasSampleList(false)
    , HasVerbosity(false)
    , m_verbosity(0)
{
    // set program details
    Options::SetProgramInfo("GenericSequenceTools SnpCall",
                            "haplotype-based SNP caller",
                            "[-region <REGION>] [-config <CONFIG>] [-compress] -bam <INPUT> -ref <FASTA> -out <OUTPUT>");

    // set option group
    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input and Output Options");
    Options::AddValueOption("-bam", "FILE", "the input BAM file", "", HasInput, InputFile, IO_Opts);
    Options::AddValueOption("-sample", "STR", "the name of sequencing sample", "", HasSample, SampleList, IO_Opts);
    Options::AddValueOption("-ref", "FILE", "the genome file", "", HasFasta, FastaFile, IO_Opts);
    Options::AddValueOption("-out", "FILE", "the output VCF file", "", HasOutput, OutputFile, IO_Opts);
    Options::AddValueOption("-config", "FILE", "model configure file", "", HasConfig, ConfigFile, IO_Opts);
    Options::AddValueOption("-down-sample","INT","sample data down to the specified depth [300]","",HasDownSample,DownSample,IO_Opts);

    // set option group
    OptionGroup* Rg_Opts = Options::CreateOptionGroup("Genome Region Options");
    Options::AddValueOption("-region", "STR", "set region of interesting, e.g. chr:p0..chr:p1 for a segment, chr:0 for a chromosome", "", HasRegion, Region, Rg_Opts);

    // set option group
    OptionGroup* SnpCall_Opts = Options::CreateOptionGroup("SNP Calling Options");
    Options::AddValueOption("-ploidy", "INT", "the ploidy of the individual [2]", "", HasPloidy, Ploidy, SnpCall_Opts);
    Options::AddValueOption("-edge-prune-leve", "INT", "prune away edges in graph if edge count below the value [1]", "", HasGraphPruneLevel, GraphPruneLevel, SnpCall_Opts);
    Options::AddValueOption("-topK", "INT", "the number of consensus sequences [3]", "", HasTopK, TopK, SnpCall_Opts);
    Options::AddValueOption("-side-size", "INT", "the size of flanking segments of scanning window [10]", "", HasFlankingSize, FlankingSize, SnpCall_Opts);
    Options::AddValueOption("-band", "FLT", "the band of alignment computation [1.0]", "", HasBand, Band, SnpCall_Opts);
    Options::AddValueOption("-prior-type", "STR", "prior type: full or flat [full]", "", HasPriorType, PriorType, SnpCall_Opts);
    Options::AddValueOption("-prior-value", "FLT", "prior probability of variant [0.001]", "", HasPrior, Prior, SnpCall_Opts);
    Options::AddValueOption("-snp-quality", "FLT", "output SNP if variant quality is greater than the value [30]", "", HasVariantQualityFilter, VariantQualityFilter, SnpCall_Opts);

    // set option group
    OptionGroup* ReadFilter_Opts = Options::CreateOptionGroup("Read Filtration Options");
    Options::AddValueOption("-min-read-len", "INT", "the minimum read length [100]", "", HasMinReadLen, MinReadLen, ReadFilter_Opts);
    Options::AddValueOption("-min-map-qual", "INT", "the minimum map quality score [10]", "", HasMinMapQual, MinMapQual, ReadFilter_Opts);
    Options::AddValueOption("-max-mismatches", "FLT", "the maximum ratio of mismatches in read [0.03]", "", HasMaxMismatchFrac, MaxMismatchFrac, ReadFilter_Opts);

    // set option group
    OptionGroup* SnpFilter_Opts = Options::CreateOptionGroup("Simple SNP Filtration Options");
    Options::AddValueOption("-min-snp-read", "INT", "the minimum number of reads supporting SNP [2]", "", HasMinSnpRead, MinSnpRead, SnpFilter_Opts);
    Options::AddValueOption("-min-snp-frac", "FLT", "the minimum fraction of reads supporting SNP [0.1]", "", HasMinSnpFrac, MinSnpFrac, SnpFilter_Opts);

    // set option group
    OptionGroup* Verbose_Opts = Options::CreateOptionGroup("Verbosity Options");
    Options::AddValueOption("-verbose", "INT", "verbose level", "", HasVerbosity, Verbosity, Verbose_Opts);
}

int IndividualSnpCallTool::Help()
{
    Options::DisplayHelp();
    return 0;
}


int IndividualSnpCallTool::Run(int argc, char *argv[])
{
    // just GenericSequenceTools SnpCall
    if (argc==2)
        return Help();

    // parse the command line options
    Options::Parse(argc, argv, 1);

    // BAM file
    if (!HasInput){
        cerr << "GenericSequenceTools SnpCall ERROR: BAM file not provided... Aborting" << endl;
        return false;
    }

    // open input files
    if (!m_bamReader.Open(InputFile)){
        cerr << "GenericSequenceTools SnpCall ERROR: could not open input BAM file... Aborting" << endl;
        return false;
    }

    // if input is not stdin and a region is provided, look for index files
    if (HasInput && HasRegion){
        if (!m_bamReader.LocateIndex()){
            cerr << "GenericSequenceTools SnpCall ERROR: could not locate index file... Aborting" << endl;
            return false;
        }
    }

    // set region if specified
    if (HasRegion){
        if (BamTools::Utilities::ParseRegionString(Region, m_bamReader, m_bamRegion)){
            if (m_bamReader.HasIndex()){
                if (!m_bamReader.SetRegion(m_bamRegion)){
                    cerr << "GenericSequenceTools SnpCall ERROR: set region failed... Aborting" << endl;
                    m_bamReader.Close();
                    return false;
                }
            }

        }else{
            cerr << "GenericSequenceTools SnpCall ERROR: could not parse REGION: " << Region << endl;
            m_bamReader.Close();
            return false;
        }
    }

    // Fasta file
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

    // Probabilistic Alignment Machine
    // make sure output configure file is provided
    if (!HasConfig)
    {
        cerr << "GenericSequenceTools SnpCall ERROR: no output configure file provided... Aborting" << endl;
        return 0;
    }
    GenericProbabilisticAlignment ProbAligner;
    ifstream inputConfig;
    inputConfig.open(ConfigFile);
    inputConfig >> ProbAligner;
    inputConfig.close();


    // Snp Call Options
    if (HasPloidy)
        m_ploidy = stoi(Ploidy);

    if (HasTopK)
        m_topK = stoi(TopK);

    if (HasFlankingSize)
        m_flankingSize = stoi(FlankingSize);

    if (HasGraphPruneLevel)
        m_graphPruneLevel = stoi(GraphPruneLevel);

    if (HasBand)
        m_band = stod(Band);

    if (HasPrior)
        m_prior = stod(Prior);

    if (HasVariantQualityFilter)
        m_variantQualityFilter = stod(VariantQualityFilter);

    // Read filter options
    if (HasMinReadLen)
        m_minReadLen = stoi(MinReadLen);

    if (HasMinMapQual)
        m_minMapQual = stoi(MinMapQual);

    if (HasMaxMismatchFrac)
        m_maxMismatchFrac = stod(MaxMismatchFrac);

    // Snp filter options
    if (HasMinSnpRead)
        m_minSnpRead = stoi(MinSnpRead);

    if (HasMinSnpFrac)
        m_minSnpFrac = stod(MinSnpFrac);

    // verbosity
    if (HasVerbosity)
        m_verbosity = stoi(Verbosity);

    // executation
    vector<GenericVariant> SnpCallResults;

    VariantCallSetting SnpCallSetting;
    SnpCallSetting.m_ploidy               = m_ploidy;
    SnpCallSetting.m_priorType            = PriorType;
    SnpCallSetting.m_prior                = m_prior;
    SnpCallSetting.m_variantQualityFilter = m_variantQualityFilter;
    SnpCallSetting.m_topK                 = m_topK;
    SnpCallSetting.m_flankingSize         = m_flankingSize;
    SnpCallSetting.m_graphPruneLevel      = m_graphPruneLevel;
    SnpCallSetting.m_band                 = m_band;

    if (HasSample)
    {
        for (int i=0; i<SampleList.size(); i++)
        {
            SnpCallSetting.m_sampleList.insert(SampleList[i]);
        }
    }

    GenericIndividualSnpCall IndividualSnpCaller(m_downsample, m_minReadLen, m_minMapQual, m_maxMismatchFrac, m_minSnpRead, m_minSnpFrac, m_verbosity);
    IndividualSnpCaller.call(m_fasta, m_bamReader, m_bamRegion, ProbAligner, SnpCallSetting, SnpCallResults);

    // save SNP calling results to VCF output
    RefVector chromosomes = m_bamReader.GetReferenceData();
    GenericVcfTools::write(chromosomes, OutputFile, SnpCallSetting, SnpCallResults);

    return 1;
}
