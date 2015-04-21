#include "GenericGraphTools.h"
#include "GenericSequenceGlobal.h"
using namespace GenericSequenceTools;

#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>

#include <cstdlib>
#include <climits>
#include <sstream>
using namespace std;

#include <api/BamAux.h>
#include <api/BamReader.h>
#include <utils/bamtools_fasta.h>
using namespace BamTools;

inline void Verbose(string msg)
{
    cerr << "[" << CurrentTime() << " PyroTools-LocalStrainCall] " << msg << endl;
}

// check the existence of a file
// reference: http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
inline bool file_exist(const string& filename)
{
    struct stat buffer;
    return (stat (filename.c_str(), &buffer) == 0);
}


LocalStrainCallTool::LocalStrainCallTool()
    :HasGenomeFile(false)
    ,HasBamAlnFile(false)
    ,HasRegion(false)
    ,topk(30)
    ,topK(50)
    ,edgePruneLevel(2)
    ,minUniqHit(2)
    ,minReadLength(100)
    ,useDiffRead(false)
    ,vcfOut(false)
    ,skipIndel(true)
    ,dissim(0.017)
    ,lambda(1.0)
    ,useGenome(false)
    ,prefix("local_strains")
    ,verbose(0)
{

}

int LocalStrainCallTool::commandOptionParser(int argc, char *argv[]){

    auto isHexString = [](string s)
    {
        for (auto c : s)
        {
            if (c=='x') return true;
            if (c=='X') return true;
        }
        return false;
    };

    int c;

    while (1)
    {
        static struct option long_options[] =
        {
            {"top-k",      required_argument, 0, 'k'},
            {"edge-prune", required_argument, 0, 'e'},
            {"len",        required_argument, 0, 'l'},
            {"mq",         required_argument, 0, 'q'},
            {"ff",         required_argument, 0, 0},
            {"diff-read",  no_argument,       0, 0},
            {"vcf",        no_argument,       0, 0},
            {"uniq-hit",   required_argument, 0, 'n'},
            {"top-K",      required_argument, 0, 'K'},
            {"indel",      no_argument,       0, 'I'},
            {"dissim",     required_argument, 0, 'd'},
            {"lambda",     required_argument, 0, 0},
            {"use-genome", no_argument,       0, 'g'},
            {"prefix",     required_argument, 0, 'p'},
            {"verbose",    required_argument, 0, 'v'},
            {"help",       no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "k:e:l:q:n:K:Id:gp:v:h", long_options, &option_index);

        // detect the end of the options
        if (c==-1) break;

        switch(c)
        {
        case 0:

            switch(option_index)
            {
            case 4:
                if (isHexString(optarg)){
                    filterFlag |= stoul(optarg,nullptr,16);
                }else{
                    filterFlag |= stoi(optarg);
                }
                break;
            case 5:
                useDiffRead = true;
                break;
            case 6:
                vcfOut = true;
                break;
            case 11:
                lambda = stod(optarg);
                break;
            }

            break;

        case 'k':
            topk=stoi(optarg);
            break;

        case 'K':
            topK=stoi(optarg);
            break;

        case 'e':
            edgePruneLevel=stoi(optarg);
            break;


        case 'l':
            minReadLength=stoi(optarg);
            break;

        case 'q':
            minMapQuality=stoi(optarg);
            break;

        case 'n':
            minUniqHit=stoi(optarg);
            break;

        case 'I':
            skipIndel = false;
            break;

        case 'd':
            dissim = stod(optarg);
            break;

        case 'g':
            useGenome = true;
            break;

        case 'p':
            prefix = optarg;
            break;

        case 'v':
            verbose=stoi(optarg);
            break;

        case 'h':
            Help();
            exit(EXIT_SUCCESS);
            break;

        default:
            abort();
        }
    }

    // genome file
    if (optind<argc)
    {
        HasGenomeFile = true;
        GenomeFile=argv[optind++];

        // check the existence of the genome file
        if (!file_exist(GenomeFile))
        {
            cerr << "[PyroTools-LocalStrainCall] error: "
                 << GenomeFile << " not existed" << endl;
            exit(EXIT_FAILURE);
        }

        // check the existence of the genome index file
        if (!file_exist(GenomeFile+".fai"))
        {
            cerr << "[PyroTools-LocalStrainCall] error: "
                 << (GenomeFile+".fai") << " not existed" << endl;
            exit(EXIT_FAILURE);
        }
    }
    if (!HasGenomeFile){
        cerr << "[PyroTools-LocalStrainCall] error: "
             << "genome file is not specified" << endl;
        exit(EXIT_FAILURE);
    }

    // BAM file
    if (optind<argc)
    {
        HasBamAlnFile = true;
        BamAlnFile = argv[optind++];

        // check the existence of the bam file
        BamReader bamReader;
        if (!bamReader.Open(BamAlnFile))
        {
            cerr << "[PyroTools-LocalStrainCall] error: "
                 << BamAlnFile << " not existed or invalid" << endl;
            exit(EXIT_FAILURE);
        }
    }
    if (!HasBamAlnFile){
        cerr << "[PyroTools-LocalStrainCall] error: "
             << "bam file is not specified" << endl;
        exit(EXIT_FAILURE);
    }

    // Region
    if (optind<argc)
    {
        HasRegion = true;
        RegionOfInterest = argv[optind++];
    }
    if (!HasRegion)
    {
        cerr << "[PyroTools-LocalStrainCall] error: "
             << "region of interest is not specified" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

int LocalStrainCallTool::Help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    PyroTools StrainCall [OPTIONS] <GENOME> <BAM> <REGION>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    It is to reconstruct the strains in a local region, and detect the SNPs and/or Indels in the strains." << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "--options for graph searching" << endl;
    cerr << "    -k,--top-k         the number of reporting paths of GraphConsensus (default:30) [INT]" << endl;
    cerr << "    -K,--top-K         the number of enumerating paths internally used by GraphConsensus (default:50) [INT]" << endl;
    cerr << "--options for graph construction" << endl;
    cerr << "    -e,--edge-prune    graph edge pruning based on the read hit counts (default:2) [INT]" << endl;
    cerr << "    -n,--uniq-hit      minimal duplicates of the unique read (default:2) [INT]" << endl;
    cerr << "    -l,--len           minimal read length (default:100) [INT]" << endl;
    cerr << "    -q,--mq            minimal map quality (default:15) [INT]" << endl;
    cerr << "    --ff               filter reads with the specific SAM/BAM flag(s) [INT]" << endl;
    cerr << "    --diff-read        use reads harboring differences" << endl;
    cerr << "--options for vcf output" << endl;
    cerr << "    --vcf              toggle on vcf output" << endl;
    cerr << "    -I,--indel         inclusively print indels" << endl;
    cerr << "--options for strain clustering" << endl;
    cerr << "    -d,--dissim        specify the threshold of strain dissimilarity (default:0.017) [FLT]" << endl;
    cerr << "    --lambda           specify the lambda value (default:1.0) [FLT]" << endl;
    cerr << "--options for strain inference" << endl;
    cerr << "    -g,--use-genome    include the genome sequence" << endl;
    cerr << "--options for output" << endl;
    cerr << "    -p,--prefix        specify the prefix of output files (default:local_strains) [STR]" << endl;
    cerr << "    -v,--verbose       the verbosity level (default:0) [INT]" << endl;
    cerr << "    -h,--help          print help message" << endl;
    return 0;
}

int LocalStrainCallTool::Run(int argc, char *argv[]){
    // print help message if no argument provided
    if (argc==2){
        Help();
        exit(EXIT_SUCCESS);
    }

    // parse the command line
    commandOptionParser(argc-1, argv+1);

    // run the local strain inference routine
    if (LocalStrainInference()==EXIT_FAILURE)
        return EXIT_FAILURE;


    return 0;
}


int LocalStrainCallTool::LocalStrainInference(){
    stringstream ssCMD;
    string CMD;

    // get the working directory
    size_t CWD_SIZE = PATH_MAX;
    char CWD[CWD_SIZE];
    getcwd(CWD, CWD_SIZE);

    // make a temporary directory
    mode_t TWD_MODE = ACCESSPERMS;
    string TWD = string(CWD)+"/temporary";
    mkdir(TWD.c_str(), TWD_MODE);
    // change the working directory to the temporary directory
    chdir(TWD.c_str());

    // the path to pyrotools
    string PATH_TO_PYROTOOLS_DIR = GetExecPath();
    string PYROTOOLS = PATH_TO_PYROTOOLS_DIR + "/" + "PyroTools";

    // extract local reads in the temporary directory
    if (verbose>=1) Verbose("Extract out the local reads");
    // command line for CropBam
    ssCMD.clear();
    ssCMD.str("");
    ssCMD << PYROTOOLS << " " << "CropBam" << " ";
    ssCMD << "--format" << " " << "fasta" << " ";
    ssCMD << "--roi" << " " << RegionOfInterest << " ";
    ssCMD << "--mq" << " " << minMapQuality << " ";
    ssCMD << "--len" << " " << minReadLength << " ";
    ssCMD << "--slen" << " " << 30 << " ";
    ssCMD << "--iden" << " " << 0.1 << " ";
    ssCMD << "--ff" << " " << filterFlag << " ";
    ssCMD << "--uniq" << " ";
    ssCMD << "--freq" << " " << "temp_reads.freq" << " ";
    ssCMD << "--freq-thres" << " " << minUniqHit << " ";
    ssCMD << BamAlnFile << " " << ">" << " " << "temp_reads.fas";
    CMD = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    // execute CropBam command
    system(CMD.c_str());

    // invoke GraphConsensus to reconstruct the local strains and call variants
    if (verbose>=1) Verbose("Reconstruct the local strains");
    // command line for GraphConsensus
    ssCMD << PYROTOOLS << " " << "GraphConsensus" << " ";
    ssCMD << "--score" << " " << "assign" << " ";
    ssCMD << "--top-k" << " " << topk << " ";
    ssCMD << "--top-K" << " " << topK << " ";
    ssCMD << "--edge-prune" << " " << edgePruneLevel << " ";
    ssCMD << "--edge-frac" << " " << 0.005 << " ";
    if (useDiffRead) ssCMD << "--diff-read" << " ";
    ssCMD << "--uniq-hit" << " " << minUniqHit << " ";
    ssCMD << "--len" << " " << minReadLength << " ";
    ssCMD << "--mq" << " " << minMapQuality << " ";
    ssCMD << "--vcf" << " " << "temp_local_strains.vcf" << " ";
    if (!skipIndel) ssCMD << "--indel" << " ";
    ssCMD << "--out-format" << " " << "fasta" << " ";
    ssCMD << GenomeFile << " " << BamAlnFile << " " << RegionOfInterest << " ";
    ssCMD << ">" << " " << "temp_local_strains.fas" << " ";
    CMD = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    // execute GraphConsensus command
    system(CMD.c_str());

    // use genome in inference
    if (useGenome)
    {
        ssCMD << "samtools faidx " << GenomeFile << " " << RegionOfInterest << " > temp_reference.fas";
        CMD = ssCMD.str();
        ssCMD.clear();
        ssCMD.str("");
        system(CMD.c_str());
        // merge two files
        ssCMD << "cat temp_reference.fas >> temp_local_strains.fas";
        CMD = ssCMD.str();
        ssCMD.clear();
        ssCMD.str("");
        system(CMD.c_str());
    }

    // modify sequence name
    ssCMD << "sed -i -e 's/^>/>" << prefix << "_/' temp_local_strains.fas";
    CMD = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    system(CMD.c_str());

    // XGSUTILS tools
    string STRAINCLUSTER = PATH_TO_PYROTOOLS_DIR + "/../xgsutils/faxutils/haplo_cluster";
    string STRAINSEARCH  = PATH_TO_PYROTOOLS_DIR + "/../xgsutils/faxutils/haplo_search";
    string SSWALIGN      = PATH_TO_PYROTOOLS_DIR + "/../xgsutils/alnutils/ssw_aligner";

    // clustering
    if (verbose>=1) Verbose("Local strain clustering");
    // command line for clustering
    ssCMD << STRAINCLUSTER << " ";
    ssCMD << "-t" << " " << dissim << " ";
    ssCMD << "temp_local_strains.fas" << " ";
    ssCMD << ">" << " " << "temp_local_strains_centroid.fas";
    CMD = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    system(CMD.c_str());

    /////////////////////////////////////////////
    // changed by Feng Zeng @ Feb 8 2015
//    // bowtie2 mapping
//    if (verbose>=1) Verbose("Bowtie2 dictionary construction");
//    // make the sequence dictionary
//    ssCMD << "bowtie2-build" << " ";
//    ssCMD << "temp_local_strains_centroid.fas" << " ";
//    ssCMD << "temp_bowtie2_dict" << " ";
//    ssCMD << ">/dev/null 2>&1";
//    CMD = ssCMD.str();
//    ssCMD.clear();
//    ssCMD.str("");
//    system(CMD.c_str());

    // changed by Feng Zeng @ Feb 8 2015
//    // read mapping
//    if (verbose>=1) Verbose("Bowtie2 read mapping");
//    // command line for read mapping
//    ssCMD << "bowtie2" << " ";
//    ssCMD << "-f" << " ";
//    ssCMD << "--end-to-end" << " ";
//    ssCMD << "--ignore-quals" << " ";
//    //ssCMD << "--norc" << " ";
//    //ssCMD << "--mp" << " " << 3 << " ";
//    //ssCMD << "--rdg" << " " << 8 << "," << 5 << " ";
//    ssCMD << "--all" << " ";
//    ssCMD << "-x" << " " << "temp_bowtie2_dict" << " ";
//    ssCMD << "-U" << " " << "temp_reads.fas" << " ";
//    ssCMD << "-S" << " " << "temp_reads.sam" << " ";
//    ssCMD << ">/dev/null 2>&1";
//    CMD = ssCMD.str();
//    ssCMD.clear();
//    ssCMD.str("");
//    system(CMD.c_str());
    //////////////////////////////////////////

    // compute all alignments
    if (verbose>=1) Verbose("Compute all alignments");
    // command line for alignment computation
    ssCMD << SSWALIGN << " ";
    ssCMD << "-f sam" << " ";
    ssCMD << "temp_local_strains_centroid.fas" << " ";
    ssCMD << "temp_reads.fas" << " ";
    ssCMD << "2>/dev/null" << " ";
    ssCMD << ">temp_reads.sam";
    CMD = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    system(CMD.c_str());

    // BIC-based model selection
    if (verbose>=1) Verbose("BIC model selection");
    // command line for model selection
    ssCMD << STRAINSEARCH << " ";
    ssCMD << "-f" << " " << "temp_local_strains_centroid.fas" << " ";
    ssCMD << "-l" << " " << lambda << " ";
    ssCMD << "temp_reads.sam" << " ";
    ssCMD << "temp_reads.freq" <<" ";
    ssCMD << ">" << " " << "temp_local_strains_bic.fas";
    CMD = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    system(CMD.c_str());

    // faidx
    ssCMD << "samtools faidx" << " ";
    ssCMD << "temp_local_strains_bic.fas";
    CMD = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    system(CMD.c_str());

    ssCMD << "cut -f1 temp_local_strains_bic.fas.fai > temp_local_strains_bic.list";
    CMD = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    system(CMD.c_str());

    // vcftools filter
    ssCMD << "vcftools --vcf temp_local_strains.vcf --keep temp_local_strains_bic.list --recode --stdout";
    ssCMD << " > temp_local_strains_bic.vcf";
    CMD = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    system(CMD.c_str());

    // report
    ssCMD << "cp temp_local_strains_bic.fas ";
    ssCMD << CWD << "/" << prefix << ".fas";
    CMD = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    system(CMD.c_str());

    if (vcfOut){
        ssCMD << "cp temp_local_strains_bic.vcf ";
        ssCMD << CWD << "/" << prefix << ".vcf";
        CMD = ssCMD.str();
        ssCMD.clear();
        ssCMD.str("");
        system(CMD.c_str());
    }

    // clean everything
    chdir(CWD);
    ssCMD << "rm -rf" << " ";
    ssCMD << TWD;
    CMD = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    system(CMD.c_str());

    return 0;
}
