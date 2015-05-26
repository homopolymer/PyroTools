#include "GenericGraphTools.h"
#include "GenericSequenceGlobal.h"
using namespace GenericSequenceTools;

#include <tuple>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
using namespace std;

#include <api/BamAux.h>
#include <api/BamReader.h>
#include <utils/bamtools_fasta.h>
using namespace BamTools;

#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>


inline void SB_Verbose(string msg)
{
    cerr << "[" << CurrentTime() << " PyroTools-StrainBinning] " << msg << endl;
}

// check the existence of a file
// reference: http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
inline bool file_exist(const string& filename)
{
    struct stat buffer;
    return (stat (filename.c_str(), &buffer) == 0);
}


// constructor
StrainBinning::StrainBinning()
    :HasGenomeFile(false)
    ,HasShortFile(false)
    ,HasLongFile(false)
    ,CR_CORES(1)
    ,RS_width(250)
    ,RS_space(1000)
    ,RS_diff(5)
    ,LSC_topk(10)
    ,LSC_edgePruneLevel(100)
    ,LSC_minUniqHit(2)
    ,LSC_minReadLength(100)
    ,LSC_minMapQuality(15)
    ,LSC_filterFlag(4)
    ,LSC_dissim(0.017)
    ,SVC_minReadLength(100)
    ,SVC_minMapQuality(15)
    ,SVC_filterFlag(4)
    ,SVC_minIdentity(0.9)
    ,SVC_minBaseQuality(10)
    ,SVC_minSnpHit(2)
    ,SVC_minSnpRatio(0.05)
    ,GC_editdist(5)
    ,MS_verbose(0)
{

}


// help message
int StrainBinning::Help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    PyroTools StrainBinning [OPTIONS] <GENOME> <BAM> <LONG_SEQS>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    It is to aggregate long reads or assembled contigs/scaffolds into strain bins" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "--positional options" << endl;
    cerr << "    <GENOME>           specify the fasta file of reference genome" << endl;
    cerr << "    <BAM>              specify the BAM file of short reads against genome" << endl;
    cerr << "    <LONG_SEQS>        specify the fasta file of long reads or assembled contigs/scaffolds" << endl;
    cerr << "--options for computational resource" << endl;
    cerr << "    -c,--cores         specify the number of cores to be used (default:1) [INT]" << endl;
    cerr << "    --bt2              specify bowtie2 indexing file of LONG_SEQS" << endl;
    cerr << "--options for regional specification" << endl;
    cerr << "    -w,--width         the width of genomic regions for the reconstruction of local strain sequences (default:250) [INT]" << endl;
    cerr << "    -s,--space         the space between two local reconstructing localtions (default:1000) [INT]" << endl;
    cerr << "    -x,--diff-num      minimal number of differences in a location (default:5) [INT]" << endl;
    cerr << "    -r,--roi           specify the region of interest for StrainBinning [STR]" << endl;
    cerr << "--options for StrainCall" << endl;
    cerr << "    -k,--top-k         the number of reporting paths of GraphConsensus (default:10) [INT]" << endl;
    cerr << "    -e,--edge-prune    graph edge pruning based on the read hit counts (default:100) [INT]" << endl;
    cerr << "    -n,--uniq-hit      minimal duplicates of the unique read (default:2) [INT]" << endl;
    cerr << "    -l,--len           minimal read length (default:100) [INT]" << endl;
    cerr << "    -q,--mq            minimal map quality (default:15) [INT]" << endl;
    cerr << "    --ff               filter reads with the specific SAM/BAM flag(s) (default:4) [INT]" << endl;
    cerr << "    --diff-read        merely use reads harboring differences" << endl;
    cerr << "    -g,--use-genome    include the genome sequence" << endl;
    cerr << "    -d,--dissim        specify the threshold of strain dissimilarity (default:0.017) [FLT]" << endl;
    cerr << "--options for SimpleVarCall" << endl;
    cerr << "    -l,--len           the same as previously defined" << endl;
    cerr << "    -q,--mq            the same as previously defined" << endl;
    cerr << "    --ff               the same as previously defined" << endl;
    cerr << "    --iden             minimal alignment identity (default:0.9) [FLT]" << endl;
    cerr << "    --bq               filter read alleles with base quality less than the value (default:10) [INT]" << endl;
    cerr << "    --snp-hit          filter SNPs with the hits less than the value (default:2) [INT]" << endl;
    cerr << "    --snp-ratio        filter SNPs with the ratio less than the value (default:0.05) [FLT]" << endl;
    cerr << "--options for graph clustering" << endl;
    cerr << "    --edit-dist        minimal edit distance between local strains and long sequences (default:5) [INT]" << endl;
    cerr << "--options for message" << endl;
    cerr << "    -v,--verbose       the verbosity level (default:0) [INT]" << endl;
    cerr << "    -h,--help          print help message" << endl;
    return 0;
}

// parse command line options
int StrainBinning::commandOptionParser(int argc, char *argv[]){
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
            // options for computational resource
            {"cores",      required_argument, 0, 'c'}, // 0
            {"bt2",        required_argument, 0, 0},   // 1
            // options for regional specification
            {"width",      required_argument, 0, 'w'}, // 2
            {"space",      required_argument, 0, 's'}, // 3
            {"diff-num",   required_argument, 0, 'x'}, // 4
            {"roi",        required_argument, 0, 'r'}, // 5
            // options for StrainCall
            {"top-k",      required_argument, 0, 'k'}, // 6
            {"edge-prune", required_argument, 0, 'e'}, // 7
            {"len",        required_argument, 0, 'l'}, // 8
            {"mq",         required_argument, 0, 'q'}, // 9
            {"ff",         required_argument, 0, 0},   // 10
            {"diff-read",  no_argument,       0, 0},   // 11
            {"uniq-hit",   required_argument, 0, 'n'}, // 12
            {"dissim",     required_argument, 0, 'd'}, // 13
            {"use-genome", no_argument,       0, 'g'}, // 14
            // options for SimpleVarCall
            {"iden",       required_argument, 0, 0},   // 15
            {"bq",         required_argument, 0, 0},   // 16
            {"snp-hit",    required_argument, 0, 0},   // 17
            {"snp-ratio",  required_argument, 0, 0},   // 18
            // options for graph clustering
            {"edit-dist",  required_argument, 0, 0},   // 19
            // options for message
            {"verbose",    required_argument, 0, 'v'}, // 20
            {"help",       no_argument,       0, 'h'}, // 21
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "c:w:s:x:r:k:e:l:q:n:d:gv:h", long_options, &option_index);

        // detect the end of the options
        if (c==-1) break;

        //
        switch(c)
        {
        case 0:

            switch(option_index)
            {
            case 1:
                CR_BT2LONG = optarg;
                break;
            case 10:
                if (isHexString(optarg)){
                    LSC_filterFlag |= stoul(optarg,nullptr,16);
                    SVC_filterFlag |= stoul(optarg,nullptr,16);
                }else{
                    LSC_filterFlag |= stoi(optarg);
                    SVC_filterFlag |= stoi(optarg);
                }
                break;
            case 11:
                LSC_useDiffRead = true;
                break;
            case 15:
                SVC_minIdentity = stod(optarg);
                break;
            case 16:
                SVC_minBaseQuality = stoi(optarg);
                break;
            case 17:
                SVC_minSnpHit = stoi(optarg);
                break;
            case 18:
                SVC_minSnpRatio = stod(optarg);
                break;
            case 19:
                GC_editdist = stoi(optarg);
                break;
            }

            break;

        case 'c':
            CR_CORES = stoi(optarg);
            break;

        case 'w':
            RS_width = stoi(optarg);
            break;

        case 's':
            RS_space = stoi(optarg);
            break;

        case 'x':
            RS_diff = stoi(optarg);
            break;

        case 'r':
            RS_roi = optarg;
            break;

        case 'k':
            LSC_topk = stoi(optarg);
            break;

        case 'e':
            LSC_edgePruneLevel = stoi(optarg);
            break;

        case 'l':
            LSC_minReadLength = stoi(optarg);
            SVC_minReadLength = stoi(optarg);
            break;

        case 'q':
            LSC_minMapQuality = stoi(optarg);
            SVC_minMapQuality = stoi(optarg);
            break;

        case 'n':
            LSC_minUniqHit = stoi(optarg);
            break;

        case 'd':
            LSC_dissim = stod(optarg);
            break;

        case 'g':
            LSC_useGenome = true;
            break;

        case 'v':
            MS_verbose = stoi(optarg);
            break;

        case 'h':
            Help();
            exit(EXIT_SUCCESS);
            break;

        default:
            abort();
        }
    }

    // genome
    // genome file
    if (optind<argc)
    {
        HasGenomeFile = true;
        GenomeFile=argv[optind++];

        // check the existence of the genome file
        if (!file_exist(GenomeFile))
        {
            cerr << "[PyroTools-StrainBinning] error: "
                 << GenomeFile << " not existed" << endl;
            exit(EXIT_FAILURE);
        }

        // check the existence of the genome index file
        if (!file_exist(GenomeFile+".fai"))
        {
            cerr << "[PyroTools-StrainBinning] error: "
                 << (GenomeFile+".fai") << " not existed" << endl;
            exit(EXIT_FAILURE);
        }
    }
    if (!HasGenomeFile){
        cerr << "[PyroTools-StrainBinning] error: "
             << "genome file is not specified" << endl;
        exit(EXIT_FAILURE);
    }

    // short -eq file
    if (optind<argc)
    {
        HasShortFile = true;
        ShortSeqFile = argv[optind++];

        // check the existence of the short-seq file
        BamReader bamReader;
        if (!bamReader.Open(ShortSeqFile))
        {
            cerr << "[PyroTools-StrainBinning] error: "
                 << ShortSeqFile << " not existed or invalid" << endl;
            exit(EXIT_FAILURE);
        }else{
            GenomeDict = bamReader.GetReferenceData();
        }
    }
    if (!HasShortFile){
        cerr << "[PyroTools-StrainBinning] error: "
             << "short-seq file is not specified" << endl;
        exit(EXIT_FAILURE);
    }

    // long seq file
    if (optind<argc)
    {
        HasLongFile = true;
        LongSeqFile = argv[optind++];

        // check the existence of the long-seq file
        if (!file_exist(LongSeqFile))
        {
            cerr << "[PyroTools-StrainBinning] error: "
                 << LongSeqFile << " not existed or invalid" << endl;
            exit(EXIT_FAILURE);
        }
    }
    if (!HasLongFile){
        cerr << "[PyroTools-StrainBinning] error: "
             << "long-seq file is not specified" << endl;
        exit(EXIT_FAILURE);
    }
}

// Running interface
int StrainBinning::Run(int argc, char *argv[]){
    // print help message if no argument provide
    if (argc==2){
        Help();
        exit(EXIT_SUCCESS);
    }

    // parse the commandline option
    commandOptionParser(argc-1,argv+1);

    // run the internal executation program
    Binning();

    return 0;
}

int StrainBinning::Binning()
{
    string CMD, Result, res;
    stringstream ssCMD, ssResult;

    ////////////////////////////////////////////////////////////////////
    // get the characteristics of long read
    if (MS_verbose>=1) SB_Verbose("collect long read metrics");

    int    LR_number = 0;
    double LR_mean_length = 0;

    if (!file_exist(LongSeqFile+".fai")){
        ssCMD << "samtools faidx " << LongSeqFile;
        CMD = ssCMD.str();
        ssCMD.clear();
        ssCMD.str("");
        system(CMD.c_str());
    }

    ifstream LongSeqFaiInput(LongSeqFile+".fai");
    while(getline(LongSeqFaiInput, res)){
        ssResult << res;
        int l;
        string dummy;
        ssResult >> dummy >> l;
        LR_number ++;
        LR_mean_length += l;
        ssResult.clear();
        ssResult.str("");
    }
    LongSeqFaiInput.close();

    LR_mean_length /= LR_number;

    /////////////////////////////////////////////////////////////////
    // get the working directory
    size_t CWD_SIZE = PATH_MAX;
    char CWD[CWD_SIZE];
    getcwd(CWD, CWD_SIZE);

    // make the temporary directory
    srand(time(0));
    mode_t TWD_MODE = ACCESSPERMS;
    string TWD = string(CWD)+"/binning_temporary_"+to_string(rand());
    mkdir(TWD.c_str(), TWD_MODE);
    // change the working directory to the temporary directory
    chdir(TWD.c_str());

    // the path to pyrotools
    string PATH_TO_PYROTOOLS_DIR = GetExecPath();
    string PYROTOOLS = PATH_TO_PYROTOOLS_DIR + "/" + "PyroTools";

    // the path to xgsutils
    string PATH_TO_XGSUTILS_DIR = PATH_TO_PYROTOOLS_DIR + "/../xgsutils";

    //////////////////////////////////////////////////////////////////
    map<int,bool> GenomeEmpty;
    map<int,string> GenomeTempDir;

    // check whether there is mapping on i-th genome
    for (int i=0; i<GenomeDict.size(); i++){
        ssCMD << "samtools view " << ShortSeqFile << " " << GenomeDict[i].RefName << " | head -n1";
        CMD = ssCMD.str();
        ssCMD.clear();
        ssCMD.str("");
        RunCmdGetReturn(CMD, Result);
        if (Result.empty()){
            GenomeEmpty[i] = true;
        }else{
            GenomeEmpty[i] = false;
        }
    }

    // creat a temporary directory for i-th genome
    for (int i=0; i<GenomeDict.size(); i++){
        string g_twd = TWD + "/" + to_string(i) + "_" + GenomeDict[i].RefName;
        mkdir(g_twd.c_str(), TWD_MODE);
        GenomeTempDir[i] = g_twd;
    }

    //////////////////////////////////////////////////////////////////
    // 1. find locations for local reconstruction
    if (MS_verbose>=1) SB_Verbose("find locations for local reconstruction");
    // loop over genomes
    for (int i=0; i<GenomeDict.size(); i++){
        if (MS_verbose>=1) SB_Verbose("invoke SimpleVarCall on genome " + GenomeDict[i].RefName);

        // if there is no mapping on this genome
        if (GenomeEmpty[i]){
            if (MS_verbose>=1) SB_Verbose("there is no mapping");
            continue;
        }

        // move into genome-specific temporary directory
        chdir(GenomeTempDir[i].c_str());

        // 1.1 call SNPs
        // simply call variant
        ssCMD << PYROTOOLS << " SimpleVarCall ";
        ssCMD << "--roi " << GenomeDict[i].RefName << " ";
        ssCMD << "--mq " << SVC_minMapQuality << " ";
        ssCMD << "--len " << SVC_minReadLength << " ";
        ssCMD << "--iden " << SVC_minIdentity << " ";
        ssCMD << "--ff " << SVC_filterFlag << " ";
        ssCMD << "--bq " << SVC_minBaseQuality << " ";
        ssCMD << "--snp-hit " << SVC_minSnpHit << " ";
        ssCMD << "--disable-indel" << " ";
        ssCMD << "--num-cores " << CR_CORES << " ";
        ssCMD << "-v " << MS_verbose << " ";
        ssCMD << GenomeFile << " " << ShortSeqFile << " ";
        ssCMD << "2>/dev/null 1>" << i << "_" << GenomeDict[i].RefName << "_snp.bed;";
        ssCMD << "sort -nk2,2 " << i << "_" << GenomeDict[i].RefName << "_snp.bed > "
              << i << "_" << GenomeDict[i].RefName << "_snp.sort.bed;";
        ssCMD << "awk -v r=" << SVC_minSnpRatio << " "
              << "'{if($6>=$7*r){print}}'" << " "
              << i << "_" << GenomeDict[i].RefName << "_snp.sort.bed" << " > "
              << i << "_" << GenomeDict[i].RefName << "_snp.sort.flt.bed;";
        ssCMD << "rm " << i << "_" << GenomeDict[i].RefName << "_snp.bed;";
        ssCMD << "rm " << i << "_" << GenomeDict[i].RefName << "_snp.sort.bed;";
        ssCMD << "awk '{print $1\"\\t\"$2\"\\t\"$3}' " << i << "_" << GenomeDict[i].RefName << "_snp.sort.flt.bed "
              << "> "   << i << "_" << GenomeDict[i].RefName << "_snp.bed;";
        ssCMD << "rm " << i << "_" << GenomeDict[i].RefName << "_snp.sort.flt.bed";
        CMD = ssCMD.str();
        ssCMD.clear();
        ssCMD.str("");
        system(CMD.c_str());

        // define bed by depth
        ssCMD << PATH_TO_XGSUTILS_DIR << "/bamutils/xBamToBedByDepth" << " "
              << ShortSeqFile << " "
              << "-g 10" << " "
              << "-r " << GenomeDict[i].RefName << " "
              << "2>/dev/null" << " "
              << "| awk '{$2=$2+100;$3=$3-100;if($3-$2>300){print $1\"\t\"$2\"\t\"$3}}'" << " "
              << "> " << i << "_" << GenomeDict[i].RefName << "_depth.bed";
        CMD = ssCMD.str();
        ssCMD.str("");
        system(CMD.c_str());

        // 1.2 define locations to perform StrainCall
        // make windows on the genome
        ssCMD << "bedtools makewindows "
              //<< "-g " << GenomeFile << ".fai "
              << "-b " << i << "_" << GenomeDict[i].RefName << "_depth.bed" << " "
              << "-w " << RS_width << " "
              << "-s " << RS_space << " "
              << "| awk '{if ($2>0){print}}'" << " "
              << "> " << i << "_" << GenomeDict[i].RefName << "_windows.bed";
        CMD = ssCMD.str();
        ssCMD.clear();
        ssCMD.str("");
        system(CMD.c_str());

        // select windows dependent on the number of differences
        ssCMD << "bedtools intersect "
              << "-wa -c "
              << "-a " << i << "_" << GenomeDict[i].RefName << "_windows.bed "
              << "-b " << i << "_" << GenomeDict[i].RefName << "_snp.bed "
              << "| awk -v n=" << RS_diff << " '{if($4>=n){print}}' "
              << "> " << i << "_" << GenomeDict[i].RefName << "_windows_lsc.bed";
        CMD = ssCMD.str();
        ssCMD.clear();
        ssCMD.str("");
        system(CMD.c_str());

        // load in windows
        vector<tuple<string,int,int>> LSC_windows;

        ifstream LscWindowsInput(to_string(i)+"_"+GenomeDict[i].RefName+"_windows_lsc.bed");
        while(getline(LscWindowsInput,res)){
            ssResult << res;
            string g_id;
            int g_p0,g_p1;
            ssResult >> g_id >> g_p0 >> g_p1;
            ssResult.clear();
            ssResult.str("");

            LSC_windows.emplace_back(make_tuple(g_id,g_p0,g_p1));
        }

        // invoke StrainCall in the window
        map<int,int> windows_strain_count;
        map<int,int> windows_strain_freq;
        map<int,string> windows_prefix1;
        map<int,string> windows_prefix2;
        for (int j=0; j<LSC_windows.size(); j++){

            string g_id;
            int g_p0,g_p1;
            g_id = get<0>(LSC_windows[j]);
            g_p0 = get<1>(LSC_windows[j]);
            g_p1 = get<2>(LSC_windows[j]);

            // 0
            string prefix1 = g_id + "_" + to_string(g_p0) + "_" + to_string(g_p1);
            windows_prefix1[j] = prefix1;

            if (MS_verbose>=1) SB_Verbose("processing "+g_id+":"+to_string(g_p0)+"-"+to_string(g_p1));

            ssCMD << PYROTOOLS << " StrainCall "
                  << "-k " << LSC_topk << " -K " << (LSC_topk+30>100?LSC_topk+30:100) << " "
                  << "-e " << LSC_edgePruneLevel << " "
                  << "-n " << LSC_minUniqHit << " "
                  << "-l " << LSC_minReadLength << " "
                  << "-q " << LSC_minMapQuality << " "
                  << "--ff " << LSC_filterFlag << " "
                  << (LSC_useDiffRead?"--diff-read":"") << " "
                  << (LSC_useGenome?"-g":"") << " "
                  << "-d " << LSC_dissim << " "
                  << "-p " << prefix1 << " "
//                  << "-p " << g_id << "_" << g_p0 << "_" << g_p1 << " "
                  << "-v " << (MS_verbose>1?1:0) << " "
                  << GenomeFile << " "
                  << ShortSeqFile << " "
                  << g_id << ":" << g_p0 << "-" << g_p1;
            CMD = ssCMD.str();
            ssCMD.clear();
            ssCMD.str("");
            system(CMD.c_str());

            // 1 : shift
            int shift = 100;
            int t_ndl,t_ndr;

            ssCMD << "echo -e \"" << g_id << "\\t" << (g_p0-100) << "\\t" << (g_p1-100) << "\"" << " "
                  << "| bedtools intersect -wa -c -a stdin -b " << i << "_" << GenomeDict[i].RefName << "_snp.bed" << " "
                  << "| awk '{print $4}'";
            CMD = ssCMD.str();
            ssCMD.clear();
            ssCMD.str("");
            RunCmdGetReturn(CMD, Result);
            t_ndl = stoi(Result);

            ssCMD << "echo -e \"" << g_id << "\\t" << (g_p0+100) << "\\t" << (g_p1+100) << "\"" << " "
                  << "| bedtools intersect -wa -c -a stdin -b " << i << "_" << GenomeDict[i].RefName << "_snp.bed" << " "
                  << "| awk '{print $4}'";
            CMD = ssCMD.str();
            ssCMD.clear();
            ssCMD.str("");
            RunCmdGetReturn(CMD, Result);
            t_ndr = stoi(Result);

            if (t_ndl>t_ndr) shift = -100;

            string prefix2 = g_id + "_" + to_string(g_p0+shift) + "_" + to_string(g_p1+shift);
            windows_prefix2[j] = prefix2;

            if (MS_verbose>=1) SB_Verbose("processing "+g_id+":"+to_string(g_p0+shift)+"-"+to_string(g_p1+shift));

            ssCMD << PYROTOOLS << " StrainCall "
                  << "-k " << LSC_topk << " -K " << (LSC_topk+30) << " "
                  << "-e " << LSC_edgePruneLevel << " "
                  << "-n " << LSC_minUniqHit << " "
                  << "-l " << LSC_minReadLength << " "
                  << "-q " << LSC_minMapQuality << " "
                  << "--ff " << LSC_filterFlag << " "
                  << (LSC_useDiffRead?"--diff-read":"") << " "
                  << (LSC_useGenome?"-g":"") << " "
                  << "-d " << LSC_dissim << " "
                  << "-p " << prefix2 << " "
//                  << "-p " << g_id << "_" << (g_p0+shift) << "_" << (g_p1+shift) << " "
                  << "-v " << (MS_verbose>1?1:0) << " "
                  << GenomeFile << " "
                  << ShortSeqFile << " "
                  << g_id << ":" << (g_p0+shift) << "-" << (g_p1+shift);
            CMD = ssCMD.str();
            ssCMD.clear();
            ssCMD.str("");
            system(CMD.c_str());

            // 2 : bowtie mapping

            if (MS_verbose>=1) SB_Verbose("map local strain sequences to long reads");

            ssCMD << "bowtie2 -a "
                  << "-x " << CR_BT2LONG << " "
                  << "-f -U " << prefix1 << ".fas" << " "
                  << "-U " << prefix2 << ".fas" << " "
                  << "-S " << j << "_" << prefix1 << ".sam" << " "
                  << "2>/dev/null";
            CMD = ssCMD.str();
            ssCMD.clear();
            ssCMD.str("");
            system(CMD.c_str());

            // 2.1 : filter reverse mapping
            ssCMD << "samtools view -h -F 0x10" << " "
                  << j << "_" << prefix1 << ".sam" << " "
                  << "> " << j << "_" << prefix1 << ".sam2" << ";"
                  << "mv " << j << "_" << prefix1 << ".sam2" << " "
                  << j << "_" << prefix1 << ".sam";
            CMD = ssCMD.str();
            ssCMD.clear();
            ssCMD.str("");
            system(CMD.c_str());

            // 3 : convert sam to table
            ssCMD << PATH_TO_XGSUTILS_DIR << "/gsautils/long_read_cns_table "
                  << "-d " << GC_editdist << " "
                  << "-t 0 "
                  << j << "_" << prefix1 << ".sam" << " "
                  << "> " << j << "_" << prefix1 << ".table";
            CMD = ssCMD.str();
            ssCMD.clear();
            ssCMD.str("");
            system(CMD.c_str());

            // 4 : graph laplacian spectrum and the number of putative strain
            if (MS_verbose>=1) SB_Verbose("bi-clustering local strain sequences and long reads");

            ssCMD << PATH_TO_XGSUTILS_DIR << "/gsautils/BipartiteGraphClustering.py "
                  << "--tbl=" << j << "_" << prefix1 << ".table" << " "
                  << "--p1=" << prefix1 << " "
                  << "--p2=" << prefix2 << " "
                  << "--t=" << 1e-7;
            CMD = ssCMD.str();
            ssCMD.clear();
            ssCMD.str("");
            RunCmdGetReturn(CMD, Result);

            // save the number
            int strain_count = stoi(Result);
            windows_strain_count[j] = strain_count;
            if (windows_strain_freq.count(strain_count)==0)
                windows_strain_freq[strain_count] = 1;
            else
                windows_strain_freq[strain_count] += 1;

            if (MS_verbose>=1) SB_Verbose("number of putative strains: "+to_string(strain_count));

        }

        // get the estimate of the number of putative strains
        // http://stackoverflow.com/questions/9370945/c-help-finding-the-max-value-in-a-map
        auto ptr = max_element(windows_strain_freq.begin(), windows_strain_freq.end(),
                               [](const pair<int,int>& a, const pair<int,int>& b){
                                  return a.second < b.second;});
        int estimate_strain_num = ptr->first;

        if (MS_verbose>=1) SB_Verbose("number of strains: " + to_string(estimate_strain_num));

        // re-clustering based on the estimated number of putative strains
        if (MS_verbose>=1) SB_Verbose("re-clustering long sequences based on the estimated number of putative strains");

        for (int j=0; j<LSC_windows.size(); j++){
            string prefix1 = windows_prefix1[j];
            string prefix2 = windows_prefix2[j];
            ssCMD << PATH_TO_XGSUTILS_DIR << "/gsautils/BipartiteGraphClustering.py "
                  << "--tbl=" << j << "_" << prefix1 << ".table" << " "
                  << "--p1=" << prefix1 << " "
                  << "--p2=" << prefix2 << " "
                  << "--t=" << 1e-7 << " "
                  << "--k=" << estimate_strain_num << " "
                  << "> " << j << "_" << prefix1 << ".long_read_pairs";
            CMD = ssCMD.str();
            ssCMD.clear();
            ssCMD.str("");
            system(CMD.c_str());
        }

        // move to the top temporary directory
        chdir(TWD.c_str());
    }

    ///////////////////////////////////////////////////////////////////
    // clean everything
//    chdir(CWD);
//    ssCMD << "rm -rf" << " ";
//    ssCMD << TWD;
//    CMD = ssCMD.str();
//    ssCMD.clear();
//    ssCMD.str("");
//    system(CMD.c_str());

    return 0;
}
