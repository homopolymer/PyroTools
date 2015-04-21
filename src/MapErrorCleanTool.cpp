#include "GenericBamAlignmentTools.h"
#include "GenericSequenceGlobal.h"
#include "GenericRegionTools.h"
using namespace GenericSequenceTools;

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamMultiReader.h"
#include "api/BamAlignment.h"
#include "api/BamAux.h"
using namespace BamTools;

#include <stdio.h>
#include <climits>
#include <random>
#include <tuple>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <unordered_set>
using namespace std;

#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <libgen.h>
// using openmp
#include <omp.h>

inline void Verbose(string msg)
{
    cerr << "[PyroTools-MapErrorClean] "
         << msg << endl;
}

MapErrorCleanTool::MapErrorCleanTool()
    :compressMode(1)
    ,step(1000)
    ,alnFlagMarker(0)
    ,mapQualThres(5)
    ,readLenThres(100)
    ,numThreads(1)
    ,verbose(0)
{
}

int MapErrorCleanTool::Help()
{
    cerr << "SYNOPSIS" << endl;
    cerr << "    PyroTools MapErrorClean [OPTIONS] <GENOME> <BAM>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    Mask the erroneously mapped reads" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "    -b              output the bam file" << endl;
    cerr << "    --bt2dict       specify the path to the bowtie2 indexed dictionary [STR]" << endl;
    cerr << "    --roi           specify the region of interest [STR]" << endl;
    cerr << "    --roi-list      specify the regions of interest" << endl;
    cerr << "    --mq            skip reads with mapping quality less than the value (default:5) [INT]" << endl;
    cerr << "    --ff            skip reads with the specified flags [INT]" << endl;
    cerr << "    --num-cores     specify the number of threads (default:1) [INT]" << endl;
    cerr << "    --window        specify the size of scanning window (default:1e+3) [INT]" << endl;
    cerr << "    -v,--verbose    specify the level of verbosity [INT]" << endl;
    cerr << "                    valid arguments:" << endl;
    cerr << "                        0 for quiet (default)" << endl;
    cerr << "                        1 for debug" << endl;
    cerr << "    -h,--help       print help message" << endl;
    return 0;
}

int MapErrorCleanTool::parseCommandLind(int argc, char *argv[])
{
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
            {"roi",           required_argument, 0, 0},
            {"roi-list",      required_argument, 0, 0},
            {"mq",            required_argument, 0, 0},
            {"len",           required_argument, 0, 0},
            {"ff",            required_argument, 0, 0},
            {"bt2dict",       required_argument, 0, 0},

            {"num-cores",     required_argument, 0, 0},
            {"window",        required_argument, 0, 0},
            {"verbose",       required_argument, 0, 'v'},
            {"help",          no_argument,       0, 'h'},
            {0,0,0,0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "bv:h", long_options, &option_index);

        // detect the end of the options
        if (c==-1) break;

        switch(c)
        {
        case 0:
            switch(option_index)
            {
            case 0:
                regionStrings.emplace_back(optarg);
                break;
            case 1:
                {
                    ifstream infile(optarg);
                    string line;
                    while (getline(infile,line))
                    {
                        if (!line.empty())
                            regionStrings.emplace_back(optarg);
                    }
                    infile.close();
                }
                break;
            case 2:
                mapQualThres=stoi(optarg);
                break;
            case 3:
                readLenThres=stoi(optarg);
                break;
            case 4:
                if (isHexString(string(optarg))){
                    alnFlagMarker |= stoul(optarg,nullptr,16);
                }else{
                    alnFlagMarker |= stoi(optarg);
                }
                break;
            case 5:
                bt2dict=optarg;
                break;
            case 6:
                numThreads=stoi(optarg);
                break;
            case 7:
                step=stod(optarg);
                break;
            default:
                abort();
            }
            break;
        case 'b':
            compressMode=0;
            break;
        case 'v':
            verbose=stoi(optarg);
            break;
        case 'h':
            Help();
            exit(EXIT_SUCCESS);
            break;
        case '?':
            exit(EXIT_FAILURE);
            break;
        default:
            abort();
        }
    }

    // genome file
    if (optind<argc)
    {
        genomeFile=argv[optind++];

        // check the existence of the genome file
        if (!FileExist(genomeFile))
        {
            cerr << "[PyroTools-SimpleVarCall] error: "
                   << genomeFile << " not existed" << endl;
            exit(EXIT_FAILURE);
        }

        // check the existence of the genome index file
        if (!FileExist(genomeFile+".fai"))
        {
            cerr << "[PyroTools-SimpleVarCall] error: "
                   << (genomeFile+".fai") << " not existed" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // bam file
    for (; optind<argc;)
    {
        bamFiles.emplace_back(argv[optind++]);

        // check the existence of the bam file
        auto f=*bamFiles.rbegin();
        BamReader bamReader;
        if (!bamReader.Open(f))
        {
            cerr << "[PyroTools-ConsensusGraph] error: "
                 << f << " not existed or invalid" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // canonicalize the regions
    if (regionStrings.empty())
    {
        BamMultiReader bamReader;
        bamReader.Open(bamFiles);
        RefVector genomeDict = bamReader.GetReferenceData();

        for (auto g : genomeDict)
        {
            regionStrings.emplace_back(g.RefName);
        }
    }

    return 0;
}

int MapErrorCleanTool::mapErrorClean()
{
    // open bam file
    BamMultiReader bamReader;
    bamReader.Open(bamFiles);

    // genome dictionary
    RefVector genomeDict = bamReader.GetReferenceData();

    // scanning window
    vector<tuple<int,int,int>> windows;
    int numWindows = GenericRegionTools::toScanWindow(genomeDict, regionStrings, step, windows);

    // random integer generator
    random_device rd;
    mt19937 rg(rd());
    uniform_int_distribution<> dist(1,UINT_MAX);

    // get working directory path
    size_t CWD_SIZE = PATH_MAX;
    char CWD[CWD_SIZE];
    getcwd(CWD,CWD_SIZE);

    // make a temporary directory
    mode_t temporaryMode = S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH;
    string temporaryDir;
    temporaryDir = string(CWD)+"/temporary";
    mkdir(temporaryDir.c_str(), temporaryMode);
    // change the working directory to the temporary directory
    chdir(temporaryDir.c_str());


    unordered_set<string> errMapReads;
    // loop over windows
    omp_set_dynamic(0);
    omp_set_num_threads(numThreads);
    #pragma omp parallel for shared(rg,dist,errMapReads)
    for (int i=0; i<numWindows; i++)
    {
        int wId = get<0>(windows[i]);
        int wLp = get<1>(windows[i]);
        int wRp = get<2>(windows[i]);

        if (verbose>=1) Verbose("process subregion " + genomeDict[wId].RefName + ":" + to_string(wLp+1) + "-" +to_string(wRp));
        // generate random number
        unsigned int rn = dist(rg);
        // build a temporary directory
        string TEMP_DIR = temporaryDir+"/temporary_"+to_string(rn);
        // path to pyrotools
        string PATH_TO_PYROTOOLS = GetExecPath();
        // temporary shell script
        string temporaryScript="run_temporary_" + to_string(rn) + ".sh";
        // temporary result file
        string temporaryResult="temporary_"+to_string(rn)+".txt";
        // write the pipeline to a temporary local file
        stringstream subcmd,cmd;
        subcmd << "mkdir " << TEMP_DIR << "\\n"
               << "cd " << TEMP_DIR << "\\n"
               << PATH_TO_PYROTOOLS << "\./xgsutils/xgsutils "
               << "xBamMapErrorRead "
               << "--roi " << genomeDict[wId].RefName << ":" << (wLp+1) << "-" << wRp << " "
               << "--ff " << alnFlagMarker << " "
               << "--mq " << mapQualThres << " "
               << "--bt2dict " << bt2dict << " "
               << genomeFile << " " << bamFiles[0] << " "
               << "temporary_" << rn << ".txt" << "\\n"
               << "cd .." << "\\n"
               << "if [ ! -z \"" << TEMP_DIR << "/" << temporaryResult << " \"]" << "\\n"
               << "then" << "\\n"
               << "  mv " << TEMP_DIR << "/" << temporaryResult << " " << temporaryDir << "/\\n"
               << "fi" << "\\n"
               << "rm -rf " << TEMP_DIR;
        cmd << "echo \"" << subcmd.str() << "\" > " << temporaryScript;
        system(cmd.str().c_str());
        // run the temporary shell script
        string runScript = "bash ./" + temporaryScript;
        system(runScript.c_str());

        if (!FileExist(temporaryDir+"/"+temporaryResult))
                continue;

        // read the file and load the erroneously mapped reads
        ifstream infile(temporaryDir+"/"+temporaryResult);
        string line;
        while (getline(infile,line))
        {
            if (line.empty()) continue;
            string readname;
            stringstream readstream(line);
            readstream >> readname;
            errMapReads.emplace(readname);
        }
        infile.close();

        // remove the temporary file
        system(string("rm -rf "+temporaryDir+"/"+temporaryResult).c_str());
        system(string("rm -rf "+temporaryDir+"/"+temporaryScript).c_str());
    }

    // change back from the temporary directory to the parent directory
    chdir(CWD);
    // remove the temporary working directory
    system(string("rm -rf "+temporaryDir).c_str());

    if (verbose>=1) Verbose("output the bam alignments");

    // output new bam file
    BamAlignment aln;
    // bam reader
    bamReader.Open(bamFiles);
    // bam writer
    BamWriter bamWriter;
    if (compressMode==0)
        bamWriter.SetCompressionMode(BamWriter::CompressionMode::Compressed);
    if (compressMode==1)
        bamWriter.SetCompressionMode(BamWriter::CompressionMode::Uncompressed);
    bamWriter.Open(string("stdout"),bamReader.GetHeaderText(),genomeDict);
    // retrieve the alignment
    while(bamReader.GetNextAlignmentCore(aln))
    {
        string aname = aln.Name;
        if (aln.IsPaired()){
            if (aln.IsFirstMate())
                aname += ".1";
            else
                aname += ".2";
        }

        auto ptr = errMapReads.find(aname);
        if (ptr!=errMapReads.end()){
            aln.SetIsMapped(false);
        }
        bamWriter.SaveAlignment(aln);
    }

    // close all
    bamReader.Close();
    bamWriter.Close();


    return 0;
}

int MapErrorCleanTool::Run(int argc, char *argv[])
{
    // print help message if no arguments provided
    if (argc==2){ Help(); exit(0); }

    // parse the command line
    parseCommandLind(argc-1,argv+1);

    // run the program
    if (mapErrorClean()!=0)
        return EXIT_FAILURE;

    return 0;
}
