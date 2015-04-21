#include "GenericReadBins.h"
#include "GenericTool.h"
#include "GenericDataStatistics.h"
#include "DarkAlignmentsTool.h"
#include "GenericProbabilisticAlignment.h"
#include "GenericIndividualSnpCall.h"
#include "GenericGraphTools.h"
#include "SimpleVariantCallTool.h"
#include "GenericBamAlignmentTools.h"
using namespace GenericSequenceTools;

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

// PyroTools subtool names
static const string REALIGNER   = "GraphReAlign";
static const string SNPCALLER   = "SnpCall";
static const string INDELCALLER = "IndelCall";
static const string ERRORSTAT   = "ErrStat";
static const string CNSGRAPH    = "GraphConsensus";
static const string SIMPLEVC    = "SimpleVarCall";
static const string CROPBAM     = "CropBam";
static const string MAPERRCLN   = "MapErrorClean";
static const string LOCALSC     = "StrainCall";
static const string STRAINBIN   = "StrainBinning";

// PyroTools help/version constants
static const string HELP          = "help";
static const string LONG_HELP     = "--help";
static const string SHORT_HELP    = "-h";
static const string VERSION       = "version";
static const string LONG_VERSION  = "--version";
static const string SHORT_VERSION = "-v";

// determine if string is a help constant
static bool IsHelp(char* str){
    return (str==HELP ||
            str==LONG_HELP ||
            str==SHORT_HELP);
}

//// determine if string is a version constant
//static bool IsVersion(char* str){
//    return (str==VERSION ||
//            str==LONG_VERSION ||
//            str==SHORT_VERSION);
//}

// subtool factory method
GenericAbstractTool* CreateTool(const string& arg){
    // determine tool type based on arg
    if (arg == ERRORSTAT)
    {
        GenericReadBins bins(0, 15, 5);
        return new DataStatisticsTool(bins);
    }

    // realigner
    if (arg == REALIGNER) return new ProbabilisticAlignmentTool;

    // snp calling
    if (arg == SNPCALLER) return new IndividualSnpCallTool;
    // simple variant calling
    if (arg == SIMPLEVC)  return new SimpleVariantCallTool;

    // consensus graph
    if (arg == CNSGRAPH)  return new ConsensusGraphTool;

    // chop bam file
    if (arg == CROPBAM)   return new CropBamTool;

    // clean map error
    if (arg == MAPERRCLN) return new MapErrorCleanTool;

    // local strain call
    if (arg == LOCALSC) return new LocalStrainCallTool;

    // strain binning
    if (arg == STRAINBIN) return new StrainBinning;

    // unknown arg
    return 0;
}

// print help info
int Help(int argc, char* argv[]){
    // check for "PyroTools help command" to print tool-specific help message
    if (argc>2){
        // determine desired subtool
        GenericAbstractTool* tool = CreateTool(argv[2]);

        // if tool known, print its help on screen
        if (tool) return tool->Help();
    }

    // print general PyroTools help message
    cerr << "SYNOPSIS" << endl;
    cerr << "    PyroTools COMMAND [ARGS]" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    a toolkit to process the high throughput sequencing data" << endl;
    cerr << "" << endl;
    cerr << "COMMANDS" << endl;
    cerr << "--Variant calling programs" << endl;
    cerr << "    SimpleVarCall      profile the variants by using the heuristic method, parallel mode supported" << endl;
    cerr << "    SnpCall            detect the SNPs after alignment adjust and mapping error clean" << endl;
    cerr << "    IndelCall          detect the indels by by using the haplotype-based method" << endl;
    cerr << "" << endl;
    cerr << "--Strain calling programs" << endl;
    cerr << "    StrainBinning      aggregate long reads or assembled contigs/scaffolds into strain-" << endl;
    cerr << "                       specific bin" << endl;
    cerr << "    StrainCall         reconstruct the strains in a local region by using the viral/microbial" << endl;
    cerr << "                       sequencing data, and inclusively output the SNPs and/or Indels in the strains," << endl;
    cerr << "                       the size of the local region is recommended not larger than 500bp" << endl;
    cerr << "--Bias correction programs" << endl;
    cerr << "    ReAligner          correct the alignment bias" << endl;
    cerr << "    MapErrorClean      clean mapping errors" << endl;
    cerr << "" << endl;
    cerr << "--Haplotype reconstruction programs" << endl;
    cerr << "    GraphConsensus     output the top-k paths in the consensus graph" << endl;
    cerr << "    HaplotypeBuild     reconstruct the haplotype/consensus sequences" << endl;
    cerr << "" << endl;
    cerr << "--Error statistics programs" << endl;
    cerr << "    ErrStat            count the errors in data" << endl;
    cerr << "" << endl;
    cerr << "--BAM processing programs" << endl;
    cerr << "    CropBam            extract reads that might be partial in a region" << endl;
    cerr << endl;
    cerr << "See 'PyroTools help COMMAND' for more information on a specific command." << endl;
    return EXIT_SUCCESS;
}

//// print version info
//int Version(void){
//    stringstream versionStream("");
//    versionStream << SINGLECELLTOOLS_VERSION_MAJOR << "."
//                  << SINGLECELLTOOLS_VERSION_MINOR << "."
//                  << SINGLECELLTOOLS_VERSION_BUILD;

//    cout << endl;
//    cout << "SingleCellTools" << versionStream.str() << endl;
//    cout << "Single cell sequencing data toolkit" << endl;
//    cout << "Author: Feng Zeng" << endl;
//    cout << "(c) 2008-2014 Automation Department, Tsinghua University" << endl;
//    cout << "(c) 2014-now  Automation Department, Xiamen Uninversity" << endl;
//    cout << endl;

//    return EXIT_SUCCESS;
//}

int main(int argc, char* argv[])
{
    // just 'PyroTools'
    if (argc==1) return Help(argc, argv);

    // 'PyroTools help', 'PyroTools --help', or 'PyroTools -h'
    if (IsHelp(argv[1])) return Help(argc, argv);

//    // 'PyroTools version', 'PyroTools'
//    if (IsVersion(argv[1])) return Version();

    // determine desired subtool, run if found
    GenericAbstractTool* tool = CreateTool(argv[1]);
    if (tool) return tool->Run(argc, argv);

    // no tool matched, show help
    return Help(argc, argv);
}

