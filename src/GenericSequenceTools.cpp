#include "GenericReadBins.h"
#include "GenericTool.h"
#include "GenericDataStatistics.h"
#include "DarkAlignmentsTool.h"
#include "GenericProbabilisticAlignment.h"
#include "GenericIndividualSnpCall.h"
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

    if (arg == REALIGNER) return new ProbabilisticAlignmentTool;

    if (arg == SNPCALLER) return new IndividualSnpCallTool;

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
    cerr << endl;
    cerr << "usage: PyroTools [--help] COMMAND [ARGS]" << endl;
    cerr << endl;
    cerr << "Available PyroTools commands:" << endl;
    cerr << "" << endl;
    cerr << "\tGraphReAlign     Local realignment around homopolymers, dinucleotide repeats and MNP sites" << endl;
    cerr << "\tSnpCall          Haplotype-based SNP calling algorithm for sequencing data of an individual sample [coming soon]" << endl;
    cerr << "\tIndelCall        Haplotype-based Indel calling algorithm for sequencing data of an individual sample [coming soon]" << endl;
    cerr << "\tErrStat          Count the frequency of errors in data" << endl;
    cerr << endl;
    cerr << "See 'PyroTools help COMMAND' for more information on a specific command." << endl;
    cerr << endl;
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

