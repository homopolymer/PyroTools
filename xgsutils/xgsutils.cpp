#include <string>
#include <iostream>
#include <limits.h>
#include <unordered_map>
using namespace std;

#include <unistd.h>
#ifdef __APPLE__
    #include <mach-o/dyld.h>
#endif

#include <stdlib.h>
#include <libgen.h>

// subprogram path
unordered_map<string,string> CmdReg=
{
    // alnutils
    {"xEmbossNeedle",                   "alnutils/xEmbossNeedle"},
    {"xEmbossWater",                    "alnutils/xEmbossWater"},
    {"xPairwiseIdentity",               "alnutils/xPairwiseIdentity"},

    // asmutils
    {"xCeleraAssemblerUnitig",          "asmutils/xCeleraAssemblerUnitig"},
    {"xCeleraAssemblerUnitigIllumina",  "asmutils/xCeleraAssemblerUnitigIllumina"},
    {"xAssemblyContigN50",              "asmutils/xAssemblyContigN50"},
    {"xAssemblyContigN90",              "asmutils/xAssemblyContigN90"},
    {"xAssemblyComplexity",             "asmutils/xAssemblyComplexity"},

    // bamutils
    {"xBamAverageDepth",                "bamutils/xBamAverageDepth"},
    {"xBamAverageDepthEachChr",         "bamutils/xBamAverageDepthEachChr"},
    {"xBamExtractPairEndReadIntoFastq", "bamutils/xBamExtractPairEndReadIntoFastq"},
    {"xBamExtractReadToFasta",          "bamutils/xBamExtractReadToFasta"},
    {"xBamExtractReadToFastq",          "bamutils/xBamExtractReadToFastq"},
    {"xBamFilterByNM",                  "bamutils/xBamFilterByNM"},
    {"xBamFromSamToBam",                "bamutils/xBamFromSamToBam"},
    {"xBamMapErrorRead",                "bamutils/xBamMapErrorRead"},
    {"xBamReadCountEachChr",            "bamutils/xBamReadCountEachChr"},
    {"xBamToBedByDepth",                "bamutils/xBamToBedByDepth"},
    {"xBamToBedByRead",                 "bamutils/xBamToBedByRead"},
    {"xBamMarkClipSite",                "bamutils/xBamMarkClipSite"},
    {"xBamExtractClipSite",             "bamutils/xBamExtractClipSite"},
    {"xBamQuickFindLargeDelete",        "bamutils/xBamQuickFindLargeDelete"},
    
    // bedutils
    {"xBedConcatenate",                 "bedutils/xBedConcatenate"},
    {"xBedSortByGenomeCoordinate",      "bedutils/xBedSortByGenomeCoordinate"},

    // faxutils
    {"xFastaSequenceCount",             "faxutils/xFastaSequenceCount"},
    {"xFastaSequenceLength",            "faxutils/xFastaSequenceLength"},
    {"xFastaMutate",                    "faxutils/xFastaMutate"},
    {"xFastqSequenceCount",             "faxutils/xFastqSequenceCount"},
    {"xFastqFilterBlacklist",           "faxutils/xFastqFilterBlacklist"},
    {"xFastqKeepWhitelist",             "faxutils/xFastqKeepWhitelist"},

    // gtrutils
    {"xTrFinderMicroSatellite",         "gtrutils/xTrFinderMicroSatellite"},
    {"xTrFinderMiniSatellite",          "gtrutils/xTrFinderMicroSatellite"},
    {"xTrFinderSatellite",              "gtrutils/xTrFinderMicroSatellite"}
};

// help message
int help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    xgsutils COMMAND [OPTIONS]" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    xgsutils is a routine toolkit of sequencing data processing" << endl;
    cerr << "" << endl;
    cerr << "COMMANDS" << endl;
    cerr << "--alnutils" << endl;
    cerr << "    xEmbossNeedle                      run Emboss needle alignment" << endl;
    cerr << "    xEmbossWater                       run Emboss water alignment" << endl;
    cerr << "    xPairwiseIdentity                  compute the identity of the pairwise alignment" << endl;
    cerr << "" << endl;
    cerr << "--asmutils" << endl;
    cerr << "    xCeleraAssemblerUnitig             run Celera Assembler" << endl;
    cerr << "    xCeleraAssemblerUnitigIllumina     run Celera Assembler" << endl;
    cerr << "    xAssemblyContigN50                 compute the N50 statistics of the contigs" << endl;
    cerr << "    xAssemblyContigN90                 compute the N90 statistics of the contigs" << endl;
    cerr << "    xAssemblyComplexity                estimate the complexity of assembly in a region" << endl;
    cerr << "" << endl;
    cerr << "--bamutils" << endl;
    cerr << "    xBamAverageDepth                   compute the average sequencing depth" << endl;
    cerr << "    xBamAverageDepthEachChr            compute the average sequencing depth for each chromosome/contig" << endl;
    cerr << "    xBamExtractPairEndReadIntoFastq    extrac sequencing reads out of the bam file(s)" << endl;
    cerr << "    xBamExtractReadToFasta             extract sequencing reads out of the bam file(s)" << endl;
    cerr << "    xBamExtractReadToFastq             extract sequencing reads out of the bam file(s)" << endl;
    cerr << "    xBamFilterByNM                     filter sequencing reads by NM field" << endl;
    cerr << "    xBamFromSamToBam                   convert a sam file to bam format" << endl;
    cerr << "    xBamMapErrorRead                   list the erroneously mapped reads" << endl;
    cerr << "    xBamReadCountEachChr               count the number of sequencing reads on each chromosome/contig" << endl;
    cerr << "    xBamToBedByDepth                   convert the covered regions in the bam file to bed format" << endl;
    cerr << "    xBamToBedByRead                    convert sequencing reads in the bed file" << endl;
    cerr << "    xBamMarkClipSite                   mark the hard/soft clip sites on the genome" << endl;
    cerr << "    xBamExtractClipSite                extract the hard/soft clipping sequences" << endl;
    cerr << "    xBamQuickFindLargeDelete           quickly find the large deletions in the region of interest" << endl;
    cerr << "" << endl;
    cerr << "--bedutils" << endl;
    cerr << "    xBedConcatenate                    concatenate multiple bed files" << endl;
    cerr << "    xBedSortByGenomeCoordinate         sort a bed file by genomic coordinate" << endl;
    cerr << "" << endl;
    cerr << "--faxutils" << endl;
    cerr << "    xFastaSequenceCount                count the number of fasta sequences in the file" << endl;
    cerr << "    xFastaSequenceLength               list the length of fasta sequences in the file" << endl;
    cerr << "    xFastaMutate                       mutate the genome sequence" << endl;
    cerr << "    xFastqSequenceCount                count the number of fastq sequences in the file" << endl;
    cerr << "    xFastqFilterBlacklist              remove the reads in the blacklist" << endl;
    cerr << "    xFastqKeepWhitelist                print the reads in the whitelist" << endl;
    cerr << "" << endl;
    cerr << "--gtrutils" << endl;
    cerr << "    xTrFinderMicroSatellite            annotate the microsatellites in the sequence(s)" << endl;
    cerr << "    xTrFinderMiniSatellite             annotate the minisatellites in the sequence(s)" << endl;
    cerr << "    xTrFinderSatellite                 annotate the satellites in the sequence(s)" << endl;
    cerr << "" << endl;
    cerr << "Use xgsutils help command or xgsutils command -h to learn the specific command" << endl;
    return 0;
}

// get the path direct to xgsutils
// reference: http://stackoverflow.com/questions/143174/how-do-i-get-the-directory-that-a-program-is-running-from
string getexepath()
{
    string path;
    char result[PATH_MAX];
    uint32_t size=sizeof(result);
    #ifdef __linux
        ssize_t count = readlink( "/proc/self/exe", result, size );
        path = string( result, (count > 0) ? count : 0 );
    #elif __APPLE__
        _NSGetExecutablePath( result, &size );
        path = string(result);
    #endif
    return string(dirname((char*)path.c_str()));
}

int runSubCommand(string pathToXgsutils, string cmd, string arguments)
{
    string fullcmd;
    fullcmd=pathToXgsutils + "/" + CmdReg[cmd.c_str()] + arguments;
    system(fullcmd.c_str());
    
    return 0;
}

int main(int argc, char *argv[])
{
    // print help message if no command is specified
    if (argc==1) {help();exit(0);}
    

    // path to xgsutils
    string pathToXgsutils=getexepath();
    
    // set the arguments
    string arguments;
    for (int i=2; i<argc; i++)
    {
        arguments += " ";
        arguments += argv[i];
    }
    
    // run the command
    runSubCommand(pathToXgsutils, string(argv[1]), arguments);

    return 0;
}
