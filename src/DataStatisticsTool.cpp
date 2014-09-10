#include "GenericBamAlignmentTools.h"
#include "GenericDataStatistics.h"
using namespace GenericSequenceTools;

#include "api/BamAlignment.h"
#include "utils/bamtools_options.h"
#include "utils/bamtools_utilities.h"
using namespace BamTools;

#include <iostream>
using namespace std;


DataStatisticsTool::DataStatisticsTool(GenericReadBins &bins)
    : HasInputFiles(false)
    , HasInputFileList(false)
    , HasGenomeFilename(false)
    , HasRegion(false)
    , output(cout.rdbuf())
    , statistics(bins)
{
    // set program details
    BamTools::Options::SetProgramInfo("GenericSequenceTools errstat",
                                      "count the frequency of errors in BAM file(s)",
                                      "[-region <REGION>] -ref <REFERENCE> -bam <filename> -bam <filename> ...");

    // set options
    BamTools::OptionGroup* IO_Opts = BamTools::Options::CreateOptionGroup("Input Options");
    BamTools::Options::AddValueOption("-bam", "BAM filename", "the input BAM file(s)", "", HasInputFiles, InputFiles, IO_Opts);
    BamTools::Options::AddValueOption("-list", "filename", "the input BAM file list, one file per line", "", HasInputFileList, InputFileList, IO_Opts);
    BamTools::Options::AddValueOption("-ref", "FASTA filename", "FASTA reference file", "", HasGenomeFilename, GenomeFilename, IO_Opts);
    BamTools::Options::AddValueOption("-region", "REGION", "genomic region of interest", "", HasRegion, Region, IO_Opts);

}

int DataStatisticsTool::Run(int argc, char *argv[])
{
    // parse command line arguments
    BamTools::Options::Parse(argc, argv, 1);

    // -------------------------------
    // BAM file(s) input

    // BAM file(s) must be provided
    if (!HasInputFiles && !HasInputFileList){
        cerr << "GenericSequenceTools errstat ERROR: BAM file(s) not provided... Aborting" << endl;
        return false;
    }

    // add files in the filelist to input file list
    if (HasInputFileList){

        ifstream filelist(InputFileList.c_str(), ios::in);
        if (!filelist.is_open()){
            cerr << "GenericSequenceTools errstat ERROR: could not open input BAM file list... Aborting" << endl;
            return false;
        }

        string line;
        while (getline(filelist, line))
            InputFiles.push_back(line);
    }

    // open input files
    if (!bamReader.Open(InputFiles)){
        cerr << "GenericSequenceTools errstat ERROR: could not open input BAM file(s)... Aborting" << endl;
        return false;
    }

    // if input is not stdin and a region is provided, look for index files
    if ((HasInputFiles || HasInputFileList) && HasRegion){
        if (!bamReader.LocateIndexes()){
            cerr << "GenericSequenceTools errstat ERROR: could not locate index file(s)... Aborting" << endl;
            return false;
        }
    }

    // retrieve reference data
    genomeReferences = bamReader.GetReferenceData();

    // set region if specified
    BamTools::BamRegion bamRegion;
    if (HasRegion){
        if (BamTools::Utilities::ParseRegionString(Region, bamReader, bamRegion)){

            if (bamReader.HasIndexes()){
                if (!bamReader.SetRegion(bamRegion)){
                    cerr << "GenericSequenceTools errstat ERROR: set region failed... Aborting" << endl;
                    bamReader.Close();
                    return false;
                }
            }

        }else{
            cerr << "GenericSequenceTools errstat ERROR: could not parse REGION: " << Region << endl;
            bamReader.Close();
            return false;
        }
    }

    //---------------------------------------------
    // Fasta file input

    // Fasta file must be provided
    if (GenomeFilename.empty()){
        cerr << "GenericSequenceTools errstat ERROR: Fasta file not provided... Aborting" << endl;
        return false;
    }

    // check Fasta index
    string GenomeIndexFilename = "";
    if (BamTools::Utilities::FileExists(GenomeFilename+".fai"))
        GenomeIndexFilename = GenomeFilename + ".fai";

    // open Fasta file
    if (!genomeFasta.Open(GenomeFilename, GenomeIndexFilename)){
        cerr << "GenericSequenceTools errstat ERROR: " << GenomeFilename << " "
             << "can not be open... Aborting" << endl;
        return false;
    }


    // execution
    if (Execute())
        return 1;

    return 0;
}

int DataStatisticsTool::Execute()
{
    // iterate over reads in BAM file(s)
    BamAlignment alignObj;
    while(bamReader.GetNextAlignment(alignObj))
    {
        if (alignObj.IsDuplicate()) continue;
        if (alignObj.IsFailedQC()) continue;
        if (!alignObj.IsMapped()) continue;
        if (!alignObj.IsPrimaryAlignment()) continue;
        if (alignObj.IsPaired() && !alignObj.IsProperPair()) continue;
        if (alignObj.IsPaired() && !alignObj.IsMateMapped()) continue;
        if (!alignObj.HasTag("MD")) continue;

//        // debug
//        GenericBamAlignmentTools::printBamAlignmentCigar(alignObj);
//        GenericBamAlignmentTools::printBamAlignmentMD(alignObj);

        // shift InDel
        GenericBamAlignmentTools::leftShiftInDel(alignObj);

//        // debug
//        GenericBamAlignmentTools::printBamAlignmentCigar(alignObj);
//        GenericBamAlignmentTools::printBamAlignmentMD(alignObj);

        // get the alignment sequences
        string alignRead;
        string alignGenome;
        GenericBamAlignmentTools::getAlignmentSequences(alignObj, alignRead, alignGenome);

        // update the statistics
        statistics.update(alignRead, alignGenome);
    }


    // print to screen
    cout << statistics << endl;
//    statistics.printMatchMismatch();

    // close BAM reader
    bamReader.Close();

    // close Fasta
    genomeFasta.Close();

    return 1;
}

int DataStatisticsTool::Help()
{
    BamTools::Options::DisplayHelp();
    return 0;
}
