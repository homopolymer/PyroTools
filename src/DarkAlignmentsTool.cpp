#include "DarkAlignmentsTool.h"
#include "GenericBamAlignmentTools.h"
using namespace GenericSequenceTools;

#include "utils/bamtools_options.h"
#include "utils/bamtools_utilities.h"
using namespace BamTools;

#include <string>
using namespace std;

DarkAlignmentsTool::DarkAlignmentsTool()
    : HasInput(false)
    , HasOutput(false)
    , HasFasta(false)
    , HasCompress(false)
    , HasQualityThreshold(false)
    , InputFile(string())
    , OutputFile(string())
    , FastaFile(string())
    , QualityThreshold(0)
{
    // set program details
    BamTools::Options::SetProgramInfo("GenericSequenceTools zeroaln",
                                      "print the alignments with zero or small mapping quality score",
                                      "[OPTIONS] -ref <REFERENCE> -bam <INPUT_FILE> -out <OUTPUT_FILE>");

    // set options
    BamTools::OptionGroup* IO_Opts = BamTools::Options::CreateOptionGroup("Input and Output Options");
    BamTools::Options::AddValueOption("-bam", "BAM filename", "the input BAM file", "Z", HasInput, InputFile, IO_Opts);
    BamTools::Options::AddValueOption("-ref", "FASTA filename", "FASTA reference file", "", HasFasta, FastaFile, IO_Opts);
    BamTools::Options::AddValueOption("-out", "Output filename", "the output BAM file", "", HasOutput, OutputFile, IO_Opts);

    // set options
    BamTools::OptionGroup* OT_Opts = BamTools::Options::CreateOptionGroup("Additional Options");
    BamTools::Options::AddOption("-c", "output in compressed format", HasCompress, OT_Opts);
    BamTools::Options::AddValueOption("-q", "Quality score", "the threshold of mapping quality", "", HasQualityThreshold, strQualityThreshold, OT_Opts);
}


int DarkAlignmentsTool::Run(int argc, char *argv[])
{
    // parse the command line arguments
    BamTools::Options::Parse(argc, argv, 1);

    // ------------------------------------------
    // input BAM file

    // make sure that input is provided
    if (!HasInput || InputFile.empty())
    {
        cerr << "GenericSequenceTools zeroaln ERROR: input file not set... Aborting" << endl;
        return EXIT_FAILURE;
    }

    // open input file
    if (!Reader.Open(InputFile))
    {
        cerr << "GenericSequenceTools zeroaln ERROR: could not open input file... Aborting" << endl;
        return EXIT_FAILURE;
    }

    // open input index file
    if (!Reader.LocateIndex())
    {
        cerr << "GenericSequenceTools zeroaln ERROR: could not open input index file... Aborting" << endl;
        return EXIT_FAILURE;
    }

    // retrieve reference data
    GenomeSequences = Reader.GetReferenceData();

    //---------------------------------------------------
    // Reference Fasta file

    // make sure that Fasta file is provided
    if (!HasFasta || FastaFile.empty())
    {
        cerr << "GenericSequenceTools zerolan ERROR: Fasta file not set... Aborting" << endl;
        return EXIT_FAILURE;
    }

    // check Fasta index file
    string FastaIndexFile = string();
    if (BamTools::Utilities::FileExists(FastaFile+".fai"))
    {
        FastaIndexFile = FastaFile+".fai";
    }

    if (!GenomeFasta.Open(FastaFile, FastaIndexFile))
    {
        cerr << "GenericSequenceTools zeroaln ERROR: could not open Fasta file... Aborting" << endl;
        return EXIT_FAILURE;
    }

    //-----------------------------------------------------
    // Output file

    // make sure that output file is provided
    if (!HasOutput || OutputFile.empty())
    {
        cerr << "GenericSequenceTools zeroaln ERROR: output file not set... Aborting" << endl;
        return EXIT_FAILURE;
    }

    // open output file
    if (!Writer.Open(OutputFile, Reader.GetHeaderText(), GenomeSequences))
    {
        cerr << "GenericSequenceTools zeroaln ERROR: failed to open output file... Aborting" << endl;
        return EXIT_FAILURE;
    }

    if (HasCompress)
    {
        Writer.SetCompressionMode(BamWriter::Compressed);
    }else
    {
        Writer.SetCompressionMode(BamWriter::Uncompressed);
    }

    // map quality score threshold
    if (HasQualityThreshold)
        QualityThreshold = stoi(strQualityThreshold);

    // execution
    if (Execute())
        return 1;

    return 0;
}

int DarkAlignmentsTool::Help()
{
    BamTools::Options::DisplayHelp();
    return 0;
}


int DarkAlignmentsTool::Execute()
{
    BamAlignment alignObj;

    while(Reader.GetNextAlignment(alignObj))
    {
        // skip if map quality is higher than threshold
        if (alignObj.MapQuality > QualityThreshold)
            continue;

        // skip if read has no mismatch or too much mismatch
        int numMatch    = GenericBamAlignmentTools::numBamAlignmentMatches(alignObj);
        int numMismatch = GenericBamAlignmentTools::numBamAlignmentMismatches(alignObj);

        if (numMismatch == 0)
            continue;
        if (numMismatch >= 0.25*numMatch)
            continue;

        // save into output
        Writer.SaveAlignment(alignObj);
    }

    // close input
    Reader.Close();

    // close Fasta
    GenomeFasta.Close();

    // close output
    Writer.Close();

    return 0;
}
