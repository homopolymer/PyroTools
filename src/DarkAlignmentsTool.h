#ifndef DARKALIGNMENTSTOOL_H
#define DARKALIGNMENTSTOOL_H

#include "GenericTool.h"
using namespace GenericSequenceTools;

#include <string>
using namespace std;

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "utils/bamtools_fasta.h"
using namespace BamTools;

namespace GenericSequenceTools
{

class DarkAlignmentsTool : public GenericAbstractTool
{
    public:
        DarkAlignmentsTool();

    public:
        int Run(int argc, char *argv[]);
        int Help();

    public:
        int Execute();

    public:
        bool HasInput;
        bool HasOutput;
        bool HasFasta;
        bool HasCompress;
        bool HasQualityThreshold;

        string InputFile;
        string OutputFile;
        string FastaFile;
        int    QualityThreshold;
        string strQualityThreshold;

    public:
        BamReader Reader;
        BamWriter Writer;
        RefVector GenomeSequences;
        Fasta     GenomeFasta;

};

}
#endif // DARKALIGNMENTSTOOL_H
