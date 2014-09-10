#ifndef GENERICDATASTATISTICS_H
#define GENERICDATASTATISTICS_H

#include <string>
#include <vector>
#include <iostream>
using namespace std;

#include "api/BamMultiReader.h"
#include "utils/bamtools_fasta.h"
using namespace BamTools;

#include "GenericSequenceGlobal.h"
#include "GenericReadBins.h"
#include "GenericTool.h"
namespace GenericSequenceTools{


class GenericDataStatistics
{
    // ctor & dtor
    public:
        GenericDataStatistics(GenericReadBins &bins);

    public:
        void update(const string& alignRead, const string& alignGenome);
        void updateMatchMismatch(const string& alignRead, const string& alignGenome);
        void updateHomopolymerGap(const string& alignRead, const string& alignGenome);
        void updateDelete(const string& alignRead, const string& alignGenome);
        void updateInsert(const string& alignRead, const string& alignGenome);

    // I/O functions
    public:
        friend ostream& operator <<(ostream& output, GenericDataStatistics& gds)
        {
            gds.printMatchMismatch(output);
            gds.printHomopolymerGap(output);
            gds.printDelete(output);
            gds.printInsert(output);

            return output;
        }

        void printMatchMismatch(ostream& output);
        void printHomopolymerGap(ostream& output);
        void printDelete(ostream& output);
        void printInsert(ostream& output);

    public:
        // binning information
        GenericReadBins m_bins;

        //-------------------------------------------------------
        // error rate (exclude insertion/deletion)

        // mismatch counts at each bin
        VectorDouble m_binMismatchCount;
        // total counts including match and mismatch at each bin
        VectorDouble m_binCount;

        //-------------------------------------------------------
        // homopolymer gap

        // the number of homopolymers that are with deletions
        VectorDouble m_homopolymerDelete;
        // the number of homopolymers that are with insertions
        VectorDouble m_homopolymerInsert;
        // the number of homopolymers that are in data
        VectorDouble m_homopolymerCount;

        //--------------------------------------------------------
        // deletion in data

        // the number of homopolymer-in-sites that are deleted
        Matrix2Double m_homopolymerSiteDelete;
        // the number of homopolymer-in-sites that are in data
        Matrix2Double m_homopolymerSiteCount;

        //--------------------------------------------------------
        // insertion in data

        // the number of homopolymers that are adjacent to insertions,
        // and have the same nucleotide
        // homopolymer is in read
        VectorDouble m_homopolymerNextToInsert;

        // the number of homopolymers that are closed to insertions
        // at left, but separated by an interseptal homopolymer
        // homopolymer is in reference
        VectorDouble m_homopolymerLeftCloseToInsert;
        // the number of homopolymers that are closed to insertions
        // at right, but separated by an interseptal homopolymer
        // homopolymer is in reference
        VectorDouble m_homopolymerRightCloseToInsert;


};


class DataStatisticsTool: public GenericAbstractTool
{
    // ctor
    public:
        DataStatisticsTool(GenericReadBins &bins);

    // interface
    public:
        int Run(int argc, char* argv[]);
        int Help();

    public:
        int Execute();

    public:
        bool HasInputFiles;
        vector<string> InputFiles;

        bool HasInputFileList;
        string InputFileList;

    public:
        bool HasGenomeFilename;
        string GenomeFilename;

    public:
        bool HasRegion;
        string Region;

    public:
        ostream        output;
        BamMultiReader bamReader;
        BamRegion      bamRegion;
        RefVector      genomeReferences;
        Fasta          genomeFasta;
        GenericDataStatistics statistics;
};

}

#endif // DATASTATISTICS_H
