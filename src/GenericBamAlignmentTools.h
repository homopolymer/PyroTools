#ifndef GENERICBAMALIGNMENTTOOLS_H
#define GENERICBAMALIGNMENTTOOLS_H

#include "utils/bamtools_fasta.h"
#include "api/BamAux.h"
#include "api/BamAlignment.h"
using namespace BamTools;

#include <vector>
#include <iostream>
using namespace std;

#include "GenericSequenceGlobal.h"
#include "GenericTool.h"



namespace GenericSequenceTools{


typedef vector<CigarOp> Cigar;
typedef vector<Cigar> VectorCigar;

typedef string BamMD;
typedef char   BamMdType;
static BamMdType BamMdMatch    = 'M';
static BamMdType BamMdMismatch = 'X';
static BamMdType BamMdDelete   = 'D';
struct BamMdOp
{
    public:
        BamMdOp(const BamMdType& t='\0', const int& l=0, const string& s=string())
            : Type(t)
            , Length(l)
        {
            if (!s.empty())
            {
                Nucleotides.assign(s.begin(), s.end());
            }
            else
            {
                Nucleotides.assign(Length, Amb);
            }
        }

    public:
        BamMdType Type;
        int       Length;
        vector<Alpha> Nucleotides;
};
typedef vector<BamMdOp> BamMdOpArray;



typedef bool RunInDelShift;
static  bool DoInDelShift   = true;
static  bool SkipInDelShift = false;



typedef Alpha BlockState;
typedef string BlockType;
static BlockType BLOCK_MATCH      = "Match";
static BlockType BLOCK_MISMATCH   = "Mismatch";
static BlockType BLOCK_OVERCALL   = "OverCall";
static BlockType BLOCK_UNDERCALL  = "UnderCall";
static BlockType BLOCK_INSERT     = "Insert";
static BlockType BLOCK_DELETE     = "Delete";
struct BamAlignmentBlock
{

    public:

        BamAlignmentBlock( const string& alignRead = string(),
                           const string& alignGenome = string(),
                           const int& readPosition = 0,
                           const int& genomePosition = 0
                         )
            : m_alignRead(alignRead)
            , m_alignGenome(alignGenome)
            , m_read(string())
            , m_genome(string())
            , m_readPosition(readPosition)
            , m_genomePosition(genomePosition)
        {

            string::iterator iter;

            if (!m_alignRead.empty())
            {
                iter = m_alignRead.begin();
                for (; iter!=m_alignRead.end(); iter++)
                {
                    if ((*iter)!='-')
                        m_read += *iter;
                }
            }

            if (!m_alignGenome.empty())
            {
                iter = m_alignGenome.begin();
                for (; iter!=m_alignGenome.end(); iter++)
                {
                    if ((*iter)!='-')
                        m_genome += *iter;
                }
            }
        }

    public:
    BlockType Type()
    {
        BlockType blockType;

        if (m_read == m_genome)
            blockType = BLOCK_MATCH;

        else
        {
            if (m_read.empty())
                blockType = BLOCK_DELETE;

            else if (m_genome.empty())
                blockType = BLOCK_INSERT;

            else if (m_read.length() > m_genome.length())
                blockType = BLOCK_OVERCALL;

            else if (m_read.length() < m_genome.length())
                blockType = BLOCK_UNDERCALL;

            else
                blockType = BLOCK_MISMATCH;
        }

        return blockType;
    }

    int genomeLength()
    {
        return m_genome.length();
    }

    int readLength()
    {
        return m_read.length();
    }

    public:
        string m_alignRead;
        string m_alignGenome;
        string m_read;
        string m_genome;
        int    m_readPosition;
        int    m_genomePosition;
};



class GenericBamAlignmentTools
{
    // constructor
    public:
        GenericBamAlignmentTools() {}

    // utilities
    public:

        //------------------------------------------------------
        // Insert/Delete Shift

        // move insertion/deletion to right
        static void leftShiftInDel(BamAlignment& alignObj);

        // shift the InDel in the given alignment
        static void shiftInDel(string& alignRead, string& alignGenome);

        // there is InDel in the alignment
        static bool hasInDel(BamAlignment& alignObj);

        // whether InDel in the alignment should be shifted
        static bool needShiftInDel(BamAlignment& alignObj);


        // move InDel to left
        static void moveInDelLeft(string& alignReadSeq, string& alignGenomeSeq);
        // move InDel to right
        static void moveInDelRight(string& alignReadSeq, string& alignGenomeSeq);



        // ------------------------------------------------
        // MD and Cigar Manipulation

        // get MD of the alignment
        static void getBamAlignmentMD(BamAlignment& alignObj, BamMdOpArray& bamMdOps);

        // get Cigar of the alignment
        static void getBamAlignmentCigar(BamAlignment& alignObj, Cigar& cigar);

        // get genome sequence from MD and Cigar
        static void getBamAlignmentGenome(BamAlignment& alignObj, string& genome);




        // compute new Cigar from the alignment
        static void calculateCigar(const string& alignRead, const string& alignGenome, Cigar& cigar);

        // compute new MD from the alignment
        static void calculateMD(const string& alignRead, const string& alignGenome, BamMD& md);


        // ---------------------------------------------------
        // Alignment Manipulation


        // construct the alignment from the given unaligned read and genome, following
        // the give cigar sequence. Exclude the soft clips
        static void getAlignmentSequences(const string& read, const string& genome,
                                          Cigar& cigar, string& alignRead, string& alignGenome);

        static void getAlignmentSequences(BamAlignment& alignObj,
                                          string& alignRead, string& alignGenome);

        // divide the alignment into homopolymer-wise blocks
        static void divideAlignmentToBlocks(const string& alignRead, const string& alignGenome,
                                            vector<BamAlignmentBlock>& blocks);



        // retrieve local alignment in the given alignment, given the start position and size
        static int getLocalAlignment(BamAlignment& alignObj, int position, int size,
                                      string& localRead, string& localGenome,
                                      Cigar& cigar=*static_cast<Cigar*>(0),
                                      BamMD& md=*static_cast<BamMD*>(0),
                                      int& numMismatch=*static_cast<int*>(0),
                                      int& numInDel=*static_cast<int*>(0));

        // retrieve local alignment in the given alignment, given the start position and size
        static int getLocalAlignment(BamAlignment& alignObj, int position, int size,
                                      string& localRead, string &localReadQual,
                                      string& localGenome);

        // retrieve local alignment in the given alignment, given the start position and size,
        // and also clip away the indel at both the head and tail ends
        static void getLocalAlignmentWithEndClip(BamAlignment &alignObj, int position, int size,
                                                 string &localRead, string &localGenome,
                                                 int &startPos, int &endPos);

        // change a part of an alignment, given the start position and size, given new alignment
        static void changeLocalAlignment(BamAlignment& alignObj, int position, int size,
                                         string& alnLocalRead, string& alnLocalGenome,
                                         BamAlignment& newAlignObj);

        // ----------------------------------------------
        // Alignment Information

        // name of the alignment
        static string getBamAlignmentName(BamAlignment& alignObj);

        // ID of the alignment
        static string getBamAlignmentID(BamAlignment& alignObj);

        // the number of matches in the alignment
        static int numBamAlignmentMatches(BamAlignment& alignObj);

        // the number of mismatches in the alignment
        static int numBamAlignmentMismatches(BamAlignment& alignObj);

        // retrieve the mismatches in the alignment
        static void getBamAlignmentMismatches(BamAlignment&   alignObj,
                                              vector<long>&   genomePositions,
                                              vector<char>&   referenceAlleles = *static_cast<vector<char>*>(0),
                                              vector<char>&   readAlleles = *static_cast<vector<char>*>(0),
                                              vector<double>& readAlleleQualities = *static_cast<vector<double>*>(0)
                                              );

        // the number of deletions in the alignment
        static int numBamAlignmentDeletes(BamAlignment& alignObj);

        // retrieve the deletions in the alignment
        static void getBamAlignmentDeletes(BamAlignment&   alignObj,
                                           vector<long>&   genomePositions,
                                           vector<string>& deletedSequences);

        // the number of insertions in the alignment
        static int numBamAlignmentInserts(BamAlignment& alignObj);

        // retrieve the insertions in the alignment
        static void getBamAlignmentInserts(BamAlignment&           alignObj,
                                           vector<long>&           genomePositions,
                                           vector<string>&         insertedSequences,
                                           vector<vector<double>>& insertedSequenceQualities);

        // the length of aligned sequence in the alignment
        static int getBamAlignmentReadLength(BamAlignment& alignObj);

        // the length of maximal gap in the alignment
        static int maxGapLength(BamAlignment& alignObj);
        static int maxGapLength(string& alignRead, string& alignGenome);
        static int maxGapLength(Cigar& cigar);

        // the number of gap in the alignment
        static int numGap(string& alignRead, string& alignGenome);
        // the number of insertion in the alignment
        static int numInsert(string& alignRead, string& alignGenome);
        // the number of deletion in the alignment
        static int numDelete(string& alignRead, string& alignGenome);

        // the number of mismatch in the alignment
        static int numMismatch(string& alignRead, string& alignGenome);
        static int numMatch(string& alignRead, string& alignGenome);

        // check it is good alignment
        static bool goodAlignment(BamAlignment& alignObj, bool keepDuplicate=false);
        static bool isDuplicateAlignment(BamAlignment& alignObj);
        static bool isUnmappedAlignment(BamAlignment& alignObj);
        static bool isSecondaryAlignment(BamAlignment& alignObj);
        static bool isSupplementaryAlignment(BamAlignment& alignObj);
        static bool isQcFailedAlignment(BamAlignment& alignObj);
        static bool isFilteredByFlag(BamAlignment& alignObj, int flag);

        // check read length
        static bool validReadLength(BamAlignment& alignObj, double minReadLength=100);
        // check map quality
        static bool validMapQuality(BamAlignment& alignObj, double minMapQuality=10);
        // check read quality
        static bool validReadQuality(BamAlignment& alignObj, double minReadQuality=10);
        // check mismatch ratio
        static bool validReadIdentity(BamAlignment& alignObj, double maxMismatchFrac=0.1);

        // BAM Cigar
        static bool isSoftClip(char op);
        static bool isHardClip(char op);
        static bool isClip(char op);
        static bool isMatch(char op);
        static bool isMismatch(char op);
        static bool isInsert(char op);
        static bool isDelete(char op);
        static bool isGap(char op);
        static bool isNotGap(char op);

        // convert cigar to string
        static void convertCigarToString(Cigar& cigar, string& cigarStr);

    // print
    public:
        static void printBamAlignmentCigar(BamAlignment& alignObj);
        static void printBamAlignmentMD(BamAlignment& alignObj);
        static void printBamAlignment(RefVector& genomeDict, BamAlignment& alignObj);

    // internal functions
    private:
        static void getAlignmentBlocks(const BamAlignment& alignObj,
                                       VectorInteger& blockReadPositions,
                                       VectorInteger& blockReadSizes,
                                       VectorInteger& blockGenomePositions,
                                       VectorInteger& blockGenomeSizes,
                                       VectorCigar&   blockCigars,
                                       vector<RunInDelShift>& blockRuns);
};

class CropBamTool : public GenericAbstractTool
{
public:
    CropBamTool();

public:
    int Help();
    int Run(int argc, char *argv[]);

public:
    int parseCommandLine(int argc, char *argv[]);
    int CropBam();

public:
    vector<string> bamFiles;
    vector<string> regionStrings;

    int mapQualThres;
    int readLenThres;
    int segmentLenThres;
    int alnFlagMarker;
    float alnIdenThres;
    string outFormat;
    bool keepClip;
    bool useUnique;
    int thresFreq;
    string outFreq;
    int numThreads;
    int verbose;

};

class MapErrorCleanTool : public GenericAbstractTool
{
public:
    MapErrorCleanTool();

public:
    int Help();
    int Run(int argc, char *argv[]);

public:
    int parseCommandLind(int argc, char *argv[]);

public:
    int mapErrorClean();

public:
    string genomeFile;
    vector<string> bamFiles;
    vector<string> regionStrings;

    int compressMode;
    int step;
    int alnFlagMarker;
    int mapQualThres;
    int readLenThres;
    string bt2dict;

    int numThreads;
    int verbose;
};
}   // namespace

#endif // GENERICBAMALIGNMENTTOOLS_H
