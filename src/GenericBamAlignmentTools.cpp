#include "GenericSequenceGlobal.h"
#include "GenericBamAlignmentTools.h"
#include "GenericFastaTools.h"
using namespace GenericSequenceTools;

#include "api/BamAux.h"
using namespace BamTools;

#include <set>
#include <tuple>
#include <string>
#include <vector>
#include <iterator>
#include <iostream>
#include <sstream>
using namespace std;

// Mark
// Right now, we assume that soft clip(s) just occurs at the begin or end of the alignment,
// and there is no soft clip(s) that occur innerly, and then divide the read into several
// alignment blocks.
// This situation is true for 454, Illumina, and Ion Torrent data.  However, it is not
// hold for PacBio data.
// This issue need to be fixed in future.
void GenericBamAlignmentTools::leftShiftInDel(BamAlignment &alignObj)
{

    if (!hasInDel(alignObj))
        return;

    if (!needShiftInDel(alignObj))
        return;

    // alignment block
    vector<BamMD> blockMDs;
    VectorCigar   blockCigars;
    VectorInteger blockReadPositions;
    VectorInteger blockReadSizes;
    VectorInteger blockGenomePositions;
    VectorInteger blockGenomeSizes;
    vector<RunInDelShift> blockRuns;

    getAlignmentBlocks(alignObj,
                       blockReadPositions, blockReadSizes,
                       blockGenomePositions, blockGenomeSizes,
                       blockCigars,
                       blockRuns);

    // the number of alignment block
    int numAlignmentBlocks = blockCigars.size();

    // iterate over blocks
    for (int i=0; i<numAlignmentBlocks; i++)
    {
        if (blockRuns[i]==SkipInDelShift)
        {
            BamMD md = string();
            blockMDs.push_back(md);
            continue;
        }

        string blockReadSequence;
        string blockGenomeSequence;

        // retrieve the sequences
        blockReadSequence = alignObj.QueryBases.substr(blockReadPositions[i], blockReadSizes[i]);
        getBamAlignmentGenome(alignObj, blockGenomeSequence);

        // reconstruct the alignment from the cigar
        string blockAlignRead;
        string blockAlignGenome;
        getAlignmentSequences(blockReadSequence, blockGenomeSequence, blockCigars[i],
                              blockAlignRead, blockAlignGenome);

        // shift the InDel
        shiftInDel(blockAlignRead, blockAlignGenome);

        // compute new Cigar
        Cigar cigar;
        calculateCigar(blockAlignRead, blockAlignGenome, cigar);
        blockCigars[i] = cigar;

        // compute new MD
        BamMD md;
        calculateMD(blockAlignRead, blockAlignGenome, md);
        blockMDs.push_back(md);
    }

    // set the cigar of the alignment to new
    Cigar NewCigar;
    for (int i=0; i<numAlignmentBlocks; i++)
    {
        Cigar cigar = blockCigars[i];

        NewCigar.insert(NewCigar.end(), cigar.begin(), cigar.end());
    }
    alignObj.CigarData.assign(NewCigar.begin(), NewCigar.end());

    // set the MD of the alignment to new
    BamMD NewMD = "";
    for (int i=0; i<numAlignmentBlocks; i++)
        NewMD += blockMDs[i];

    alignObj.EditTag("MD", "Z", NewMD);

}

bool GenericBamAlignmentTools::hasInDel(BamAlignment &alignObj)
{
    bool HasInDel = false;

    for (int i=0; i<alignObj.CigarData.size(); i++)
    {
        if ( isDelete(alignObj.CigarData[i].Type) ||
             isInsert(alignObj.CigarData[i].Type)
           )
        {
            HasInDel = true;
            break;
        }
    }

    return HasInDel;
}

bool GenericBamAlignmentTools::needShiftInDel(BamAlignment &alignObj)
{
    bool ShiftInDel = false;

    BamMdOpArray bamMdOps;
    Cigar        bamCigar;

    getBamAlignmentMD(alignObj, bamMdOps);
    bamCigar.assign(alignObj.CigarData.begin(), alignObj.CigarData.end());

    // assign each genomic position as match, mismatch, or delete
    vector<BamMdType> genomeMdx;
    vector<Alpha>     genomeAlpha;
    for (BamMdOpArray::iterator iter=bamMdOps.begin();
         iter!=bamMdOps.end(); ++iter)
    {
        genomeMdx.insert(genomeMdx.end(), (*iter).Length, (*iter).Type);
        genomeAlpha.insert(genomeAlpha.end(), (*iter).Nucleotides.begin(), (*iter).Nucleotides.end());
    }

    // iterate from read begin to read end
    int readPointer=0;
    int genomePointer=0;
    for (Cigar::iterator iter=bamCigar.begin(); iter!=bamCigar.end(); ++iter)
    {
        if (!isDelete(iter->Type))
            readPointer += (*iter).Length;

        if (!isClip(iter->Type) && !isInsert(iter->Type))
            genomePointer += (*iter).Length;

        if (isDelete(iter->Type))
        {
            if (genomeMdx[genomePointer]==BamMdMatch &&
                    genomeAlpha[genomePointer-1]==alignObj.QueryBases[readPointer])
            {
                ShiftInDel = true;
                break;
            }
        }

        if (isInsert(iter->Type))
        {
            if (genomeMdx[genomePointer]==BamMdMatch &&
                    alignObj.QueryBases[readPointer-1]==alignObj.QueryBases[readPointer])
            {
                ShiftInDel = true;
                break;
            }
        }
    }



    return ShiftInDel;
}

void GenericBamAlignmentTools::shiftInDel(string &alignRead, string &alignGenome)
{
    char   blockState;
    string blockRead;
    string blockGenome;
    vector<tuple<string,string>> alignBlocks;

    // divide the alignment into blocks
    bool blockChange;

    blockRead   = "";
    blockGenome = "";
    for (int i=0; i<alignRead.size(); i++)
    {
        // initialize
        if (i==0)
        {
            if (alignRead[i]!='-')
                blockRead = alignRead[i];

            if (alignGenome[i]!='-')
                blockGenome = alignGenome[i];

            blockChange = false;

            if (alignGenome[i]!='-')
                blockState = alignGenome[i];
            else
                blockState = alignRead[i];

            continue;
        }

        // block change
        if ((blockState!=alignRead[i] && alignRead[i]!='-') ||
            (blockState!=alignGenome[i] && alignGenome[i]!='-'))
        {
            blockChange = true;
        }

        // save previous block
        if (blockChange)
        {
            alignBlocks.push_back(tuple<string,string>(blockRead, blockGenome));

            blockRead   = "";
            blockGenome = "";

            blockChange = false;

            if (alignGenome[i]!='-')
                blockState = alignGenome[i];
            else
                blockState = alignRead[i];
        }

        // update block sequences
        if (alignRead[i]!='-')
            blockRead += alignRead[i];

        if (alignGenome[i]!='-')
            blockGenome += alignGenome[i];
    }

    if (!blockChange){
        alignBlocks.push_back(tuple<string,string>(blockRead, blockGenome));
    }

    // link blocks to build new alignment
    alignGenome = "";
    alignRead   = "";
    for (int i=0; i<alignBlocks.size(); i++)
    {
        blockRead   = get<0>(alignBlocks[i]);
        blockGenome = get<1>(alignBlocks[i]);

        alignRead   += blockRead;
        alignGenome += blockGenome;

        for (int j=blockRead.size(); j<blockGenome.size(); j++)
            alignRead += '-';

        for (int j=blockGenome.size(); j<blockRead.size(); j++)
            alignGenome += '-';
    }
}


void GenericBamAlignmentTools::moveInDelLeft(string &alignReadSeq, string &alignGenomeSeq)
{
    // change alignment into blocks
    vector<BamAlignmentBlock> blocks;
    divideAlignmentToBlocks(alignReadSeq, alignGenomeSeq, blocks);

    // new read and genome
    string ar(""),ag("");

    // iterate from first block to last block
    for (vector<BamAlignmentBlock>::iterator iter=blocks.begin(); iter!=blocks.end(); ++iter)
    {
        if (iter->Type()==BLOCK_MATCH || iter->Type()==BLOCK_MISMATCH)
        {
            ar.insert(ar.end(), iter->m_read.begin(), iter->m_read.end());
            ag.insert(ag.end(), iter->m_genome.begin(), iter->m_genome.end());
        }

        if (iter->Type()==BLOCK_INSERT)
        {
            ar.insert(ar.end(), iter->m_read.begin(), iter->m_read.end());
            ag.insert(ag.end(), iter->m_read.length(), Spa);
        }

        if (iter->Type()==BLOCK_DELETE)
        {
            ar.insert(ar.end(), iter->m_genome.length(), Spa);
            ag.insert(ag.end(), iter->m_genome.begin(), iter->m_genome.end());
        }

        if (iter->Type()==BLOCK_OVERCALL)
        {
            ar.insert(ar.end(), iter->m_read.begin(), iter->m_read.end());
            ag.insert(ag.end(), iter->m_read.length()-iter->m_genome.length(), Spa);
            ag.insert(ag.end(), iter->m_genome.begin(), iter->m_genome.end());

        }

        if (iter->Type()==BLOCK_UNDERCALL)
        {
            ar.insert(ar.end(), iter->m_genome.length()-iter->m_read.length(), Spa);
            ar.insert(ar.end(), iter->m_read.begin(), iter->m_read.end());
            ag.insert(ag.end(), iter->m_genome.begin(), iter->m_genome.end());
        }
    }

    // save into alignReadSeq and alignGenomeSeq
    alignReadSeq   = ar;
    alignGenomeSeq = ag;
}


void GenericBamAlignmentTools::moveInDelRight(string &alignReadSeq, string &alignGenomeSeq)
{
    // change alignment into blocks
    vector<BamAlignmentBlock> blocks;
    divideAlignmentToBlocks(alignReadSeq, alignGenomeSeq, blocks);

    // new read and genome
    string ar(""),ag("");

    // iterate from first block to last block
    for (vector<BamAlignmentBlock>::iterator iter=blocks.begin(); iter!=blocks.end(); ++iter)
    {

        if (iter->Type()==BLOCK_MATCH || iter->Type()==BLOCK_MISMATCH)
        {
            ar.insert(ar.end(), iter->m_read.begin(), iter->m_read.end());
            ag.insert(ag.end(), iter->m_genome.begin(), iter->m_genome.end());
        }

        if (iter->Type()==BLOCK_INSERT)
        {
            ar.insert(ar.end(), iter->m_read.begin(), iter->m_read.end());
            ag.insert(ag.end(), iter->m_read.length(), Spa);
        }

        if (iter->Type()==BLOCK_DELETE)
        {
            ar.insert(ar.end(), iter->m_genome.length(), Spa);
            ag.insert(ag.end(), iter->m_genome.begin(), iter->m_genome.end());
        }

        if (iter->Type()==BLOCK_OVERCALL)
        {
            ar.insert(ar.end(), iter->m_read.begin(), iter->m_read.end());
            ag.insert(ag.end(), iter->m_genome.begin(), iter->m_genome.end());
            ag.insert(ag.end(), iter->m_read.length()-iter->m_genome.length(), Spa);
        }

        if (iter->Type()==BLOCK_UNDERCALL)
        {
            ar.insert(ar.end(), iter->m_read.begin(), iter->m_read.end());
            ar.insert(ar.end(), iter->m_genome.length()-iter->m_read.length(), Spa);
            ag.insert(ag.end(), iter->m_genome.begin(), iter->m_genome.end());
        }
    }

    // save into alignReadSeq and alignGenomeSeq
    alignReadSeq   = ar;
    alignGenomeSeq = ag;
}


void GenericBamAlignmentTools::calculateCigar(const string &alignRead, const string &alignGenome, Cigar &cigar)
{
    char Match  = 'M';
    char Insert = 'I';
    char Delete = 'D';

    bool blockChange;
    char blockType;
    int  blockLength;

    for (int i=0; i<alignRead.size(); i++)
    {
        // initialization
        if (i==0)
        {
            if (alignRead[i]=='-')
                blockType = Delete;
            if (alignGenome[i]=='-')
                blockType = Insert;
            if (alignRead[i]!='-' && alignGenome[i]!='-')
                blockType = Match;

            blockLength = 1;

            blockChange = false;

            continue;
        }

        // block change
        if (blockType==Match && (alignRead[i]=='-' || alignGenome[i]=='-'))
        {
            blockChange = true;
        }
        if (blockType==Delete && alignRead[i]!='-')
        {
            blockChange = true;
        }
        if (blockType==Insert && alignGenome[i]!='-')
        {
            blockChange = true;
        }

        // save previous block
        if (blockChange)
        {
            cigar.push_back(CigarOp(blockType, blockLength));

            if (alignRead[i]=='-')
                blockType = Delete;
            if (alignGenome[i]=='-')
                blockType = Insert;
            if (alignRead[i]!='-' && alignGenome[i]!='-')
                blockType = Match;

            blockLength = 0;

            blockChange = false;

        }

        // update block information
        blockLength += 1;
    }

    if (!blockChange)
        cigar.push_back(CigarOp(blockType, blockLength));
}

void GenericBamAlignmentTools::calculateMD(const string &alignRead, const string &alignGenome, BamMD &md)
{
    char Match    = 'M';
    char Mismatch = 'X';
    char Insert   = 'I';
    char Delete   = 'D';

    // the state sequence of alignment
    vector<char> alignState;
    for (int i=0; i<alignRead.length(); ++i)
    {
        if (alignRead[i]!=Spa && alignGenome[i]==Spa)
        {
            alignState.push_back(Insert);
            continue;
        }

        if (alignRead[i]==Spa && alignGenome[i]!=Spa)
        {
            alignState.push_back(Delete);
            continue;
        }

        if (alignRead[i]!=alignGenome[i])
        {
            alignState.push_back(Mismatch);
            continue;
        }

        alignState.push_back(Match);
    }

    // vector of substrings
    vector<tuple<char,string>> blockStrings;
    for (int i=0; i<alignRead.length(); ++i)
    {
        if (i==0)
        {
            blockStrings.push_back(tuple<char,string>(alignState[i], alignGenome.substr(i,1)));
            continue;
        }

        // current state is match
        if (alignState[i]==Match)
        {
            if (alignState[i-1]==alignState[i])
            {
                vector<tuple<char,string>>::reverse_iterator riter = blockStrings.rbegin();
                get<1>(*riter) += alignGenome.substr(i,1);
            }else
            {
                blockStrings.push_back(tuple<char,string>(alignState[i],alignGenome.substr(i,1)));
            }

            continue;
        }

        // current state is mismatch
        if (alignState[i]==Mismatch)
        {
            blockStrings.push_back(tuple<char,string>(alignState[i],alignGenome.substr(i,1)));

            continue;
        }

        // current state is insert
        if (alignState[i]==Insert)
        {
            if (alignState[i-1]==Insert)
            {
                vector<tuple<char,string>>::reverse_iterator riter = blockStrings.rbegin();
                get<1>(*riter) += alignGenome.substr(i,1);
            }else
            {
                blockStrings.push_back(tuple<char,string>(alignState[i],alignGenome.substr(i,1)));
            }

            continue;
        }

        // current state is delete
        if (alignState[i]==Delete)
        {
            if (alignState[i-1]==Delete)
            {
                vector<tuple<char,string>>::reverse_iterator riter = blockStrings.rbegin();
                get<1>(*riter) += alignGenome.substr(i,1);
            }else
            {
                blockStrings.push_back(tuple<char,string>(alignState[i],alignGenome.substr(i,1)));
            }

            continue;
        }
    }

    // remove insert and merge neighbor match
    vector<tuple<char,string>> newBlockStrings;
    for (int i=0; i<blockStrings.size(); ++i)
    {
        if (i==0)
        {
            newBlockStrings.push_back(blockStrings[i]);
            continue;
        }
        // match
        int j = newBlockStrings.size()-1;
        if (isMatch(get<0>(blockStrings[i])))
        {
            if (isMatch(get<0>(newBlockStrings[j])))
            {
                get<1>(newBlockStrings[j]) += get<1>(blockStrings[i]);
            }else
            {
                newBlockStrings.push_back(blockStrings[i]);
            }
            continue;
        }
        // mismatch
        if (isMismatch(get<0>(blockStrings[i])))
        {
            newBlockStrings.push_back(blockStrings[i]);
            continue;
        }
        // delete
        if (isDelete(get<0>(blockStrings[i])))
        {
            newBlockStrings.push_back(blockStrings[i]);
            continue;
        }

    }

    // MD string
    md = "";
    int i=0;
    vector<tuple<char,string>>::iterator iter=newBlockStrings.begin();
    for (; i<newBlockStrings.size(); iter++, i++)
    {
        if (isInsert(get<0>(*iter)))
            continue;

        if (isMatch(get<0>(*iter)))
        {
            md += to_string(get<1>(*iter).length());
            continue;
        }

        if (isMismatch(get<0>(*iter)))
        {
            if (i==0)
            {
                md += to_string(0) + get<1>(*iter);
            }else
            {
                vector<tuple<char,string>>::iterator iter2=prev(iter);

                if (isMatch(get<0>(*iter2)))
                {
                    md += get<1>(*iter);
                }else
                {
                    md += to_string(0) + get<1>(*iter);
                }
            }
            continue;
        }

        if (isDelete(get<0>(*iter)))
        {
            md += "^" + get<1>(*iter);
            continue;
        }
    }
}




void GenericBamAlignmentTools::getBamAlignmentMD(BamAlignment &alignObj, BamMdOpArray &bamMdOps)
{
    BamMD md;
    alignObj.GetTag("MD", md);

    bool blockChange;
    int  blockLength;
    BamMdType blockType;
    string blockSequence;
    stringstream buffer;

    blockChange = false;
    blockLength = 0;
    blockType = BamMdMatch;

    for (int i=0; i<md.length(); i++)
    {
        // block change
        if (md[i]=='^')
            blockChange = true;

        if (blockType==BamMdDelete && isdigit(md[i]))
            blockChange = true;

        if (blockType==BamMdMismatch && isdigit(md[i]))
        {
            if (md[i]=='0')
                continue;
            else
                blockChange = true;
        }

        if (blockType==BamMdMatch && !isdigit(md[i]))
            blockChange = true;

        // save previous block
        if (blockChange)
        {
            // turn off
            blockChange = false;

            // block length and sequence
            blockSequence = string();
            if (blockType==BamMdMatch)
                buffer >> blockLength;
            if (blockType==BamMdMismatch || blockType==BamMdDelete)
            {
                buffer >> blockSequence;
                blockLength = blockSequence.length();
            }

            // save into array
            bamMdOps.push_back(BamMdOp(blockType, blockLength, blockSequence));

            // reset buffer
            buffer.str(string());
            buffer.clear();

            // update block type
            if (md[i]=='^')
            {
                blockType = BamMdDelete;
                continue;
            }

            if (isdigit(md[i]))
                blockType = BamMdMatch;

            if (!isdigit(md[i]))
                blockType = BamMdMismatch;
        }

        buffer << md[i];
    }

    if (!blockChange)
    {
        blockSequence = string();
        // block length
        if (blockType==BamMdMatch)
            buffer >> blockLength;
        if (blockType==BamMdMismatch || blockType==BamMdDelete)
        {
            buffer >> blockSequence;
            blockLength = blockSequence.length();
        }

        // save into array
        bamMdOps.push_back(BamMdOp(blockType, blockLength, blockSequence));
    }
}

void GenericBamAlignmentTools::getBamAlignmentCigar(BamAlignment &alignObj, Cigar &cigar)
{
    cigar.assign(alignObj.CigarData.begin(), alignObj.CigarData.end());
}

void GenericBamAlignmentTools::getBamAlignmentGenome(BamAlignment &alignObj, string &genome)
{
    BamMD md;

    // make sure data has MD tag
    if (!alignObj.HasTag("MD"))
    {
        cerr << "GenericBamAlignmentTools ERROR: "
             << alignObj.Name << " "
             << "has no MD tag..." << " "
             << endl
             << "Use 'samtools calmd' to fix it"
             << endl;
        exit(EXIT_FAILURE);
    }

    // get MD tag
    alignObj.GetTag("MD", md);

    // parse MD string to array
    BamMdOpArray bamMdOps;
    getBamAlignmentMD(alignObj, bamMdOps);

    // genome
    genome = string();
    for (BamMdOpArray::iterator iter=bamMdOps.begin(); iter!=bamMdOps.end(); ++iter)
    {
        genome.insert(genome.end(), (*iter).Nucleotides.begin(), (*iter).Nucleotides.end());
    }

    // iterate through Cigar
    int readPointer = 0;
    int genomePointer = 0;

    for (Cigar::iterator iter=alignObj.CigarData.begin(); iter!=alignObj.CigarData.end(); ++iter)
    {
        if (isSoftClip(iter->Type)) readPointer += iter->Length;

        if (isMatch(iter->Type) || isMismatch(iter->Type))
        {
            for (int i=0; i<iter->Length; ++i, ++readPointer, ++genomePointer)
            {
                if (genome[genomePointer]==Amb)
                    genome[genomePointer] = alignObj.QueryBases[readPointer];
            }
        }

        if (isInsert(iter->Type)) readPointer += iter->Length;

        if (isDelete(iter->Type)) genomePointer += iter->Length;
    }
}




void GenericBamAlignmentTools::getAlignmentBlocks(const BamAlignment &alignObj, VectorInteger &blockReadPositions, VectorInteger &blockReadSizes, VectorInteger &blockGenomePositions, VectorInteger &blockGenomeSizes, VectorCigar &blockCigars, vector<RunInDelShift> &blockRuns)
{
    // iterate throught the cigar of the alignment
    int blockReadStart   = 0;
    int blockGenomeStart = alignObj.Position;
    int blockReadSize    = 0;
    int blockGenomeSize  = 0;
    Cigar blockCigar;
    string nonClip = "non-clip";
    string clip    = "clip";
    string blockType   = nonClip;
    bool blockChange   = false;
    for (int i=0; i<alignObj.CigarData.size(); i++){

        // check cigar operator
        if ( isClip(alignObj.CigarData[i].Type)     &&
             isMatch(alignObj.CigarData[i].Type)    &&
             isMismatch(alignObj.CigarData[i].Type) &&
             isGap(alignObj.CigarData[i].Type)
           )
        {
            stringstream CigarSequence;
            for (int j=0; j<alignObj.CigarData.size(); j++)
                CigarSequence << alignObj.CigarData[j].Length
                              << alignObj.CigarData[j].Type;
            cerr << "GenericBamAlignmentTools ERROR:" << " "
                 << "unknown Cigar operator " << alignObj.CigarData[i].Type << "\n"
                 << "read name is " << alignObj.Name << ", "
                 << "cigar is " << CigarSequence.str()
                 << endl;
            exit(0);
        }

        // the signal of block change
        if (isClip(alignObj.CigarData[i].Type) && blockType==nonClip)
        {
            blockChange = true;
        }

        // the signal of block change
        if (!isClip(alignObj.CigarData[i].Type) && blockType==clip)
        {
            blockChange = true;
        }

        //
        if (blockChange)
        {
            // turn off
            blockChange = false;

            if ((blockReadSize+blockGenomeSize)>0){
                // save the block
                blockReadPositions.push_back(blockReadStart);
                blockReadSizes.push_back(blockReadSize);
                blockGenomePositions.push_back(blockGenomeStart);
                blockGenomeSizes.push_back(blockGenomeSize);
                blockCigars.push_back(blockCigar);

                if (blockType==clip)
                    blockRuns.push_back(SkipInDelShift);
                if (blockType==nonClip)
                    blockRuns.push_back(DoInDelShift);
            }

            // update the variable
            if (blockType==clip)
            {
                blockType = nonClip;
            }
            else if (blockType==nonClip)
            {
                blockType = clip;
            }

            blockReadStart   += blockReadSize;
            blockGenomeStart += blockGenomeSize;

            blockReadSize   = 0;
            blockGenomeSize = 0;

            blockCigar.clear();
        }

        // update the variable
        if ( isSoftClip(alignObj.CigarData[i].Type)     ||
             isMatch(alignObj.CigarData[i].Type)    ||
             isMismatch(alignObj.CigarData[i].Type) ||
             isInsert(alignObj.CigarData[i].Type)
           )
        {
            blockReadSize += alignObj.CigarData[i].Length;
        }

        if ( isMatch(alignObj.CigarData[i].Type)    ||
             isMismatch(alignObj.CigarData[i].Type) ||
             isDelete(alignObj.CigarData[i].Type)
           )
        {
            blockGenomeSize += alignObj.CigarData[i].Length;
        }

        blockCigar.push_back(alignObj.CigarData[i]);
    }

    if (!blockChange){
        blockReadPositions.push_back(blockReadStart);
        blockReadSizes.push_back(blockReadSize);
        blockGenomePositions.push_back(blockGenomeStart);
        blockGenomeSizes.push_back(blockGenomeSize);
        blockCigars.push_back(blockCigar);

        if (blockType==clip)
            blockRuns.push_back(SkipInDelShift);
        if (blockType==nonClip)
            blockRuns.push_back(DoInDelShift);
    }
}

void GenericBamAlignmentTools::getAlignmentSequences(const string &read, const string &genome, Cigar &cigar, string &alignRead, string &alignGenome)
{
    alignRead   = "";
    alignGenome = "";

    int readPointer   = 0;
    int genomePointer = 0;

    for (int i=0; i<cigar.size(); i++){

        int shiftSize = cigar[i].Length;

        if (isMatch(cigar[i].Type) || isMismatch(cigar[i].Type))
        {
            for (int j=0; j<shiftSize; j++)
            {
                alignRead   += read[readPointer+j];
                alignGenome += genome[genomePointer+j];
            }

            // update pointers
            readPointer   += shiftSize;
            genomePointer += shiftSize;
        }

        if (isDelete(cigar[i].Type))
        {
            for (int j=0; j<shiftSize; j++)
            {
                alignRead   += '-';
                alignGenome += genome[genomePointer+j];
            }
            // update pointer
            genomePointer += shiftSize;
        }

        if (isInsert(cigar[i].Type))
        {
            for (int j=0; j<shiftSize; j++)
            {
                alignRead   += read[readPointer+j];
                alignGenome += '-';
            }
            // update pointer
            readPointer += shiftSize;
        }

        if (isSoftClip(cigar[i].Type))
        {
            readPointer += shiftSize;
        }
    }
}

void GenericBamAlignmentTools::getAlignmentSequences(BamAlignment &alignObj, string &alignRead, string &alignGenome)
{
    alignRead   = "";
    alignGenome = "";

    string read = alignObj.QueryBases;

    string genome;
    getBamAlignmentGenome(alignObj, genome);

    int readPointer   = 0;
    int genomePointer = 0;

    for (int i=0; i<alignObj.CigarData.size(); i++)
    {
        int shiftSize = alignObj.CigarData[i].Length;

        if ( isMatch(alignObj.CigarData[i].Type)    ||
             isMismatch(alignObj.CigarData[i].Type)
           )
        {
            for (int j=0; j<shiftSize; j++)
            {
                alignRead   += read[readPointer+j];
                alignGenome += genome[genomePointer+j];
            }
            // update pointers
            readPointer   += shiftSize;
            genomePointer += shiftSize;
        }

        if (isDelete(alignObj.CigarData[i].Type))
        {
            for (int j=0; j<shiftSize; j++)
            {
                alignRead   += '-';
                alignGenome += genome[genomePointer+j];
            }
            // update pointer
            genomePointer += shiftSize;
        }

        if (isInsert(alignObj.CigarData[i].Type))
        {
            for (int j=0; j<shiftSize; j++)
            {
                alignRead   += read[readPointer+j];
                alignGenome += '-';
            }
            // update pointer
            readPointer += shiftSize;
        }

        if (isSoftClip(alignObj.CigarData[i].Type))
        {
            // update pointer
            readPointer += shiftSize;
        }
    }
}



void GenericBamAlignmentTools::divideAlignmentToBlocks(const string &alignRead, const string &alignGenome, vector<BamAlignmentBlock> &blocks)
{
    int t_MATCH = 0;
    int t_MISMATCH = 1;
    int t_INSERT = 2;
    int t_DELETE = 3;

    BlockState blockGenomeState;
    BlockState blockReadState;
    bool blockChange;
    int blockLength;
    int blockReadLength;
    int blockGenomeLength;
    int readPointer;
    int genomePointer;
    stringstream alignReadBuf;
    stringstream alignGenomeBuf;
    string::const_iterator alignReadIter;
    string::const_iterator alignGenomeIter;
    int alignLength;

    // alignment state
    vector<int> alignState;
    for (int i=0; i<alignRead.size(); ++i)
    {
        if (alignRead[i]!=Spa && alignGenome[i]!=Spa)
        {
            if (alignRead[i]!=alignGenome[i] && alignRead[i]!=Amb && alignGenome[i]!=Amb)
            {
                alignState.push_back(t_MISMATCH);
            }else
            {
                alignState.push_back(t_MATCH);
            }
        }else if (alignRead[i]==Spa)
        {
            alignState.push_back(t_DELETE);
        }else if (alignGenome[i]==Spa)
        {
            alignState.push_back(t_INSERT);
        }
    }

    // iterate from alignment begin to alignment end
    blockChange       = false;
    blockGenomeState  = Nil;
    blockReadState    = Nil;
    blockLength       = 0;
    blockReadLength   = 0;
    blockGenomeLength = 0;

    readPointer       = 0;
    genomePointer     = 0;

    alignReadIter     = alignRead.begin();
    alignGenomeIter   = alignGenome.begin();

    alignLength = alignRead.length();
    for (int i=0; i<alignLength; ++i, ++alignReadIter, ++alignGenomeIter)
    {

        // mismatch
        if (alignState[i]==t_MISMATCH)
            blockChange = true;

        // delete
        if (alignState[i]==t_DELETE && (*alignGenomeIter)!=blockGenomeState)
            blockChange = true;

        // insert
        if (alignState[i]==t_INSERT && (*alignReadIter)!=blockReadState)
            blockChange = true;

        // match
        if (alignState[i]==t_MATCH)
        {
            if ((*alignReadIter)!=blockReadState && (*alignReadIter)!=Amb)
                blockChange = true;

            if ((*alignGenomeIter)!=blockGenomeState && (*alignReadIter)!=Amb)
                blockChange = true;
        }

        // save previous block
        if (blockChange)
        {
            // save into array
            if (blockLength>0)
            {
                if (blockGenomeLength==0)
                {
                    blocks.push_back(BamAlignmentBlock(alignReadBuf.str(), alignGenomeBuf.str(),
                                                       readPointer, genomePointer-1));
                }
                else if (blockReadLength==0)
                {
                    blocks.push_back(BamAlignmentBlock(alignReadBuf.str(), alignGenomeBuf.str(),
                                                       readPointer-1, genomePointer));
                }
                else
                {
                    blocks.push_back(BamAlignmentBlock(alignReadBuf.str(), alignGenomeBuf.str(),
                                                       readPointer, genomePointer));
                }
            }

            // reset block state
            blockChange = false;

            blockGenomeState = (*alignGenomeIter);
            if ((*alignGenomeIter)==Spa || (*alignGenomeIter)==Amb)
                blockGenomeState = (*alignReadIter);

            blockReadState   = (*alignReadIter);
            if ((*alignReadIter)==Spa || (*alignReadIter)==Amb)
                blockReadState = (*alignGenomeIter);

            readPointer   += blockReadLength;
            genomePointer += blockGenomeLength;

            blockLength       = 0;
            blockReadLength   = 0;
            blockGenomeLength = 0;

            alignReadBuf.str(string());
            alignReadBuf.clear();

            alignGenomeBuf.str(string());
            alignGenomeBuf.clear();
        }


        // update variables
        blockLength++;

        if ((*alignReadIter)!=Spa)
            blockReadLength ++;

        if ((*alignGenomeIter)!=Spa)
            blockGenomeLength ++;

        alignReadBuf << (*alignReadIter);
        alignGenomeBuf<< (*alignGenomeIter);
    }

    if (!blockChange)
    {
        // save into array
        if (blockLength>0)
        {
            if (blockGenomeLength==0)
            {
                blocks.push_back(BamAlignmentBlock(alignReadBuf.str(), alignGenomeBuf.str(),
                                                   readPointer, genomePointer-1));
            }
            else if (blockReadLength==0)
            {
                blocks.push_back(BamAlignmentBlock(alignReadBuf.str(), alignGenomeBuf.str(),
                                                   readPointer-1, genomePointer));
            }
            else
            {
                blocks.push_back(BamAlignmentBlock(alignReadBuf.str(), alignGenomeBuf.str(),
                                                   readPointer, genomePointer));
            }
        }
    }
}


void GenericBamAlignmentTools::printBamAlignmentCigar(BamAlignment &alignObj)
{
    stringstream buffer;
    for (int i=0; i<alignObj.CigarData.size(); i++)
    {
        buffer << alignObj.CigarData[i].Length;
        buffer << alignObj.CigarData[i].Type;
    }
    cout << alignObj.Name << ", "
         << alignObj.RefID+1 << ":"
         << alignObj.Position+1 << "-"
         << alignObj.GetEndPosition()+1 << ", "
         << "Cigar: "
         << buffer.str() << endl;
}

void GenericBamAlignmentTools::printBamAlignmentMD(BamAlignment &alignObj)
{
    BamMD md;
    if (alignObj.HasTag("MD"))
        alignObj.GetTag("MD", md);

    cout << alignObj.Name << ", "
         << alignObj.RefID+1 << ":"
         << alignObj.Position+1 << "-"
         << alignObj.GetEndPosition()+1 << ", "
         << "MD: "
         << md << endl;
}

bool GenericBamAlignmentTools::isSoftClip(char op)
{
    return (op=='S');
}

bool GenericBamAlignmentTools::isHardClip(char op)
{
    return (op=='H');
}

bool GenericBamAlignmentTools::isClip(char op)
{
    return (isSoftClip(op) || isHardClip(op));
}

bool GenericBamAlignmentTools::isMatch(char op)
{
    return ((op=='M') || (op=='='));
}

bool GenericBamAlignmentTools::isMismatch(char op)
{
    return (op=='X');
}

bool GenericBamAlignmentTools::isInsert(char op)
{
    return (op=='I');
}

bool GenericBamAlignmentTools::isDelete(char op)
{
    return (op=='D');
}

bool GenericBamAlignmentTools::isGap(char op)
{
    return (isDelete(op) || isInsert(op));
}

bool GenericBamAlignmentTools::isNotGap(char op)
{
    return (isMatch(op) || isMismatch(op));
}




int GenericBamAlignmentTools::numBamAlignmentMatches(BamAlignment &alignObj)
{
    int number = 0;
    bool toCount = false;

    BamMD md;
    alignObj.GetTag("MD", md);

    stringstream buffer;
    for (BamMD::iterator iter=md.begin(); iter!=md.end(); ++iter)
    {

        if (isdigit(*iter))
        {
            toCount = true;
            buffer << (*iter);
        }
        else
        {
            toCount = false;

            int length;
            buffer >> length;

            number += length;

            buffer.str(string());
            buffer.clear();
        }
    }

    if (toCount)
    {
        int length;
        buffer >> length;

        number += length;

        buffer.str(string());
        buffer.clear();
    }


    return number;
}


int GenericBamAlignmentTools::numBamAlignmentMismatches(BamAlignment &alignObj)
{
    int number = 0;

    BamMD md;
    alignObj.GetTag("MD", md);

    for (BamMD::iterator iter=md.begin()+1; iter!=md.end(); ++iter)
    {

        BamMD::iterator iterPrev = prev(iter, 1);
        if (isalpha((*iter)) && isdigit((*iterPrev)))
            number ++;
    }

    // remove ambiguous bases
    int readPointer = 0;

    for (Cigar::iterator iter=alignObj.CigarData.begin(); iter!=alignObj.CigarData.end(); ++iter)
    {
        if (isSoftClip(iter->Type) || isInsert(iter->Type))
            readPointer += iter->Length;

        if (isMatch(iter->Type) || isMismatch(iter->Type))
        {
            for (int i=0; i<iter->Length; ++i, ++readPointer)
            {
                if (alignObj.QueryBases[readPointer]==Amb)
                    number --;
            }
        }
    }

    return number;
}

void GenericBamAlignmentTools::getBamAlignmentMismatches(BamAlignment &alignObj,
                                                         vector<long> &genomePositions,
                                                         vector<char> &referenceAlleles,
                                                         vector<char> &readAlleles,
                                                         vector<double> &readAlleleQualities)
{
    int readPointer;
    int genomePointer;

    // map from genome pointer to read pointer
    vector<int> g2r;

    // iterate through Cigar sequence
    readPointer   = 0;
    for (Cigar::iterator iter=alignObj.CigarData.begin(); iter!=alignObj.CigarData.end(); ++iter)
    {
        // soft-clip and insert
        if (isSoftClip(iter->Type) || isInsert(iter->Type))
            readPointer += iter->Length;

        // match and mismatch
        if (isMatch(iter->Type) || isMismatch(iter->Type))
        {
            for (int i=0; i<iter->Length; ++i, ++readPointer)
            {
                g2r.push_back(readPointer);
            }
        }

        // delete
        if (isDelete((iter->Type)))
        {
            for (int i=0; i<iter->Length; ++i)
            {
                g2r.push_back(readPointer);
            }
        }
    }

    // retrieve MD information
    BamMdOpArray MdInfo;
    getBamAlignmentMD(alignObj, MdInfo);

    // iterate through MD sequence
    genomePointer = 0;
    for (BamMdOpArray::iterator iter=MdInfo.begin(); iter!=MdInfo.end(); ++iter)
    {
        if (iter->Type==BamMdMatch || iter->Type==BamMdDelete)
            genomePointer += iter->Length;

        if (iter->Type==BamMdMismatch)
        {
            for (int i=0; i<iter->Nucleotides.size(); ++i, ++genomePointer)
            {
                readPointer = g2r[genomePointer];

                // skip if base N
                if (iter->Nucleotides[i]==Amb) continue;
                if (alignObj.QueryBases[readPointer]==Amb) continue;

                genomePositions.push_back(alignObj.Position+genomePointer);

                if (&referenceAlleles!=0)
                    referenceAlleles.push_back(iter->Nucleotides[i]);

                if (&readAlleles!=0)
                    readAlleles.push_back(alignObj.QueryBases[readPointer]);

                if (&readAlleleQualities!=0)
                    readAlleleQualities.push_back(int(alignObj.Qualities[readPointer])-33);
            }
        }
    }


}



int GenericBamAlignmentTools::numBamAlignmentDeletes(BamAlignment &alignObj)
{
    int number = 0;

    BamMD md;
    alignObj.GetTag("MD", md);

    for (BamMD::iterator iter=md.begin(); iter!=md.end(); iter++)
    {
        if ((*iter)=='^')
            number ++;
    }

    return number;
}

void GenericBamAlignmentTools::getBamAlignmentDeletes(BamAlignment &alignObj,
                                                      vector<long> &genomePositions,
                                                      vector<string> &deletedSequences)
{
    // get MD information
    BamMdOpArray MdInfo;
    getBamAlignmentMD(alignObj, MdInfo);

    // iterate through MD sequence
    int genomePointer = 0;
    for (BamMdOpArray::iterator iter=MdInfo.begin(); iter!=MdInfo.end(); ++iter)
    {

        if (iter->Type==BamMdDelete)
        {
            genomePositions.push_back(alignObj.Position + genomePointer);
            deletedSequences.push_back(string(iter->Nucleotides.begin(), iter->Nucleotides.end()));
        }

        genomePointer += iter->Length;
    }
}


int GenericBamAlignmentTools::numBamAlignmentInserts(BamAlignment &alignObj)
{
    int number = 0;

    for (Cigar::iterator iter=alignObj.CigarData.begin(); iter!=alignObj.CigarData.end(); ++iter)
    {
        if (iter->Type=='I')
            number ++;
    }

    return number;
}

void GenericBamAlignmentTools::getBamAlignmentInserts(BamAlignment &alignObj, vector<long> &genomePositions, vector<string> &insertedSequences, vector<vector<double> > &insertedSequenceQualities)
{
    int readPointer;
    int genomePointer;

    // iterate through Cigar sequence
    readPointer   = 0;
    genomePointer = 0;
    for (Cigar::iterator iter=alignObj.CigarData.begin(); iter!=alignObj.CigarData.end(); ++iter)
    {
        if (isSoftClip(iter->Type))
            readPointer += iter->Length;

        if (isMatch(iter->Type) || isMismatch(iter->Type))
        {
            readPointer += iter->Length;
            genomePointer += iter->Length;
        }

        if (isDelete(iter->Type))
            genomePointer += iter->Length;

        if (isInsert(iter->Type))
        {
            string subseq;
            vector<double> subseqQuality;


            for (int i=0; i<iter->Length; ++i, ++readPointer)
            {
                subseq += alignObj.QueryBases[readPointer];
                subseqQuality.push_back(int(alignObj.Qualities[readPointer])-33);
            }

            genomePositions.push_back(genomePointer);
            insertedSequences.push_back(subseq);
            insertedSequenceQualities.push_back(subseqQuality);
        }
    }
}


int GenericBamAlignmentTools::getBamAlignmentReadLength(BamAlignment &alignObj)
{
    int len = 0;

    for (Cigar::iterator iter=alignObj.CigarData.begin(); iter!=alignObj.CigarData.end(); ++iter)
    {
        if (isMatch(iter->Type) || isMismatch(iter->Type) || isInsert(iter->Type))
            len += iter->Length;
    }

    return len;
}


int GenericBamAlignmentTools::maxGapLength(BamAlignment &alignObj)
{
    int len = 0;
    for (Cigar::iterator iter=alignObj.CigarData.begin(); iter!=alignObj.CigarData.end(); ++iter)
    {
        if (isGap(iter->Type))
        {
            if (len<iter->Length)
                len = iter->Length;
        }
    }

    return len;
}

int GenericBamAlignmentTools::maxGapLength(string &alignRead, string &alignGenome)
{
    Cigar cigar;
    calculateCigar(alignRead, alignGenome, cigar);

    int len = 0;
    for (Cigar::iterator iter=cigar.begin(); iter!=cigar.end(); ++iter)
    {
        if (isGap(iter->Type))
        {
            if (len<iter->Length)
                len = iter->Length;
        }
    }

    return len;
}

int GenericBamAlignmentTools::maxGapLength(Cigar &cigar)
{
    int len = 0;
    for (Cigar::iterator iter=cigar.begin(); iter!=cigar.end(); ++iter)
    {
        if (isGap(iter->Type))
        {
            if (len<iter->Length)
                len = iter->Length;
        }
    }

    return len;
}

bool GenericBamAlignmentTools::goodAlignment(BamAlignment &alignObj, bool keepDuplicate)
{
    bool good = true;

    if (!alignObj.IsMapped()
            ||
            alignObj.IsFailedQC()
            ||
            (!keepDuplicate && alignObj.IsDuplicate())
            ||
            !alignObj.IsPrimaryAlignment()
            )
    {
        good = false;
    }

    if (alignObj.IsPaired() && !alignObj.IsProperPair())
    {
        good = false;
    }

    return good;
}



void GenericBamAlignmentTools::getLocalAlignment(BamAlignment &alignObj, int position, int size,
                                                 string &localRead, string &localGenome,
                                                 Cigar &cigar, BamMD &md,
                                                 int &numMismatch, int &numInDel)
{
    // get the whole alignment
    string alignRead, alignGenome;
    getAlignmentSequences(alignObj, alignRead, alignGenome);

    // move insert/delete to right
    moveInDelRight(alignRead, alignGenome);

    // the shift from position to the begin of alignment
    int genomeShift = position-alignObj.Position;

    // local alignment
    string localAlignRead(""), localAlignGenome("");

    // pointer on aligned read and aligned genome
    int readPointer = -1, genomePointer = -1;

    // iterate from first align position to last align position
    for (int i=0; i<alignRead.length(); ++i)
    {
        if (alignRead[i]!=Spa)
            readPointer += 1;

        if (alignGenome[i]!=Spa)
            genomePointer += 1;

        if (genomePointer>=genomeShift && genomePointer<genomeShift+size)
        {
            localAlignRead   += alignRead[i];
            localAlignGenome += alignGenome[i];
        }
    }

    // calculate Cigar
    if (&cigar!=0)
    {
        calculateCigar(localAlignRead, localAlignGenome, cigar);
    }

    // calculate MD
    if (&md!=0)
    {
        calculateMD(localAlignRead, localAlignGenome, md);
    }

    // compute the number of mismatches
    if (&numMismatch!=0)
    {
        numMismatch = 0;
        for (int i=0; i<localAlignRead.size(); ++i)
        {
            if (localAlignRead[i]==Amb || localAlignRead[i]==Spa)
                continue;

            if (localAlignGenome[i]==Amb || localAlignGenome[i]==Spa)
                continue;

            if (localAlignRead[i]!=localAlignGenome[i])
                numMismatch += 1;
        }
    }

    // compute the number of indels
    if (&numInDel!=0)
    {
        numInDel = 0;
        for (int i=0; i<localAlignRead.size(); ++i)
        {
            if (localAlignRead[i]==Spa)
                numInDel += 1;

            if (localAlignGenome[i]==Spa)
                numInDel += 1;
        }
    }

    // remove all spaces
    localRead   = "";
    localGenome = "";
    for (int i=0; i<localAlignRead.size(); ++i)
    {
        if (localAlignRead[i]!=Spa)
            localRead += localAlignRead[i];

        if (localAlignGenome[i]!=Spa)
            localGenome += localAlignGenome[i];
    }
}

string GenericBamAlignmentTools::getBamAlignmentID(BamAlignment &alignObj)
{
    return (alignObj.Name+":"+to_string(alignObj.AlignmentFlag));
}


int GenericBamAlignmentTools::numMatch(string &alignRead, string &alignGenome)
{
    int nmt = 0;

    for (int i=0; i<alignRead.length(); ++i)
    {
        if (alignRead[i]==Spa || alignGenome[i]==Spa)
            continue;

        if (alignRead[i]==alignGenome[i] || alignRead[i]==Amb || alignGenome[i]==Amb)
        {
            nmt++;
            continue;
        }
    }

    return nmt;
}

int GenericBamAlignmentTools::numMismatch(string &alignRead, string &alignGenome)
{
    int nms = 0;

    for (int i=0; i<alignRead.length(); ++i)
    {
        if (alignRead[i]==Spa || alignGenome[i]==Spa)
            continue;

        if (alignRead[i]!=alignGenome[i] && alignRead[i]!=Amb && alignGenome[i]!=Amb)
        {
            nms++;
            continue;
        }
    }

    return nms;
}


void GenericBamAlignmentTools::changeLocalAlignment(BamAlignment &alignObj, int position, int size, string &alnLocalRead, string &alnLocalGenome, BamAlignment& newAlignObj)
{
    // the aligned read and genome of the given alignment
    string alnRead, alnGenome;
    getAlignmentSequences(alignObj, alnRead, alnGenome);
    moveInDelRight(alnRead, alnGenome);

    // the start position of the given alignment
    int alnStartPos = alignObj.Position;

    // the shift between start position of interesting and alignment start position
    int genomeShift = position-alnStartPos;

    // the start and end pointer of the local part that will be replaced with new one
    int i0=-1, i1=-1;

    // the genome pointer
    int genomePointer=-1;
    for (int k=0; k<alnRead.length(); k++)
    {
        // current genome position
        if (alnGenome[k]!=Spa)
        {
            genomePointer++;
        }

        // check whether it is the start position of local part
        if (genomePointer>=genomeShift && i0<0)
        {
            i0 = k;
        }

        // check whether it is the end position of local part
        if (genomePointer>=genomeShift && genomePointer<genomeShift+size)
        {
            i1 = k;
        }
    }

    // old sequence in genome
    string oldLocalGenome;
    oldLocalGenome=alnGenome.substr(i0, i1-i0+1);

    // replacing
    alnRead.replace(i0, i1-i0+1, alnLocalRead);
    alnGenome.replace(i0, i1-i0+1, alnLocalGenome);

    // move indel to right
    moveInDelRight(alnRead, alnGenome);

    // compute new cigar
    Cigar newCigar;
    calculateCigar(alnRead, alnGenome, newCigar);
    // compute new MD
    string newMD;
    calculateMD(alnRead, alnGenome, newMD);

    // old Cigar
    Cigar oldCigar = alignObj.CigarData;

    // soft clip at begin
    Cigar::iterator iter = oldCigar.begin();
    if (isSoftClip(iter->Type))
        newCigar.insert(newCigar.begin(), *iter);
    // soft clip at end
    iter = prev(oldCigar.end());
    if (isSoftClip(iter->Type))
        newCigar.push_back(*iter);

    // new alignment object
    newAlignObj = BamAlignment(alignObj);
    newAlignObj.CigarData = newCigar;
    newAlignObj.EditTag("MD","Z",newMD);
    if (genomeShift<0)
    {
        int l0=0;
        for (int i=0; i<oldLocalGenome.length(); ++i)
        {
            if (oldLocalGenome[i]!=Spa)
                l0++;
        }
        int l1=0;
        for (int i=0; i<alnLocalGenome.length(); ++i)
        {
            if (alnLocalGenome[i]!=Spa)
                l1++;
        }

        newAlignObj.Position = alignObj.Position+l0-l1;
    }
}


bool GenericBamAlignmentTools::validMapQuality(BamAlignment &alignObj, double minMapQuality)
{
    return (alignObj.MapQuality>=minMapQuality);
}

bool GenericBamAlignmentTools::validReadLength(BamAlignment &alignObj, double minReadLength)
{
    int len = GenericBamAlignmentTools::getBamAlignmentReadLength(alignObj);
    return (len>=minReadLength);
}

bool GenericBamAlignmentTools::validReadIdentity(BamAlignment &alignObj, double maxMismatchFrac)
{
    int nX  = numBamAlignmentMismatches(alignObj);
    int len = getBamAlignmentReadLength(alignObj);
    return (nX<len*maxMismatchFrac);
}
