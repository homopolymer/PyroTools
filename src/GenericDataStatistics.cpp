#include "GenericDataStatistics.h"
#include "GenericBamAlignmentTools.h"
using namespace GenericSequenceTools;

#include <iostream>
#include <iterator>
using namespace std;

GenericDataStatistics::GenericDataStatistics(GenericReadBins& bins)
    : m_bins(bins)
{
    //-------------------------------------------------------
    // error rate (exclude insertion/deletion)

    // mismatch counts at each bin
    m_binMismatchCount = VectorDouble(m_bins.m_binNum, 0);
    // total counts including match and mismatch at each bin
    m_binCount = VectorDouble(m_bins.m_binNum, 0);

    //-------------------------------------------------------
    // homopolymer gap

    // the number of homopolymers that are with deletions
    m_homopolymerDelete = VectorDouble(numHomopolymerSize, 0);
    // the number of homopolymers that are with insertions
    m_homopolymerInsert = VectorDouble(numHomopolymerSize, 0);
    // the number of homopolymers that are in data
    m_homopolymerCount  = VectorDouble(numHomopolymerSize, 0);

    //--------------------------------------------------------
    // deletion in data

    // the number of homopolymer-in-sites that are deleted
    m_homopolymerSiteDelete = Matrix2Double(numHomopolymerSize, VectorDouble(numHomopolymerSize, 0));
    // the number of homopolymer-in-sites that are in data
    m_homopolymerSiteCount  = Matrix2Double(numHomopolymerSize, VectorDouble(numHomopolymerSize, 0));

    //--------------------------------------------------------
    // insertion in data

    // the number of homopolymers that are adjacent to insertions,
    // and have the same nucleotide
    // homopolymer is in read
    m_homopolymerNextToInsert = VectorDouble(numHomopolymerSize, 0);

    // the number of homopolymers that are closed to insertions
    // at left, but separated by an interseptal homopolymer
    // homopolymer is in reference
    m_homopolymerLeftCloseToInsert = VectorDouble(numHomopolymerSize, 0);
    // the number of homopolymers that are closed to insertions
    // at right, but separated by an interseptal homopolymer
    // homopolymer is in reference
    m_homopolymerRightCloseToInsert = VectorDouble(numHomopolymerSize, 0);
}

void GenericDataStatistics::update(const string &alignRead, const string &alignGenome)
{
    // match and mismatch
    updateMatchMismatch(alignRead, alignGenome);

    // homopolymer gap
    updateHomopolymerGap(alignRead, alignGenome);

    // delete
    updateDelete(alignRead, alignGenome);

    // insert
    updateInsert(alignRead, alignGenome);

}


void GenericDataStatistics::updateMatchMismatch(const string &alignRead, const string &alignGenome)
{
    string::const_iterator alignReadIter;
    string::const_iterator alignGenomeIter;

    int alignLength = alignRead.length();

    // last position on read
    int readLastPosition = 0;
    for (alignReadIter=alignRead.begin(); alignReadIter!=alignRead.end(); ++alignReadIter)
    {
        if ((*alignReadIter)!=Spa)
            readLastPosition++;
    }
    readLastPosition--;

    // iterate from alignment begin to alignment end
    alignReadIter   = alignRead.begin();
    alignGenomeIter = alignGenome.begin();

    for (int i=0; i<alignLength; ++i, ++alignReadIter, ++alignGenomeIter)
    {
        if ((*alignReadIter)==Spa || (*alignGenomeIter)==Spa)
            continue;

        int x;
        if (2*i>readLastPosition)
            x = readLastPosition-i;
        else
            x = i;

        int binIndex = m_bins.binIndex(x);

        if ((*alignReadIter)!=(*alignGenomeIter))
            m_binMismatchCount[binIndex] += 1;

        m_binCount[binIndex] += 1;
    }
}


void GenericDataStatistics::updateHomopolymerGap(const string &alignRead, const string &alignGenome)
{
    vector<BamAlignmentBlock> alignBlocks;
    GenericBamAlignmentTools::divideAlignmentToBlocks(alignRead, alignGenome, alignBlocks);

    // iterator from first block to last block
    vector<BamAlignmentBlock>::iterator iter;
    for (iter=alignBlocks.begin(); iter!=alignBlocks.end(); ++iter)
    {
        if (iter->Type()==BLOCK_INSERT || iter->Type()==BLOCK_DELETE)
            continue;

        int hlen = iter->genomeLength();

        if (hlen<minHomopolymerSize || hlen>maxHomopolymerSize)
            continue;

        if (iter->Type()==BLOCK_OVERCALL)
            m_homopolymerInsert[hlen] += 1;

        if (iter->Type()==BLOCK_UNDERCALL)
            m_homopolymerDelete[hlen] += 1;

        m_homopolymerCount[hlen] += 1;

    }
}


void GenericDataStatistics::updateDelete(const string &alignRead, const string &alignGenome)
{
    vector<BamAlignmentBlock> alignBlocks;
    GenericBamAlignmentTools::divideAlignmentToBlocks(alignRead, alignGenome, alignBlocks);

    // iterator from first block to last block
    vector<BamAlignmentBlock>::iterator iter;
    for (iter=alignBlocks.begin(); iter!=alignBlocks.end(); ++iter)
    {
        if (iter->Type()==BLOCK_INSERT || iter->Type()==BLOCK_DELETE)
            continue;

        int hlen = iter->genomeLength();

        if (hlen<minHomopolymerSize || hlen>maxHomopolymerSize)
            continue;

        for (int j=0; j<hlen; j++)
        {
            if (j>=iter->readLength())
                m_homopolymerSiteDelete[hlen][j] += 1;

            m_homopolymerSiteCount[hlen][j] += 1;
        }
    }

}


void GenericDataStatistics::updateInsert(const string &alignRead, const string &alignGenome)
{
    vector<BamAlignmentBlock> alignBlocks;
    GenericBamAlignmentTools::divideAlignmentToBlocks(alignRead, alignGenome, alignBlocks);

    // iterator from first block to last block
    vector<BamAlignmentBlock>::iterator iter;
    vector<BamAlignmentBlock>::iterator iterPrev;
    vector<BamAlignmentBlock>::iterator iterNext;
    int i=0;
    for (iter=alignBlocks.begin(); iter!=alignBlocks.end(); ++iter, ++i)
    {

        if (iter->Type()!=BLOCK_INSERT && iter->Type()!=BLOCK_OVERCALL)
            continue;

        int hlen = iter->genomeLength();


        if (iter->Type()==BLOCK_OVERCALL)
        {
            for (int j=hlen; j<iter->readLength(); j++)
            {
                if (j>=minHomopolymerSize && j<=maxHomopolymerSize)
                    m_homopolymerNextToInsert[j] += 1;
            }
        }


        if (iter->Type()==BLOCK_INSERT)
        {
            if (distance(alignBlocks.begin(), iter)>=2)
            {
                iterPrev = prev(iter, 2);

                if ( iterPrev->Type()==BLOCK_MATCH ||
                     iterPrev->Type()==BLOCK_OVERCALL ||
                     iterPrev->Type()==BLOCK_UNDERCALL
                   )
                {
                    if ( iterPrev->m_genome[0]==iter->m_read[0] &&
                         iterPrev->genomeLength()>=minHomopolymerSize &&
                         iterPrev->genomeLength()<=maxHomopolymerSize
                       )
                    {
                        m_homopolymerLeftCloseToInsert[iterPrev->genomeLength()] += 1;
                    }

                }
            }

            if (distance(iter, alignBlocks.end())>2)
            {
                iterNext = next(iter, 2);
                if ( iterNext->Type()==BLOCK_MATCH ||
                     iterNext->Type()==BLOCK_OVERCALL ||
                     iterNext->Type()==BLOCK_UNDERCALL
                   )
                {
                    if ( iterNext->m_genome[0]==iter->m_read[0] &&
                         iterNext->genomeLength()>=minHomopolymerSize &&
                         iterNext->genomeLength()<=maxHomopolymerSize
                       )
                    {
                        m_homopolymerRightCloseToInsert[iterNext->genomeLength()] += 1;
                    }
                }
            }

        }
    }

}


void GenericDataStatistics::printMatchMismatch(ostream& output)
{
    output << "[Mismatch/Total]" << endl;

    for (int i=0; i<m_bins.m_binLabels.size(); i++)
    {
        output
               << m_bins.m_binLabels[i] << " "
               << m_binMismatchCount[i] << "(" << (m_binMismatchCount[i]/m_binCount[i]) << ")" << " "
               << m_binCount[i] << endl;
    }
}


void GenericDataStatistics::printHomopolymerGap(ostream &output)
{
    output << "[Undercall/Overcall/Total]" << endl;

    for (int i=minHomopolymerSize; i<=maxHomopolymerSize; i++)
    {
        output
               << i << " "
               << m_homopolymerDelete[i] << "(" << (m_homopolymerDelete[i]/m_homopolymerCount[i]) << ")" << " "
               << m_homopolymerInsert[i] << "(" << (m_homopolymerInsert[i]/m_homopolymerCount[i]) << ")" << " "
               << m_homopolymerCount[i] << endl;
    }
}


void GenericDataStatistics::printDelete(ostream &output)
{
    output << "[Delete/Count(Site in Homopolymer)]" << endl;
    for (int i=minHomopolymerSize; i<=maxHomopolymerSize; i++)
    {
        output
                << i << " ";
        for (int j=0; j<i; j++)
        {
            output << m_homopolymerSiteDelete[i][j] << "/"
                   << m_homopolymerSiteCount[i][j]
                   << "(" << (m_homopolymerSiteDelete[i][j]/(m_homopolymerSiteCount[i][j]+1e-9)) << ")"
                   << " ";
        }

        output << endl;
    }
}


void GenericDataStatistics::printInsert(ostream &output)
{
    long double hntZ = 0;
    long double lciZ = 0;
    long double rciZ = 0;
    long double hZ   = 0;

    for (VectorDouble::iterator iter=m_homopolymerNextToInsert.begin();
         iter!=m_homopolymerNextToInsert.end(); ++iter)
        hntZ += *iter;

    for (VectorDouble::iterator iter=m_homopolymerLeftCloseToInsert.begin();
         iter!=m_homopolymerLeftCloseToInsert.end(); ++iter)
        lciZ += *iter;

    for (VectorDouble::iterator iter=m_homopolymerRightCloseToInsert.begin();
         iter!=m_homopolymerRightCloseToInsert.end(); ++iter)
        rciZ += *iter;

    for (VectorDouble::iterator iter=m_homopolymerCount.begin();
         iter!=m_homopolymerCount.end(); ++iter)
        hZ += *iter;


    output << "[Insert/Count(Homopolymer Context)]" << endl;
    for (int i=minHomopolymerSize; i<=maxHomopolymerSize; i++)
    {
        output << i << " "
               << m_homopolymerNextToInsert[i] << "(" << m_homopolymerNextToInsert[i]/hntZ << ")" << " "
               << m_homopolymerLeftCloseToInsert[i] << "(" << m_homopolymerLeftCloseToInsert[i]/lciZ << ")" << " "
               << m_homopolymerRightCloseToInsert[i] << "(" << m_homopolymerRightCloseToInsert[i]/rciZ << ")" << " "
               << m_homopolymerCount[i] << "(" << m_homopolymerCount[i]/hZ << ")"
               << endl;
    }
}
