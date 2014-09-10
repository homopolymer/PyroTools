#include "GenericFastaTools.h"
#include "GenericSequenceGlobal.h"
using namespace GenericSequenceTools;

#include <set>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
using namespace std;


GenericFastaTools::GenericFastaTools()
{
}


void GenericFastaTools::markPrecedeHomopolymer(string &seq, vector<int> &len)
{
    len.assign(seq.length(), 0);

    for (int i=1; i<seq.length(); ++i)
    {
        if (seq[i]==seq[i-1])
            len[i] = len[i-1] + 1;
    }
}

void GenericFastaTools::markPrecedeHomopolymer(vector<int> &seq, vector<int> &len)
{
    len.assign(seq.size(), 0);

    for (int i=1; i<seq.size(); ++i)
    {
        if (seq[i]==seq[i-1])
            len[i] = len[i-1] + 1;
    }
}

void GenericFastaTools::markSuccedeHomopolymer(string &seq, vector<int> &len)
{
    string rev_seq(seq.rbegin(), seq.rend());
    vector<int> rev_len;
    markPrecedeHomopolymer(rev_seq, rev_len);

    len.assign(rev_len.rbegin(), rev_len.rend());
}

void GenericFastaTools::markSuccedeHomopolymer(vector<int> &seq, vector<int> &len)
{
    vector<int> rev_seq(seq.rbegin(), seq.rend());
    vector<int> rev_len;
    markPrecedeHomopolymer(rev_seq, rev_len);

    len.assign(rev_len.rbegin(), rev_len.rend());
}

void GenericFastaTools::markLeftNeighborHomopolymer(string &seq, vector<int> &len)
{
    // precede homopolymer
    vector<int> precedeHomopolymerLength;
    markPrecedeHomopolymer(seq, precedeHomopolymerLength);

    for (vector<int>::iterator iter=precedeHomopolymerLength.begin();
         iter!=precedeHomopolymerLength.end(); ++iter)
    {
        (*iter) += 1;
    }

    // neighbor homopolymer
    vector<int>  neighborHomopolymerLength(seq.length(), 0);
    vector<char> neighborHomopolymarNucl(seq.length(), Nil);

    for (int i=1; i<seq.length(); ++i)
    {
        if (seq[i]==seq[i-1])
        {
            neighborHomopolymarNucl[i] = neighborHomopolymarNucl[i-1];
            neighborHomopolymerLength[i] = neighborHomopolymerLength[i-1];
        }else
        {
            neighborHomopolymarNucl[i] = seq[i-1];
            neighborHomopolymerLength[i] = precedeHomopolymerLength[i-1];
        }
    }

    // Markovian homopolymer
    len.assign(seq.length(), 0);
    for (int i=1; i<seq.length(); ++i)
    {
        if (seq[i]==seq[i-1])
        {
            len[i] = len[i-1];
        }else
        {
            if (seq[i]==neighborHomopolymarNucl[i-1])
                len[i] = neighborHomopolymerLength[i-1];
        }
    }
}

void GenericFastaTools::markRightNeighborHomopolymer(string &seq, vector<int> &len)
{
    string rev_seq(seq.rbegin(), seq.rend());
    vector<int> rev_len;

    markLeftNeighborHomopolymer(rev_seq, rev_len);

    len.assign(rev_len.rbegin(), rev_len.rend());
}

int GenericFastaTools::numAmbiguousBase(string &seq)
{
    int num = 0;
    for (string::iterator iter=seq.begin(); iter!=seq.end(); ++iter)
    {
        if ((*iter)==Amb)
            num++;
    }

    return num;
}


int regionInclusion(int rA0, int rAl, int rB0, int rBl)
{
    return ((rA0<rB0) && (rA0+rAl>rB0+rBl));
}

int regionDistance(int rA0, int rAl, int rB0, int rBl)
{
    return (rB0-rA0-rAl);
}

struct PositionOrder
{
    bool operator()(const tuple<int,int>& a, const tuple<int,int>& b) const
    {
        return (get<0>(a)<get<0>(b));
    }
};

void sortAndMergeRepeats(vector<tuple<int,int>>& homopolymerList)
{
    sort(homopolymerList.begin(), homopolymerList.end(), PositionOrder());

    vector<tuple<int,int>> mergeHomopolymerList;

    for (int i=0; i<homopolymerList.size(); i++)
    {
        if (i==0)
        {
            mergeHomopolymerList.push_back(homopolymerList[i]);
            continue;
        }

        int j = mergeHomopolymerList.size()-1;
        int rA0 = get<0>(mergeHomopolymerList[j]);
        int rAl = get<1>(mergeHomopolymerList[j]);

        int rB0 = get<0>(homopolymerList[i]);
        int rBl = get<1>(homopolymerList[i]);

        if (regionInclusion(rA0, rAl, rB0, rBl))
        {
             continue;
        }

        if (regionDistance(rA0, rAl, rB0, rBl)<5)
        {
            get<1>(mergeHomopolymerList[j])   = rB0+rBl-rA0;
            continue;
        }

        mergeHomopolymerList.push_back(homopolymerList[i]);
    }


    homopolymerList.clear();
    homopolymerList.assign(mergeHomopolymerList.begin(), mergeHomopolymerList.end());
}

void uniqueRepeats(vector<int>& inRepeatPos, vector<int>& inRepeatSize, vector<int>& outRepeatPos, vector<int>& outRepeatSize)
{
    // sort
    vector<tuple<int,int>> homopolymerList(inRepeatPos.size());
    for (int i=0; i<inRepeatPos.size(); ++i)
    {
        homopolymerList[i] = tuple<int,int>(inRepeatPos[i], inRepeatSize[i]);
    }

    sortAndMergeRepeats(homopolymerList);

    // unique
    set<int> hasExist;
    for (int i=0; i<homopolymerList.size(); ++i)
    {
        int p,l;
        p = get<0>(homopolymerList[i]);
        l = get<1>(homopolymerList[i]);

        if (l<0 || l>1000)
            continue;

        if (hasExist.find(p)==hasExist.end())
        {
            outRepeatPos.push_back(p);
            outRepeatSize.push_back(l);
            hasExist.insert(p);
        }
    }
}

int GenericFastaTools::findTandemRepeatRegions(string &seq, int minRepeatLength, vector<int> &repeatStartPos, vector<int> &repeatLength)
{
    // mononucleotide repeat regions
    int numMonoRepeat;
    vector<int> monoRepeatPos, monoRepeatSize;
    numMonoRepeat = findMononucleotideRepeat(seq, minRepeatLength, monoRepeatPos, monoRepeatSize);

    // dinucleotide repeat regions
    int numDiRepeat;
    vector<int> diRepeatPos, diRepeatSize;
    numDiRepeat = findDinucleotideRepeat(seq, minRepeatLength, diRepeatPos, diRepeatSize);


    // merge together
    if (numMonoRepeat>0)
    {
        repeatStartPos.insert(repeatStartPos.end(), monoRepeatPos.begin(), monoRepeatPos.end());
        repeatLength.insert(repeatLength.end(), monoRepeatSize.begin(), monoRepeatSize.end());
    }

    if (numDiRepeat>0)
    {
        repeatStartPos.insert(repeatStartPos.end(), diRepeatPos.begin(), diRepeatPos.end());
        repeatLength.insert(repeatLength.end(), diRepeatSize.begin(), diRepeatSize.end());
    }

    // sort, merge and unique
    vector<int> t_repeatPos = repeatStartPos;
    vector<int> t_repeatSize = repeatLength;
    repeatStartPos.clear();
    repeatLength.clear();
    uniqueRepeats(t_repeatPos, t_repeatSize, repeatStartPos, repeatLength);

    // return value
    return repeatStartPos.size();
}

int GenericFastaTools::findMononucleotideRepeat(string &seq, int minRepeatLength, vector<int> &repeatPos, vector<int> &repeatSize)
{
    // homopolymer precede current position
    vector<int> precedeHomopolymerLength;
    markPrecedeHomopolymer(seq, precedeHomopolymerLength);

    // homopolymer succede current position
    vector<int> succedeHomopolymerLength;
    markSuccedeHomopolymer(seq, succedeHomopolymerLength);

    // homopolymer left close to current position
    vector<int> leftCloseHomopolymerLength;
    markLeftNeighborHomopolymer(seq, leftCloseHomopolymerLength);

    // homopolymer right close to current position
    vector<int> rightCloseHomopolymerLength;
    markRightNeighborHomopolymer(seq, rightCloseHomopolymerLength);

    // homopolymer
    vector<int> position;
    vector<int> length;

    // homopolymer region length
    vector<int> homopolymerLength(seq.length(), 0);
    for (int i=0; i<seq.length(); ++i)
    {
        homopolymerLength[i] = 1 + precedeHomopolymerLength[i] + succedeHomopolymerLength[i];
    }

    // iterate from first region to last region
    for (int i=0; i<seq.length(); ++i)
    {
        if (homopolymerLength[i]>=minRepeatLength && precedeHomopolymerLength[i]==0)
        {
            position.push_back(i);
            length.push_back(homopolymerLength[i]);
        }
    }

    // extend left
    vector<int> positionLeft(position.begin(), position.end());
    vector<int> lengthLeft(length.size(), 0);
    for (int j=0; j<position.size(); ++j)
    {
        int i = position[j];
        while (leftCloseHomopolymerLength[i]>1)
        {
            int t_i;
            t_i = i - homopolymerLength[i-1] - leftCloseHomopolymerLength[i];
            i =  t_i;
        }

        positionLeft[j] =  i;
        lengthLeft[j]   =  position[j]-i;
    }

    // extend right
    vector<int> lengthRight(length.size(), 0);
    for (int j=0; j<position.size(); ++j)
    {
        int i = position[j];
        while (rightCloseHomopolymerLength[i]>1)
        {
            int t_i;
            t_i = i + homopolymerLength[i] + homopolymerLength[i + homopolymerLength[i]];
            i = t_i;
        }

        lengthRight[j] = i-position[j]+homopolymerLength[i]-homopolymerLength[position[j]];
    }

    // combine left and right results
    for (int i=0; i<position.size(); ++i)
    {
        position[i] = positionLeft[i];
        length[i]   = length[i] + lengthLeft[i] + lengthRight[i];
    }

    // sort, merge and unique
    vector<int> t_repeatPos(position);
    vector<int> t_repeatSize(length);
    uniqueRepeats(t_repeatPos, t_repeatSize, repeatPos, repeatSize);

    return (repeatPos.size());
}


map<string,int> DinucleotideRepeatPatterns = {{"AC",0},{"AG",1},{"AT",2},
                                              {"CA",3},{"CG",4},{"CT",5},
                                              {"GA",6},{"GC",7},{"GT",8},
                                              {"TA",9},{"TC",10},{"TG",11},
                                              {"N",12}};

int GenericFastaTools::findDinucleotideRepeat(string &seq, int minRepeatLength, vector<int> &repeatPos, vector<int> &repeatSize)
{
    vector<int> repeatPosFrameZero, repeatPosFrameOne;
    vector<int> repeatSizeFrameZero, repeatSizeFrameOne;

    // shift - 0
    int maxRepeatSizeFrameZero = findDinucleotideRepeatByFrameshift(seq, minRepeatLength, 0, repeatPosFrameZero, repeatSizeFrameZero);
    // shift - 1
    int maxRepeatSizeFrameOne  = findDinucleotideRepeatByFrameshift(seq, minRepeatLength, 1, repeatPosFrameOne, repeatSizeFrameOne);

    // the largest repeat as the pick
    if (maxRepeatSizeFrameZero>maxRepeatSizeFrameOne)
    {
        repeatPos  = repeatPosFrameZero;
        repeatSize = repeatSizeFrameZero;
    }else
    {
        repeatPos  = repeatPosFrameOne;
        repeatSize = repeatSizeFrameOne;
    }

    return repeatPos.size();
}

int GenericFastaTools::findDinucleotideRepeatByFrameshift(string &seq, int minRepeatLength, int shift, vector<int> &repeatPos, vector<int> &repeatSize)
{
    int maxRepeatSize;

    // transform original sequence to new form
    vector<int> DiSeq;
    for (int i=shift; i<seq.length(); i+=2)
    {
        string di = seq.substr(i,2);
        if (di[0]=='N' || di[1]=='N')
        {
            DiSeq.push_back(DinucleotideRepeatPatterns["N"]);
        }else
        {
            DiSeq.push_back(DinucleotideRepeatPatterns[di]);
        }
    }

    // mark precede homopoymer in dinucleotide repeat space
    vector<int> DiSeqHomopolymerSize;
    markSuccedeHomopolymer(DiSeq, DiSeqHomopolymerSize);

    // save repeats
    maxRepeatSize = 0;
    for (int i=0; i<DiSeqHomopolymerSize.size(); )
    {
        if (DiSeqHomopolymerSize[i]+1>=minRepeatLength)
        {
            int repeatPos0  = shift + 2*i;
            int repeatSize0 = 2*(DiSeqHomopolymerSize[i]+1);

            repeatPos.push_back(repeatPos0);
            repeatSize.push_back(repeatSize0);

            if (maxRepeatSize<repeatSize0)
                maxRepeatSize = repeatSize0;
        }

        i += DiSeqHomopolymerSize[i]+1;
    }

    return maxRepeatSize;
}

string GenericFastaTools::sequenceSignature(string &seq)
{
    string sig="";

    for (int i=0; i<seq.length(); ++i)
    {
        if (i==0)
            sig += seq[i];
        else
        {
            if (sig[sig.length()-1]!=seq[i])
                sig += seq[i];
        }
    }

    return sig;
}

int GenericFastaTools::maxHomopolymerLength(string &seq)
{
    vector<int> hl;
    markPrecedeHomopolymer(seq, hl);

    int ml=0;
    for (vector<int>::iterator iter=hl.begin(); iter!=hl.end(); ++iter)
    {
        if (ml<(*iter))
            ml = *iter;
    }

    ml++;
    return ml;
}
