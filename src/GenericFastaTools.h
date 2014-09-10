/*
 * This is a set of tools in extracting information of a sequence.
 *
 * TODO:
 * (1) A generic method to find universal short tandem repeats in a fasta object.
 * (2) Compute the entropy of a fasta object.
 *
 * Feng Zeng @2014-08-23
*/

#ifndef GENERICFASTATOOLS_H
#define GENERICFASTATOOLS_H

#include <string>
#include <tuple>
#include <vector>
using namespace std;

#include "GenericSequenceGlobal.h"
using namespace GenericSequenceTools;

namespace GenericSequenceTools
{

class GenericFastaTools
{
    public:
        GenericFastaTools();

    public:
        // mark the length of the homopolymer preceding current base
        static void markPrecedeHomopolymer(string& seq, vector<int>& len);
        static void markPrecedeHomopolymer(vector<int>& seq, vector<int>& len);

        // mark the length of the homopolymer succeding current base
        static void markSuccedeHomopolymer(string& seq, vector<int>& len);
        static void markSuccedeHomopolymer(vector<int>& seq, vector<int>& len);

        // mark the length of the homopolymer at left neighborhood
        static void markLeftNeighborHomopolymer(string& seq, vector<int>& len);

        // mark the length of the homopolymer at right neighborhood
        static void markRightNeighborHomopolymer(string& seq, vector<int>& len);



    public:
        // find mononucleotide repeat in sequence
        static int findMononucleotideRepeat(string& seq, int minRepeatLength, vector<int>& repeatPos, vector<int>& repeatSize);

        // find dinucleotide repeat in sequence
        static int findDinucleotideRepeat(string& seq, int minRepeatLength, vector<int>& repeatPos, vector<int>& repeatSize);
        static int findDinucleotideRepeatByFrameshift(string &seq, int minRepeatLength, int shift, vector<int> &repeatPos, vector<int> &repeatSize);

        // find trinucleotide repeat in sequence
        static int findTrinucleotideRepeat(string& seq, vector<int>& repeatPos, vector<int>& repeatSize);

        // find tandem repeat regions in sequence
        static int findTandemRepeatRegions(string& seq, int minRepeatLength, vector<int>& repeatStartPos, vector<int>& repeatLength);

    public:
        // the number of ambiguous base
        static int numAmbiguousBase(string& seq);

    public:
        // sequence signature
        static string sequenceSignature(string& seq);

        // maximum homopolymer length
        static int maxHomopolymerLength(string& seq);
};

}

#endif // GENERICFASTATOOLS_H
