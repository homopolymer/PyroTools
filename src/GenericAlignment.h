#ifndef GENERICALIGNMENT_H
#define GENERICALIGNMENT_H

#include <string>
using namespace std;

namespace GenericSequenceTools
{

class GenericAlignment
{
    public:
        GenericAlignment();

    public:
        static void NeedlemanWunsch(string& seqA, string& seqB,
                                    string& alnSeqA, string& alnSeqB,
                                    double matchScore=2., double mismatchScore=-2.,
                                    double gapOpenScore=-3., double gapExtendScore=-1.);
};

}

#endif // GENERICALIGNMENT_H
