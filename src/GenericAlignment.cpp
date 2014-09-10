#include "GenericAlignment.h"
#include "GenericSequenceGlobal.h"
using namespace GenericSequenceTools;

#include <sstream>
using namespace std;

GenericAlignment::GenericAlignment()
{
}


int NW_STATE_SIZE = 4;
int NW_MATCH      = 0;
int NW_MISMATCH   = 1;
int NW_INSERT     = 2;
int NW_DELETE     = 3;

double NeedlemanWunschLocalScore(string& seqA, string& seqB, int i, int j, int s, int t, double matchScore, double mismatchScore, double gapOpenScore, double gapExtendScore)
{
    double x;

    if (s==NW_MATCH && seqA[i]!=seqB[j])
        return DOUBLE_NEGATIVE_INFINITY;

    if (s==NW_MISMATCH && seqA[i]==seqB[j])
        return DOUBLE_NEGATIVE_INFINITY;

    if (s==NW_MATCH && seqA[i]==seqB[j])
        x = matchScore;

    if (s==NW_MISMATCH && seqA[i]!=seqB[j])
        x = mismatchScore;

    if (s==NW_INSERT)
    {
        if (t==NW_INSERT)
            x = gapExtendScore;
        else
            x = gapOpenScore+gapExtendScore;

    }

    if (s==NW_DELETE)
    {
        if (t==NW_DELETE)
            x = gapExtendScore;
        else
            x = gapOpenScore+gapExtendScore;
    }

    return x;
}

void GenericAlignment::NeedlemanWunsch(string &seqA, string &seqB, string &alnSeqA, string &alnSeqB, double matchScore, double mismatchScore, double gapOpenScore, double gapExtendScore)
{
    int m = seqA.length();
    int n = seqB.length();

    // two vectors for boundaries
    VectorDouble dpRowBoundary(m);
    VectorDouble dpColBoundary(n);

    // compute scores at the boundary
    for (int i=0; i<m; ++i)
    {
        if (i==0)
        {
            dpRowBoundary[i] = gapOpenScore + gapExtendScore;
        }else
        {
            dpRowBoundary[i] = dpRowBoundary[i-1] + gapExtendScore;
        }
    }

    for (int j=0; j<n; ++j)
    {
        if (j==0)
        {
            dpColBoundary[j] = gapOpenScore + gapExtendScore;
        }else
        {
            dpColBoundary[j] = dpColBoundary[j-1] + gapExtendScore;
        }
    }

    // matrices for dynamic programming
    Matrix3Double  NWAlignmentScore;
    setMatrixValue(NWAlignmentScore, m, n, NW_STATE_SIZE, DOUBLE_NEGATIVE_INFINITY);

    Matrix3Integer NWPreviousState;
    setMatrixValue(NWPreviousState, m, n, NW_STATE_SIZE, -1);

    Matrix3Integer NWPreviousI;
    setMatrixValue(NWPreviousI, m, n, NW_STATE_SIZE, -1);

    Matrix3Integer NWPreviousJ;
    setMatrixValue(NWPreviousJ, m, n, NW_STATE_SIZE, -1);

    // dynamic programming computation
    for (int i=0; i<m; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            for (int s=NW_MATCH; s<NW_STATE_SIZE; s++)
            {
                int mi,mj,mt;
                double mx = DOUBLE_NEGATIVE_INFINITY;

                if (s==NW_MATCH && seqA[i]!=seqB[j])
                    continue;

                if (s==NW_MISMATCH && seqA[i]==seqB[j])
                    continue;

                // the positions of the preceding alignment
                int pi,pj;
                if (s==NW_MATCH || s==NW_MISMATCH)
                {
                    pi = i-1;
                    pj = j-1;
                }

                if (s==NW_INSERT)
                {
                    pi = 0;
                    pj = j-1;
                }

                if (s==NW_DELETE)
                {
                    pi = i-1;
                    pj = 0;
                }

                // scoring the alignment
                for (int t=NW_MATCH; t<NW_STATE_SIZE; t++)
                {
                    double x = 0;

                    if (pi<0 && pj<0)
                    {
                        if (t!=NW_MATCH)
                            continue;

                        x = NeedlemanWunschLocalScore(seqA, seqB, i, j, s, t, matchScore, mismatchScore, gapOpenScore, gapExtendScore);
                    }else if (pi<0 && pj>=0)
                    {
                        if (t!=NW_INSERT)
                            continue;

                        x = dpColBoundary[pj] + NeedlemanWunschLocalScore(seqA, seqB, i, j, s, t, matchScore, mismatchScore, gapOpenScore, gapExtendScore);
                    }else if (pi>=0 && pj<0)
                    {
                        if (t!=NW_DELETE)
                            continue;

                        x = dpRowBoundary[pi] + NeedlemanWunschLocalScore(seqA, seqB, i, j, s, t, matchScore, mismatchScore, gapOpenScore, gapExtendScore);
                    }else
                    {
                        x = NWAlignmentScore[pi][pj][t] + NeedlemanWunschLocalScore(seqA, seqB, i, j, s, t, matchScore, mismatchScore, gapOpenScore, gapExtendScore);
                    }

                    if (mx<x)
                    {
                        mx = x;
                        mi = pi;
                        mj = pj;
                        mt = t;
                    }
                }

                NWAlignmentScore[i][j][s] = mx;
                NWPreviousState[i][j][s]  = mt;
                NWPreviousI[i][j][s]      = mi;
                NWPreviousJ[i][j][s]      = mj;
            }
        }
    }

    // backtracing
    double X=DOUBLE_NEGATIVE_INFINITY;
    int I,J,S;
    for (int s=NW_MATCH; s<NW_STATE_SIZE; ++s)
    {
        if (X<NWAlignmentScore[m-1][n-1][s])
        {
            X = NWAlignmentScore[m-1][n-1][s];
            S = s;
        }
    }

    I = m-1;
    J = n-1;
    bool trace=true;
    while (trace)
    {
        if (S==NW_MATCH || S==NW_MISMATCH)
        {
            stringstream alnSeqABuf, alnSeqBBuf;

            alnSeqABuf << seqA[I] << alnSeqA;
            alnSeqA = alnSeqABuf.str();

            alnSeqBBuf << seqB[J] << alnSeqB;
            alnSeqB = alnSeqBBuf.str();
        }

        if (S==NW_INSERT)
        {
            stringstream alnSeqABuf, alnSeqBBuf;

            alnSeqABuf << "-" << alnSeqA;
            alnSeqA = alnSeqABuf.str();

            alnSeqBBuf << seqB[J] << alnSeqB;
            alnSeqB = alnSeqBBuf.str();
        }

        if (S==NW_DELETE)
        {
            stringstream alnSeqABuf, alnSeqBBuf;

            alnSeqABuf << seqA[I] << alnSeqA;
            alnSeqA = alnSeqABuf.str();

            alnSeqBBuf << "-" << alnSeqB;
            alnSeqB = alnSeqBBuf.str();
        }

        int tI = NWPreviousI[I][J][X];
        int tJ = NWPreviousJ[I][J][X];
        int tX = NWPreviousState[I][J][X];

        if (tI<0 && tJ<0)
            trace = false;
        if (tI<0 && tX!=NW_INSERT)
            trace = false;
        if (tJ<0 && tX!=NW_DELETE)
            trace = false;
    }
}
