#ifndef SIMPLEVARIANTCALLTOOL_H
#define SIMPLEVARIANTCALLTOOL_H

#include "GenericTool.h"
using namespace GenericSequenceTools;

#include <vector>
#include <string>
using namespace std;

namespace GenericSequenceTools{

class SimpleVariantCallTool : public GenericAbstractTool
{
public:
    SimpleVariantCallTool();

public:
    int Run(int argc, char *argv[]);
    int Help();

public:
    int SimpleVariantCall();

public:
    int commandLineParser(int argc, char *argv[]);

public:
    string         genomeFile;
    vector<string> bamFiles;

    vector<string> regionStringsOfInterest;
    int            mapQualThres;
    int            readLenThres;
    float          alnIdentityThres;
    int            numAmbiguousThres;
    int            alnFlagMarker;

    int            baseQualThres;
    int            snpHitThres;
    int            indelHitThres;
    bool           ignoreSnp;
    bool           ignoreIndel;
    string         outputFormat;

    int            numThreads;
    int            windowSize;

    int            verbose;
};
}   // namespace

#endif // SIMPLEVARIANTCALLTOOL_H
