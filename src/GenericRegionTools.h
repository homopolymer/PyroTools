/*
 * A set of tools to parse the region option.
 * Support region format
 *      (1) chr
 *      (2) chr:pos
 *      (3) chr:leftPos..rightPos
 *      (4) chr:leftPos-rightPos
 *
 * Created by Feng Zeng @20140824
*/

#ifndef GENERICREGIONTOOLS_H
#define GENERICREGIONTOOLS_H

#include <tuple>
#include <vector>
#include <string>
using namespace std;

#include "api/BamAux.h"
using namespace BamTools;

namespace GenericSequenceTools
{

class GenericRegionTools
{
public:
    GenericRegionTools();

public:
    // parse a region string
    static void parse(string& Roi, string& chrName, int& chrLeftPosition, int& chrRightPosition);

    // parse a region file
    static void parse(string& RoiFile, vector<string>& chrNames, vector<int>& chrLeftPositions, vector<int>& chrRightPositions);

    // convert a set of region strings to genome positions
    static int toScanWindow(RefVector& Genome, vector<string>& RoiSet, int WindowSize, vector<tuple<int,int,int>>& WindowSet);
    static int toScanWindow(RefVector& Genome, vector<string>& RoiSet, vector<tuple<int,int,int>>& WindowSet);

};

}   // namespace

#endif // GENERICREGIONTOOLS_H
