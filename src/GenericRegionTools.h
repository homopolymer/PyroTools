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

#include <vector>
#include <string>
using namespace std;

namespace GenericSequenceTools
{

class GenericRegionTools
{
public:
    GenericRegionTools();

public:
    // parse a region string
    static void parse(string& Roi, int& chrID, int& chrLeftPosition, int& chrRightPosition);

    // parse a region file
    static void parse(string& RoiFile, vector<int>& chrIDs, vector<int>& chrLeftPositions, vector<int>& chrRightPositions);
};

}   // namespace

#endif // GENERICREGIONTOOLS_H
