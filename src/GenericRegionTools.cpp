#include "GenericRegionTools.h"
using namespace GenericSequenceTools;

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <parallel/algorithm>
#include <algorithm>
#include <forward_list>
#include <unordered_map>
using namespace std;

GenericRegionTools::GenericRegionTools()
{
}

void printRegionFormat(ostream& out)
{
    out << "Valid region format:" << endl;
    out << "(1) chr" << endl;
    out << "(2) chr:leftPos..rightPos" << endl;
    out << "(3) chr:leftPos-rightPos" << endl;
}

// chr
void parseRoiStyle0(string& Roi, string& chrName, int& chrLeftPosition, int& chrRightPosition)
{
    chrLeftPosition  = 0;
    chrRightPosition = -1;
    chrName = Roi;
}

// chr:leftPos..rightPos
void parseRoiStyle1(string& Roi, string& chrName, int& chrLeftPosition, int& chrRightPosition)
{
    int colonPos=-1;
    int dot0Pos=-1;
    int dot1Pos=-1;

    for (int i=0; i<Roi.length(); i++)
    {
        if (Roi[i]==':')
        {
            if (colonPos>=0)
            {
                cerr << "Invalid region string " << Roi << endl;
                printRegionFormat(cerr);
                cerr << "Program exit" << endl;
                exit(-1);
            }
            else
            {
                colonPos = i;
            }
        }

        if (Roi[i]=='.')
        {
            if (dot0Pos<0)
            {
                dot0Pos = i;
            }

            dot1Pos = i;
        }
    }

    if (dot0Pos==dot1Pos)
    {
        cerr << "Invalid region string " << Roi << endl;
        printRegionFormat(cerr);
        cerr << "Program exit" << endl;
        exit(-1);
    }

    // retrieve chr id
    string CHRNAME = Roi.substr(0, colonPos);
    chrName = CHRNAME;

    // retrieve left position
    string CHRLEFTPOS = Roi.substr(colonPos+1, dot0Pos-colonPos-1);
    chrLeftPosition = stoi(CHRLEFTPOS)-1;

    // retrieve right position
    string CHRRIGHTPOSITION = Roi.substr(dot1Pos+1, Roi.length()-dot1Pos-1);
    chrRightPosition = stoi(CHRRIGHTPOSITION);
}

// chr:leftPos-rightPos
void parseRoiStyle2(string& Roi, string& chrName, int& chrLeftPosition, int& chrRightPosition)
{
    int colonPos = -1;
    int dashPos  = -1;

    for (int i=0; i<Roi.length(); i++)
    {
        if (Roi[i]==':')
        {
            if (colonPos>=0)
            {
                cerr << "Invalid region string " << Roi << endl;
                printRegionFormat(cerr);
                cerr << "Program exit" << endl;
                exit(-1);
            }
            else
            {
                colonPos = i;
            }
        }

        if (Roi[i]=='-')
        {
            if (dashPos>=0)
            {
                cerr << "Invalid region string " << Roi << endl;
                printRegionFormat(cerr);
                cerr << "Program exit" << endl;
                exit(-1);
            }
            else
            {
                dashPos = i;
            }
        }
    }

    // retrieve chr id
    string CHRNAME = Roi.substr(0, colonPos);
    chrName = CHRNAME;

    // retrieve chr leftpos
    string CHRLEFTPOSITION = Roi.substr(colonPos+1, dashPos-colonPos-1);
    chrLeftPosition = stoi(CHRLEFTPOSITION)-1;

    // retrieve chr rightpos
    string CHRRIGHTPOSITION = Roi.substr(dashPos+1, Roi.length()-dashPos-1);
    chrRightPosition = stoi(CHRRIGHTPOSITION);
}

// parse the singleton ROI
void GenericRegionTools::parse(string &Roi, string& chrName, int& chrLeftPosition, int& chrRightPosition)
{
    bool hasColon  = false;
    bool hasDotDot = false;
    bool hasDash   = false;

    if (Roi.find(":")!=string::npos)
        hasColon = true;

    if (Roi.find("..")!=string::npos)
        hasDotDot = true;

    if (Roi.find("-")!=string::npos)
        hasDash = true;

    if (!hasColon)
        parseRoiStyle0(Roi, chrName, chrLeftPosition, chrRightPosition);

    if (hasColon && hasDotDot)
        parseRoiStyle1(Roi, chrName, chrLeftPosition, chrRightPosition);

    if (hasColon && hasDash)
        parseRoiStyle2(Roi, chrName, chrLeftPosition, chrRightPosition);
}

// parse the BED file
void GenericRegionTools::parse(string &RoiFile, vector<string> &chrNames, vector<int> &chrLeftPositions, vector<int> &chrRightPositions)
{
    // open the region file
    ifstream infile(RoiFile.c_str());

    // read line by line
    for (string line; getline(infile,line);)
    {
        istringstream iss(line);
        string chrName;
        int chrLeftPos,chrRightPos;
        if (!(iss >> chrName >> chrLeftPos >> chrRightPos))
        {
            break;  // error
        }
        // save the region information
        chrNames.emplace_back(chrName);
        chrLeftPositions.emplace_back(chrLeftPos-1);
        chrRightPositions.emplace_back(chrRightPos);
    }

    // close the region file
    infile.close();
}

int GenericRegionTools::toScanWindow(RefVector &Genome, vector<string> &RoiSet, int WindowSize, vector<tuple<int, int, int> > &WindowSet)
{
    int numWindows = 0;

    // convert region string to region info
    forward_list<tuple<int,int,int>> regions;
    // lambda function for parsing ROI
    auto RoiParser = [&](string r)
    {
        int rid=0,rlp,rrp;
        string rname;
        GenericRegionTools::parse(r, rname, rlp, rrp);
        for(auto g : Genome){ if (g.RefName==rname) break; rid++; }
        if (rrp<0) rrp=Genome[rid].RefLength;
        regions.emplace_front(tuple<int,int,int>(rid,rlp,rrp));
    };
    // parsing ROI
    // serial processing
    if (RoiSet.empty())
    {
        for (size_t i=0; i<Genome.size(); i++)
        {
            regions.emplace_front(tuple<int,int,int>(i,0,Genome[i].RefLength));
        }
    }
    else if (RoiSet.size()<100 && !RoiSet.empty())
    {
        std::for_each(RoiSet.begin(), RoiSet.end(), RoiParser);
    }
    // parallel processing
    else
    {
        __gnu_parallel::for_each(RoiSet.begin(), RoiSet.end(), RoiParser);
    }

    // loop over regions, and divide region into a set of windows
    size_t windowSize=WindowSize;
    for (auto roi : regions)
    {
        int roiChrID = get<0>(roi);
        int roiChrLp = get<1>(roi);
        int roiChrRp = get<2>(roi);
        for (int a=roiChrLp+1; a<roiChrRp; a+=windowSize)
        {
            int b = a+windowSize;
            if (b > roiChrRp) b = roiChrRp;
            WindowSet.emplace_back(tuple<int,int,int>(roiChrID,a-1,b));

            numWindows++;
        }
    }


    return numWindows;

}


int GenericRegionTools::toScanWindow(RefVector &Genome, vector<string> &RoiSet, vector<tuple<int, int, int> > &WindowSet)
{
    int numWindows = 0;

    // convert region string to region info
    forward_list<tuple<int,int,int>> regions;

    // lambda function for parsing ROI
    auto RoiParser = [&](string r)
    {
        int rid=0,rlp,rrp;
        string rname;
        GenericRegionTools::parse(r, rname, rlp, rrp);
        for(auto g : Genome){ if (g.RefName==rname) break; rid++; }
        if (rrp<0) rrp=Genome[rid].RefLength;
        regions.emplace_front(tuple<int,int,int>(rid,rlp,rrp));
    };

    // parsing ROI
    // serial processing
    if (RoiSet.empty())
    {
        for (size_t i=0; i<Genome.size(); i++)
        {
            regions.emplace_front(tuple<int,int,int>(i,0,Genome[i].RefLength));
        }
    }
    else if (RoiSet.size()<100 && !RoiSet.empty())
    {
        std::for_each(RoiSet.begin(), RoiSet.end(), RoiParser);
    }
    // parallel processing
    else
    {
        __gnu_parallel::for_each(RoiSet.begin(), RoiSet.end(), RoiParser);
    }

    // loop over regions, and divide region into a set of windows
    for (auto roi : regions)
    {
        int roiChrID = get<0>(roi);
        int roiChrLp = get<1>(roi);
        int roiChrRp = get<2>(roi);

        WindowSet.emplace_back(tuple<int,int,int>(roiChrID,roiChrLp,roiChrRp));

        numWindows++;
    }


    return numWindows;

}
