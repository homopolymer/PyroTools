#include "GenericRegionTools.h"
using namespace GenericSequenceTools;

#include <iostream>
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
void parseRoiStyle0(string& Roi, int& chrID, int& chrLeftPosition, int& chrRightPosition)
{
    chrLeftPosition  = 0;
    chrRightPosition = -1;
    chrID = stoi(Roi);
}

// chr:leftPos..rightPos
void parseRoiStyle1(string& Roi, int& chrID, int& chrLeftPosition, int& chrRightPosition)
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
    string CHRID = Roi.substr(0, colonPos);
    chrID = stoi(CHRID);

    // retrieve left position
    string CHRLEFTPOS = Roi.substr(colonPos+1, dot0Pos-colonPos-1);
    chrLeftPosition = stoi(CHRLEFTPOS);

    // retrieve right position
    string CHRRIGHTPOSITION = Roi.substr(dot1Pos+1, Roi.length()-dot1Pos-1);
    chrRightPosition = stoi(CHRRIGHTPOSITION);
}

// chr:leftPos-rightPos
void parseRoiStyle2(string& Roi, int& chrID, int& chrLeftPosition, int& chrRightPosition)
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
    string CHRID = Roi.substr(0, colonPos);
    chrID = stoi(CHRID);

    // retrieve chr leftpos
    string CHRLEFTPOSITION = Roi.substr(colonPos+1, dashPos-colonPos-1);
    chrLeftPosition = stoi(CHRLEFTPOSITION);

    // retrieve chr rightpos
    string CHRRIGHTPOSITION = Roi.substr(dashPos+1, Roi.length()-dashPos-1);
    chrRightPosition = stoi(CHRRIGHTPOSITION);
}

void GenericRegionTools::parse(string &Roi, int& chrID, int& chrLeftPosition, int& chrRightPosition)
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
        parseRoiStyle0(Roi, chrID, chrLeftPosition, chrRightPosition);

    if (hasColon && hasDotDot)
        parseRoiStyle1(Roi, chrID, chrLeftPosition, chrRightPosition);

    if (hasColon && hasDash)
        parseRoiStyle2(Roi, chrID, chrLeftPosition, chrRightPosition);
}
