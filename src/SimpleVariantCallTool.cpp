#include "SimpleVariantCallTool.h"
#include "GenericRegionTools.h"
#include "GenericVariant.h"
#include "GenericBamAlignmentTools.h"
#include "GenericSequenceGlobal.h"
using namespace GenericSequenceTools;

#include <getopt.h>
#include <time.h>
#include <tuple>
#include <utility>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>
#include <forward_list>
#include <unordered_map>
#include <parallel/algorithm>
#include <omp.h>
using namespace std;


#include "api/BamReader.h"
#include "api/BamMultiReader.h"
#include "api/BamAux.h"
#include "api/BamAlignment.h"
using namespace BamTools;

SimpleVariantCallTool::SimpleVariantCallTool()
    :mapQualThres(5)
    ,readLenThres(100)
    ,alnIdentityThres(0.9)
    ,numAmbiguousThres(2)
    ,alnFlagMarker(0)
    ,baseQualThres(10)
    ,snpHitThres(2)
    ,indelHitThres(2)
    ,ignoreSnp(false)
    ,ignoreIndel(false)
    ,outputFormat("bed")
    ,numThreads(1)
    ,windowSize(1e+5)
    ,verbose(0)
{
}

// verbose
inline void Verbose(string MSG)
{
    cerr << "[" << CurrentTime() << " PyroTools-SimpleVarCall] "  << MSG << endl;
}

// check the existence of a file
// reference: http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
inline bool file_exist(const string& filename)
{
    struct stat buffer;
    return (stat (filename.c_str(), &buffer) == 0);
}

// help message
int SimpleVariantCallTool::Help()
{
    cerr << "SYNOPSIS" << endl;
    cerr << "    PyroTools SimpleVarCall [OPTIONS] <GENOME> <BAM> [<BAM>...]" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    SimpleVarCall is to simply profile all suspect variants in the genome given sequencing data." << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "Read Filtration" << endl;
    cerr << "    --roi              filter reads not hit the region of interest" << endl;
    cerr << "    --roi-list         filter reads not hit the regions of interest" << endl;
    cerr << "    --mq               filter reads with mapping quality less than the value (default:5) [INT]" << endl;
    cerr << "    --len              filter reads with length less than the value (default:100) [INT]" << endl;
    cerr << "    --iden             filter reads with the identity less than the value (default:0.9) [FLT]" << endl;
    cerr << "    --na               filter reads with many ambiguous base (default:2) [INT]" << endl;
    cerr << "    --ff               filter reads with the specified flags, flag info could be found in SAM manual [INT]" << endl;
    cerr << "Variant Filtration" << endl;
    cerr << "    --bq               filter read alleles with base quality less than the value (default:10) [INT]" << endl;
    cerr << "    --snp-hit          filter SNPs with the hits less than the value (default:2) [INT]" << endl;
    cerr << "    --indel-hit        filter INDELs with the hits less than the value (default:2) [INT]" << endl;
    cerr << "Variant Output" << endl;
    cerr << "    --disable-snp      ignore SNPs" << endl;
    cerr << "    --disable-indel    ignore INDELs" << endl;
    cerr << "    --out-format       specify the format of output, valid arguments: bed or vcf (default:bed) [STR]" << endl;
    cerr << "Miscellaneous" << endl;
    cerr << "    --window           specify the size of the scanning window (default:1e+5) [INT]" << endl;
    cerr << "    --num-cores        specify the number of cpu cores being used (default:1) [INT]" << endl;
    cerr << "    -v,--verbose       specify the level of verbosity, valid arguments: 0 for quiet, 1 for debug (default:0) [INT]" << endl;
    cerr << "    -h,--help          print help message" << endl;
    return 0;
}

int SimpleVariantCallTool::commandLineParser(int argc, char *argv[])
{
    auto isHexString = [](string s)
    {
        for (auto c : s)
        {
            if (c=='x') return true;
            if (c=='X') return true;
        }
        return false;
    };

    int c;

    while (1)
    {
        static struct option long_options[] =
        {
            {"roi",           required_argument, 0, 0},
            {"roi-list",      required_argument, 0, 0},
            {"mq",            required_argument, 0, 0},
            {"len",           required_argument, 0, 0},
            {"iden",          required_argument, 0, 0},
            {"na",            required_argument, 0, 0},
            {"ff",            required_argument, 0, 0},

            {"bq",            required_argument, 0, 0},
            {"snp-hit",       required_argument, 0, 0},
            {"indel-hit",     required_argument, 0, 0},

            {"disable-snp",   no_argument,       0, 0},
            {"disable-indel", no_argument,       0, 0},
            {"out-format",    required_argument, 0, 0},

            {"num-cores",     required_argument, 0, 0},
            {"window",        required_argument, 0, 0},
            {"verbose",       required_argument, 0, 'v'},
            {"help",          no_argument,       0, 'h'},
            {0,0,0,0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "v:h", long_options, &option_index);

        // detect the end of the options
        if (c==-1) break;

        switch(c)
        {
        case 0:
            switch(option_index)
            {
            case 0:
                regionStringsOfInterest.emplace_back(optarg);
                break;
            case 1:
                {
                    ifstream infile(optarg);
                    string line;
                    while (getline(infile,line))
                    {
                        if (!line.empty())
                            regionStringsOfInterest.emplace_back(optarg);
                    }
                    infile.close();
                }
                break;
            case 2:
                mapQualThres=stoi(optarg);
                break;
            case 3:
                readLenThres=stoi(optarg);
                break;
            case 4:
                alnIdentityThres=stod(optarg);
                break;
            case 5:
                numAmbiguousThres=stoi(optarg);
                break;
            case 6:
                if (isHexString(optarg)){
                    alnFlagMarker |= stoul(optarg,nullptr,16);
                }else{
                    alnFlagMarker |= stoi(optarg);
                }
                break;
            case 7:
                baseQualThres=stoi(optarg);
                break;
            case 8:
                snpHitThres=stoi(optarg);
                break;
            case 9:
                indelHitThres=stoi(optarg);
                break;
            case 10:
                ignoreSnp=true;
                break;
            case 11:
                ignoreIndel=true;
                break;
            case 12:
                outputFormat=optarg;
                break;
            case 13:
                numThreads=stoi(optarg);
                break;
            case 14:
                windowSize=stod(optarg);
                break;
            default:
                abort();
            }
            break;
        case 'v':
            verbose=stoi(optarg);
            break;
        case 'h':
            Help();
            exit(EXIT_SUCCESS);
            break;
        case '?':
            exit(EXIT_FAILURE);
            break;
        default:
            abort();
        }
    }

    // genome file
    if (optind<argc)
    {
        genomeFile=argv[optind++];

        // check the existence of the genome file
        if (!file_exist(genomeFile))
        {
            cerr << "[PyroTools-SimpleVarCall] error: "
                   << genomeFile << " not existed" << endl;
            exit(EXIT_FAILURE);
        }

        // check the existence of the genome index file
        if (!file_exist(genomeFile+".fai"))
        {
            cerr << "[PyroTools-SimpleVarCall] error: "
                   << (genomeFile+".fai") << " not existed" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // bam file
    for (; optind<argc;)
    {
        bamFiles.emplace_back(argv[optind++]);

        // check the existence of the bam file
        auto f=*bamFiles.rbegin();
        BamReader bamReader;
        if (!bamReader.Open(f))
        {

            cerr << "[PyroTools-ConsensusGraph] error: "
                 << f << " not existed or invalid" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // canonicalize the regions
    if (regionStringsOfInterest.empty())
    {
        BamMultiReader bamReader;
        bamReader.Open(bamFiles);
        RefVector genomeDict = bamReader.GetReferenceData();

        for (auto g : genomeDict)
        {
            regionStringsOfInterest.emplace_back(g.RefName);
        }
    }

    return 0;
}

// interface
int SimpleVariantCallTool::Run(int argc, char *argv[])
{
    // print help message if no arguments provided
    if (argc==2)
    {
        Help();
        exit(EXIT_SUCCESS);
    }

    // parse command line options
    commandLineParser(argc-1,argv+1);

    // run the simple variant call
    if (SimpleVariantCall()!=0)
        return 1;

    return 0;
}

// region sorting function

inline bool isValidAlignment(BamAlignment& al, int minReadLen, int minMapQual, int flag)
{
    // check read length
    if (!GenericBamAlignmentTools::validReadLength(al,minReadLen))
        return false;

    // check map quality
    if (!GenericBamAlignmentTools::validMapQuality(al,minMapQual))
        return false;

    // check the alignment flag
    if (GenericBamAlignmentTools::isFilteredByFlag(al,flag))
        return false;

    return true;
}

// workhorse
int SimpleVariantCallTool::SimpleVariantCall()
{
    // open bam file
    BamMultiReader bamReader;
    bamReader.Open(bamFiles);
    // the dictionary of the chromosomes
    RefVector genomeDict = bamReader.GetReferenceData();

    // define the scanning window
    if (verbose>=1) Verbose(string("define the scanning window"));
    vector<tuple<int,int,int>> scanWindowSet;
    int numWindow = GenericRegionTools::toScanWindow(genomeDict, regionStringsOfInterest, windowSize, scanWindowSet);

    // compute the width for cout
    int wc1=3,wc2=3,wc3=3,wc4=1,wc5=1;
    for_each(genomeDict.begin(),genomeDict.end(),[&](RefData& a){
        if (wc1<a.RefName.length()) wc1 = a.RefName.length();
        if (wc2<to_string(a.RefLength).length()) wc2 = to_string(a.RefLength).length();
        if (wc3<to_string(a.RefLength).length()) wc3 = to_string(a.RefLength).length();
    });

    // lambda funtion for output
    auto repoVar = [&](GenericAllele& allele)
    {
        // to avoid cout interleaving error in parallel mode
        stringstream buffer;
        string info = "snp";
        if (allele.m_allele=="-" && allele.m_reference!="-")
            info = "del";
        if (allele.m_allele!="-" && allele.m_reference=="-")
            info = "ins";
        buffer << setw(wc1) << left << genomeDict[allele.m_chrID].RefName << "\t"
               << setw(wc2) << left << (allele.m_chrPosition+1) << "\t"
               << setw(wc3) << left << (allele.m_chrPosition+1+1) << "\t"
               << setw(wc4) << left << allele.m_reference << "\t"
               << setw(wc5) << left << allele.m_allele << "\t"
               << setw(3)   << left << allele.m_alleleDepth << "\t"
               << setw(3)   << left << allele.m_globalDepth << "\t"
               << setw(3)   << left << info
               << "\n";
        std::cout << buffer.str() << std::flush;
    };

    clock_t allTime = 0;
    // loop over windows
    if (verbose>=1) Verbose(string("looping over the windows"));

    omp_set_dynamic(0);
    omp_set_num_threads(numThreads);
    #pragma omp parallel for shared(genomeDict,repoVar) reduction(+:allTime)
    for (int i=0; i<numWindow; i++)
    {
        BamMultiReader bamReader;
        bamReader.Open(bamFiles);

        clock_t tStart = clock();

        // define an alignment
        BamAlignment aln;

        // alleles in the window
        unordered_map<string,GenericAllele> alleles;

        // define a lambda function
        auto depoVar = [&](int tId, int tPos, char tRef, char tAlt, double tBq, int tDep)
        {
            string tKey = GenericAllele::key(tId, tPos, tRef, tAlt);
            auto tPtr = alleles.find(tKey);
            if (tPtr==alleles.end())
            {
                GenericAllele allele;
                allele.m_chrID = tId;
                allele.m_chrPosition = tPos;
                allele.m_reference = tRef;
                allele.m_globalDepth = tDep;
                allele.m_allele = tAlt;
                allele.m_alleleDepth = 1;
                allele.m_alleleMapAvgQual = aln.MapQuality;
                allele.m_alleleStrandPos = (aln.IsReverseStrand()?0:1);
                allele.m_alleleStrandNeg = (aln.IsReverseStrand()?1:0);
                allele.m_bamData.emplace_back(tuple<char,int,int,int>(tAlt,tBq,aln.MapQuality,((aln.IsReverseStrand()?-1:1))));
                alleles.emplace(make_pair(tKey, allele));
            }else
            {
                tPtr->second.m_alleleDepth += 1;
                tPtr->second.m_alleleMapAvgQual += aln.MapQuality;
                tPtr->second.m_alleleStrandPos += (aln.IsReverseStrand()?0:1);
                tPtr->second.m_alleleStrandNeg += (aln.IsReverseStrand()?1:0);
                tPtr->second.m_bamData.emplace_back(tuple<char,int,int,int>(tAlt,tBq,aln.MapQuality,((aln.IsReverseStrand()?-1:1))));
            }
        };

        // window position
        tuple<int,int,int> w = scanWindowSet[i];
        int wId = get<0>(w);
        int wLp = get<1>(w);
        int wRp = get<2>(w);

        // get sequencing depth
        if (verbose>=1) Verbose("retrieve the depths in "+genomeDict[wId].RefName+":"+to_string(wLp+1)+"-"+to_string(wRp));

        map<long,long> genomeDepth;
        string CMD, DepthResult;
        stringstream ssCMD;
        ssCMD << "samtools depth ";
        ssCMD << "-r " << genomeDict[wId].RefName+":"+to_string(wLp+1)+"-"+to_string(wRp) << " ";
        ssCMD << "-q " << baseQualThres << " ";
        ssCMD << "-Q " << mapQualThres << " ";
        ssCMD << "-l " << readLenThres << " ";
        for (auto bf : bamFiles)
            ssCMD << " " << bf;
        CMD = ssCMD.str();

        // run the command
        RunCmdGetReturn(CMD, DepthResult);
        // parse the result
        {
            string tId;
            int tPos, tDep;
            string result;
            stringstream ssResult;
            stringstream ssDepthResult(DepthResult);
            while(getline(ssDepthResult, result)){
                ssResult << result;
                ssResult >> tId >> tPos >> tDep;
                ssResult.clear();
                ssResult.str("");
                genomeDepth[tPos-1] = tDep;
            }
        }

        // get variants
        if (verbose>=1) Verbose("retrieve the variants in "+genomeDict[wId].RefName+":"+to_string(wLp+1)+"-"+to_string(wRp));

        // rewind the reader
        bamReader.Rewind();
        // set the window region
        bamReader.SetRegion(wId, wLp, wId, wRp);
        // retrieve an alignment in the window
        while (bamReader.GetNextAlignment(aln))
        {
            // skip the alignment if it doesn't overlap the window
            if (aln.Position>=wRp || aln.GetEndPosition()<=wLp)
                continue;

            // skip the invalid alignment
            if (!isValidAlignment(aln, readLenThres, mapQualThres, alnFlagMarker))
                continue;

            // skip the alignment harboring too many mismatches
            if (!GenericBamAlignmentTools::validReadIdentity(aln, 1-alnIdentityThres))
                continue;

            if (!ignoreSnp)
            {
                // retrieve the SNPs on the read
                vector<long> readSnpPos;
                vector<char> readSnpRef;
                vector<char> readSnpAlt;
                vector<double> readSnpBq;
                GenericBamAlignmentTools::getBamAlignmentMismatches(aln, readSnpPos, readSnpRef, readSnpAlt, readSnpBq);
                for (size_t j=0; j<readSnpAlt.size(); j++)
                {
                    if (readSnpBq[j]<baseQualThres) continue;
                    if (readSnpPos[j]<wLp || readSnpPos[j]>=wRp) continue;
                    depoVar(wId, readSnpPos[j], readSnpRef[j], readSnpAlt[j], readSnpBq[j], genomeDepth[readSnpPos[j]]);
                }
            }

            if (!ignoreIndel)
            {
                // retrieve the inserts on the read
                vector<long> readInsPos;
                vector<string> readInsSeq;
                vector<vector<double>> readInsBq;
                GenericBamAlignmentTools::getBamAlignmentInserts(aln, readInsPos, readInsSeq, readInsBq);
                for (size_t j=0; j<readInsPos.size(); j++)
                {
                    int ip = readInsPos[j];
                    string iseq = readInsSeq[j];
                    vector<double> ibq = readInsBq[j];
                    for (size_t k=0; k<iseq.length(); k++, ip++)
                    {
                        if (ip<wLp || ip>=wRp) continue;
                        depoVar(wId, ip, '-', iseq[k], ibq[k], genomeDepth[ip]);
                    }
                }
            }

            if (!ignoreIndel)
            {
                // retrieve the deletes on the read
                vector<long> readDelPos;
                vector<string> readDelSeq;
                GenericBamAlignmentTools::getBamAlignmentDeletes(aln, readDelPos, readDelSeq);
                for (size_t j=0; j<readDelPos.size(); j++)
                {
                    int dp = readDelPos[j];
                    string dseq = readDelSeq[j];
                    for (size_t k=0; k<dseq.length(); k++, dp++)
                    {
                        if (dp<wLp || dp>=wRp) continue;
                        depoVar(wId, dp, dseq[k], '-', 0, genomeDepth[dp]);
                    }
                }
            }
        }

        if (verbose>=1) Verbose("output the variants");
        // filter and output
        for (auto ptr : alleles)
        {
            GenericAllele allele = get<1>(ptr);

            if (ignoreSnp && allele.m_reference!="-" && allele.m_allele!="-")
                continue;
            if (ignoreIndel && (allele.m_reference=="-" || allele.m_allele=="-"))
                continue;

            if (allele.m_reference!="-" && allele.m_allele!="-")
            {
                if (allele.m_alleleDepth>snpHitThres)
                {
                    repoVar(allele);
                }
            }else
            {
                if (allele.m_alleleDepth>indelHitThres)
                {
                    repoVar(allele);
                }
            }
        }

        clock_t tEnd = clock();
        allTime += tEnd - tStart;
        if (verbose>=1) Verbose("time elapsed "+to_string((double)(tEnd-tStart)/CLOCKS_PER_SEC)+" seconds");
    }

    if (verbose>=1) Verbose("total time elapsed "+to_string((double)allTime/CLOCKS_PER_SEC)+" seconds");

    bamReader.Close();

    return 0;
}
