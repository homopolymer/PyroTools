#include "GenericBamAlignmentTools.h"
#include "GenericRegionTools.h"
using namespace GenericSequenceTools;

#include <getopt.h>
#include <omp.h>
#include <time.h>

#include <vector>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <list>
#include <unordered_set>
using namespace std;

#include "api/BamReader.h"
#include "api/BamMultiReader.h"
#include "api/BamAux.h"
#include "api/BamAlignment.h"
using namespace BamTools;

CropBamTool::CropBamTool()
    :mapQualThres(5)
    ,readLenThres(100)
    ,segmentLenThres(30)
    ,alnFlagMarker(0)
    ,alnIdenThres(0.9)
    ,outFormat("fasta")
    ,keepClip(false)
    ,useUnique(false)
    ,thresFreq(1)
    ,numThreads(1)
    ,verbose(0)
{

}

inline void Verbose(string MSG)
{
    auto display_time = [](){
        time_t now = time(nullptr);
        char *ts = asctime(localtime(&now));
        ts[strlen(ts)-1] = '\0';
        return ts;
    };


    cerr << "[" << display_time() << " PyroTools-CropBam] " << MSG << endl;
}

// help message
int CropBamTool::Help()
{
    cerr << "SYNOPSIS" << endl;
    cerr << "    PyroTools CropBam [OPTIONS] <BAM> [<BAM>...]" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    Extract reads in a given regions, and chop off the read segment outside the region" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "    --format        specify the output format, valid argument: fasta or fastq (default:fasta) [STR]" << endl;
    cerr << "    --roi           specify the region of interest [STR]" << endl;
    cerr << "    --roi-list      specify the region(s) of interest" << endl;
    cerr << "    --mq            skip reads with mapping quality less than the value (default:5) [INT]" << endl;
    cerr << "    --len           skip reads with length less than the value (default:100) [INT]" << endl;
    cerr << "    --slen          skip the incomplete read segment with length less than the value (default:30) [INT]" << endl;
    cerr << "    --iden          skip reads with identity less than the value (default:0.9) [FLT]" << endl;
    cerr << "    --ff            skip reads with the specified flags [INT]" << endl;
    cerr << "    --keep-clip     keep the clip sequence" << endl;
    cerr << "    --uniq          output the unique reads" << endl;
    cerr << "    --freq          specify the file saving the frequency of unique reads [STR]" << endl;
    cerr << "    --freq-thres    specify the frequency threshold of unique reads (default:1) [INT]" << endl;
    cerr << "    --num-cores     specify the number of threads (default:1) [INT]" << endl;
    cerr << "    -v,--verbose    specify the level of verbosity, valid arguments: 0 for quiet, 1 for debug (default:0) [INT]" << endl;
    cerr << "    -h,--help       print help message" << endl;
    return 0;
}

int CropBamTool::parseCommandLine(int argc, char *argv[])
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
            {"slen",          required_argument, 0, 0},
            {"iden",          required_argument, 0, 0},
            {"ff",            required_argument, 0, 0},

            {"num-cores",     required_argument, 0, 0},
            {"format",        required_argument, 0, 0},
            {"keep-clip",     required_argument, 0, 0},
            {"uniq",          no_argument,       0, 0},
            {"freq",          required_argument, 0, 0},
            {"freq-thres",    required_argument, 0, 0},
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
                regionStrings.emplace_back(optarg);
                break;
            case 1:
                {
                    ifstream infile(optarg);
                    string line;
                    while (getline(infile,line))
                    {
                        if (!line.empty())
                            regionStrings.emplace_back(optarg);
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
                segmentLenThres=stoi(optarg);
                break;
            case 5:
                alnIdenThres=stod(optarg);
                break;
            case 6:
                if (isHexString(optarg)){
                    alnFlagMarker |= stoul(optarg,nullptr,16);
                }else{
                    alnFlagMarker |= stoi(optarg);
                }
                break;
            case 7:
                numThreads=stoi(optarg);
                break;
            case 8:
                outFormat=optarg;
                break;
            case 9:
                keepClip=true;
                break;
            case 10:
                useUnique = true;
                break;
            case 11:
                outFreq = optarg;
                break;
            case 12:
                thresFreq = stoi(optarg);
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

    // bam file
    for (; optind<argc;)
    {
        bamFiles.emplace_back(argv[optind++]);

        // check the existence of the bam file
        auto f=*bamFiles.rbegin();
        BamReader bamReader;
        if (!bamReader.Open(f))
        {

            cerr << "[PyroTools-CropBam] error: "
                 << f << " not existed or invalid" << endl;
            exit(EXIT_FAILURE);
        }
    }
    return 0;

}

// interface
int CropBamTool::Run(int argc, char *argv[])
{
    // print help message if no arguments provided
    if (argc==2)
    {
        Help();
        exit(0);
    }

    // parse the command line
    parseCommandLine(argc-1,argv+1);

    // run the program
    if (CropBam()!=0)
        return EXIT_FAILURE;

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


int CropBamTool::CropBam()
{
    // open bam files
    BamMultiReader bamReader;
    bamReader.Open(bamFiles);

    // the dictionary of chromosomes
    RefVector genome = bamReader.GetReferenceData();

    // get the scanning window
    vector<tuple<int,int,int>> windows;
    int numWindows = GenericRegionTools::toScanWindow(genome, regionStrings, windows);

    unordered_set<string> readpool;

    // temporary struct for sequence object
    typedef struct {
        string name;
        int head_soft_clip;
        int tail_soft_clip;
        string seq;
        string qual;
    }cropbam_seq_t;

    // temporary struct for unique seqs
    map<string,list<cropbam_seq_t>> uniqueSeqPool;

    // lambda expression for output
    auto Output = [this](cropbam_seq_t &a){
          if (this->outFormat=="fasta"){
              cout << ">" << a.name << "\t"
                   << "head_soft_clip=" << a.head_soft_clip << "\t"
                   << "tail_soft_clip=" << a.tail_soft_clip << "\t"
                   << endl
                   << a.seq << endl;
          }
          if (this->outFormat=="fastq"){
              cout << "@" << a.name << "\t"
                   << "head_soft_clip=" << a.head_soft_clip << "\t"
                   << "tail_soft_clip=" << a.tail_soft_clip << "\t"
                   << endl
                   << a.seq << endl;
              cout << "+" << endl
                   << a.qual << endl;
          }
    };

    // loop over windows
    omp_set_dynamic(0);
    omp_set_num_threads(numThreads);
    #pragma omp parallel for shared(genome)
    for (int i=0; i<numWindows; i++)
    {
        clock_t tStart = clock();

        bamReader.Open(bamFiles);

        int wId = get<0>(windows[i]);
        int wLp = get<1>(windows[i]);
        int wRp = get<2>(windows[i]);

        if (verbose>=1) Verbose("process the window " + genome[wId].RefName + ":" + to_string(wLp+1) + "-" + to_string(wRp));

        // rewind the bam reader
        bamReader.Rewind();
        // set the region
        bamReader.SetRegion(wId, wLp, wId, wRp);

        int numReads = 0;
        // retrieve the alignment
        BamAlignment aln;
        while (bamReader.GetNextAlignment(aln))
        {
            // skip the alignment if it doesn't overlap the window
            if (aln.Position>=wRp || aln.GetEndPosition()<=wLp)
                continue;

            // skip the invalid alignment
            if (!isValidAlignment(aln, readLenThres, mapQualThres, alnFlagMarker))
                continue;

            // skip the alignment harboring too many mismatches
            if (!GenericBamAlignmentTools::validReadIdentity(aln, 1-alnIdenThres))
                continue;

            stringstream keyss;
            keyss << GenericBamAlignmentTools::getBamAlignmentName(aln) << "-"
                  << wId << "-" << wLp << "-" << wRp;
            string key = keyss.str();
            auto ptr = readpool.find(key);
            if (ptr!=readpool.end())
                continue;
            readpool.emplace(key);

            // get the partial read
            string readSegment, readQualSegment, genomeSegment;
            GenericBamAlignmentTools::getLocalAlignment(aln, wLp, wRp-wLp, readSegment, readQualSegment, genomeSegment);

            // add soft clip
            int hsc=0;
            auto ptr0 = aln.CigarData.begin();
            if (aln.Position>=wLp && (ptr0->Type=='S' || ptr0->Type=='H'))
            {
                stringstream headClipSeq, headClipQual;
                for (int i=0; i<ptr0->Length; i++)
                {
                    headClipSeq << aln.QueryBases[i];
                    headClipQual << aln.Qualities[i];
                }

                if (keepClip)
                {
                    readSegment=headClipSeq.str()+readSegment;
                    readQualSegment=headClipQual.str()+readQualSegment;
                }

                hsc += ptr0->Length;
            }
            int tsc=0;
            auto ptr1 = aln.CigarData.rbegin();
            if (aln.GetEndPosition()<wRp && (ptr1->Type=='S' || ptr1->Type=='H'))
            {
                string ss="", qs="";
                auto str=aln.QueryBases.rbegin();
                auto qtr=aln.Qualities.rbegin();
                for (int i=0; i<ptr1->Length; i++,str++,qtr++)
                {
                    ss=(*str)+ss;
                    qs=(*qtr)+qs;
                }
                if (keepClip)
                {
                    readSegment=readSegment+ss;
                    readQualSegment=readQualSegment+qs;
                }
                tsc += ptr1->Length;
            }

            if (readSegment.length()>=segmentLenThres)
            {
                cropbam_seq_t a;
                a.name = GenericBamAlignmentTools::getBamAlignmentName(aln);
                a.head_soft_clip = hsc;
                a.tail_soft_clip = tsc;
                a.seq = readSegment;
                a.qual = readQualSegment;
                if (uniqueSeqPool.count(a.seq)==0)
                    uniqueSeqPool[a.seq] = list<cropbam_seq_t>(1,a);
                else
                    uniqueSeqPool[a.seq].emplace_back(a);
//                if (outFormat=="fasta"){
//                    cout << ">" << GenericBamAlignmentTools::getBamAlignmentName(aln) << "\t"
//                         << "head_soft_clip=" << hsc << "\t"
//                         << "tail_soft_clip=" << tsc << "\t"
//                         << endl
//                         << readSegment << endl;
//                }

//                if (outFormat=="fastq"){
//                    cout << "@" << GenericBamAlignmentTools::getBamAlignmentName(aln) << "\t"
//                         << "head_soft_clip=" << hsc << "\t"
//                         << "tail_soft_clip=" << tsc << "\t"
//                         << endl
//                         << readSegment << endl;
//                    cout << "+" << endl
//                         << readQualSegment << endl;
//                }

                numReads++;
            }
        }

        numReads = 0;
        if (useUnique){
            ofstream of;
            of.open(outFreq);
            for (auto a : uniqueSeqPool){
                if (a.second.size()>=thresFreq){
                    Output(*a.second.begin());
                    of << a.second.begin()->name << "\t" << a.second.size() << endl;
                    numReads ++;
                }
            }
            of.close();
        }else{
            for (auto a : uniqueSeqPool){
                for (auto b : a.second){
                    Output(b);
                    numReads ++;
                }
            }
        }

        clock_t tEnd = clock();

        if (verbose>=1) Verbose("retrieve " + to_string(numReads) + " reads");
        if (verbose>=1) Verbose("time elapsed " + to_string((double)(tEnd-tStart)/CLOCKS_PER_SEC) + " seconds");
    }

    return 0;
}


