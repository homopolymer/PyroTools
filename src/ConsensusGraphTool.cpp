#include "GenericGraphTools.h"
#include "GenericRegionTools.h"
#include "GenericBamAlignmentTools.h"
#include "GenericVariant.h"
#include "GenericVcfTools.h"
using namespace GenericSequenceTools;

#include <sys/stat.h>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <forward_list>
#include <utility>
#include <string.h>
#include <algorithm>
#include <random>
#include <ctime>
using namespace std;

#include <api/BamAux.h>
#include <api/BamReader.h>
#include <utils/bamtools_fasta.h>
using namespace BamTools;

inline void Verbose(string msg)
{
    cerr << "[" << CurrentTime() << " PyroTools-ConsensusGraph] " << msg << endl;
}

enum ConsensusGraphOutputFormat{BED,FASTA,DOT};

ConsensusGraphTool::ConsensusGraphTool()
    :HasGenomeFile(false)
    ,HasBamAlnFile(false)
    ,HasRegion(false)
    ,topK(5)
    ,topK2(30)
    ,edgePruneLevel(1)
    ,minReadLength(100)
    ,minMapQuality(5)
    ,filterFlag(0)
    ,edgePruneFrac(0.1)
    ,minUniqHit(1)
    ,useDiffRead(false)
    ,useMarkovian(false)
    ,scoreMethod("likeli")
    ,skipIndel(true)
    ,verbose(0)
    ,OutputFormat(BED)
{
    ;
}

int ConsensusGraphTool::Help()
{
    cerr << "SYNOPSIS" << endl;
    cerr << "    PyroTools GraphConsensus [OPTIONS] <GENOME> <BAM> <REGION>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    GraphConsensus outputs the top-k paths in the consensus graph" << endl;
    cerr << "    built from the multiple alignments in the BAM file." << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "--options for graph searching" << endl;
    cerr << "    -s,--score         specify the scoring method, valid arguments: likeli and assign [STR]" << endl;
    cerr << "                           \"likeli\" is for the single sample data (default)" << endl;
    cerr << "                           \"assign\" is for the population data" << endl;
    cerr << "    -m,--markov        use the first order markovian probability between nearby mismatches (experimental)" << endl;
    cerr << "    -k,--top-k         the number of graph paths being output, all paths will be output" << endl;
    cerr << "                       if the number of paths is less than the value (default:5) [INT]" << endl;
    cerr << "    -K,--top-K         the number of graph paths being counted internally (default:30) [INT]" << endl;
    cerr << "--options for graph construction" << endl;
    cerr << "    -e,--edge-prune    graph pruning based on the edge hit counts (default:1) [INT]" << endl;
    cerr << "    -c,--edge-frac     graph pruning based on the edge fraction (default:0.1) [FLT]" << endl;
    cerr << "    -d,--diff-read     use only the reads harboring differences" << endl;
    cerr << "    -n,--uniq-hit      minimal hit number of unique read (default:1) [INT]" << endl;
    cerr << "    -l,--len           minimal read length (default:100) [INT]" << endl;
    cerr << "    -q,--mq            minimal mapping quality (default:5) [INT]" << endl;
    cerr << "    --ff               filter reads with the specific flag(s) [INT]" << endl;
    cerr << "    --skip-read        skip the reads named in the file [STR]" << endl;
    cerr << "--options for variant information" << endl;
    cerr << "    --vcf              filename to save the variants in the top-k paths [STR]" << endl;
    cerr << "    -I,--indel         inclusively print indel" << endl;
    cerr << "--options for output" << endl;
    cerr << "    -f,--out-format    set the file format of output, valid arguments: bed, fasta, dot" << endl;
    cerr << "                       (default:bed) [STR]" << endl;
    cerr << "    -v,--verbose       the verbose level (default:0) [INT]" << endl;
    cerr << "    -h,--help          print help message" << endl;
    return 0;
}

// check the existence of a file
// reference: http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
inline bool file_exist(const string& filename)
{
    struct stat buffer;
    return (stat (filename.c_str(), &buffer) == 0);
}

// reference
// http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html
int ConsensusGraphTool::commandOptionsParser(int argc, char *argv[])
{
    ifstream infile;
    string line;

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
            {"top-k",      required_argument, 0, 'k'},
            {"edge-prune", required_argument, 0, 'e'},
            {"out-format", required_argument, 0, 'f'},
            {"len",        required_argument, 0, 'l'},
            {"mq",         required_argument, 0, 'q'},
            {"ff",         required_argument, 0, 0},
            {"skip-read",  required_argument, 0, 0},
            {"vcf",        required_argument, 0, 0},
            {"edge-frac",  required_argument, 0, 'c'},
            {"uniq-hit",   required_argument, 0, 'n'},
            {"top-K",      required_argument, 0, 'K'},
            {"diff-read",  no_argument,       0, 'd'},
            {"score",      required_argument, 0, 's'},
            {"markov",     no_argument,       0, 'm'},
            {"indel",      no_argument,       0, 'I'},
            {"verbose",    required_argument, 0, 'v'},
            {"help",       no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "k:e:f:l:q:c:n:K:s:v:mIdh", long_options, &option_index);

        // detect the end of the options
        if (c==-1) break;

        switch(c)
        {
        case 0:

            switch(option_index)
            {
            case 5:
                if (isHexString(optarg)){
                    filterFlag |= stoul(optarg,nullptr,16);
                }else{
                    filterFlag |= stoi(optarg);
                }
                break;
            case 6:
                infile.open(optarg);
                while (getline(infile,line))
                {
                    if (!line.empty())
                    {
                        skipReadList.emplace(line);
                    }
                }
                infile.close();
                break;
            case 7:
                vcfOut = optarg;
                break;
            }

            break;

        case 'k':
            topK=stoi(optarg);
            break;

        case 'K':
            topK2=stoi(optarg);
            break;

        case 'e':
            edgePruneLevel=stoi(optarg);
            break;

        case 'f':
            if (strcmp(optarg,"fasta")==0) OutputFormat=FASTA;
            if (strcmp(optarg,"dot")==0) OutputFormat=DOT;
            break;

        case 'l':
            minReadLength=stoi(optarg);
            break;

        case 'q':
            minMapQuality=stoi(optarg);
            break;

        case 'c':
            edgePruneFrac=stod(optarg);
            break;

        case 'n':
            minUniqHit=stoi(optarg);
            break;

        case 'd':
            useDiffRead=true;
            break;

        case 's':
            scoreMethod=optarg;
            break;

        case 'm':
            useMarkovian = true;
            break;

        case 'I':
            skipIndel = false;
            break;

        case 'v':
            verbose=stoi(optarg);
            break;

        case 'h':
            Help();
            exit(EXIT_SUCCESS);
            break;

        default:
            abort();
        }
    }

    // genome file
    if (optind<argc)
    {
        HasGenomeFile = true;
        GenomeFile=argv[optind++];

        // check the existence of the genome file
        if (!file_exist(GenomeFile))
        {
            cerr << "[PyroTools-ConsensusGraph] error: "
                   << GenomeFile << " not existed" << endl;
            exit(EXIT_FAILURE);
        }

        // check the existence of the genome index file
        if (!file_exist(GenomeFile+".fai"))
        {
            cerr << "[PyroTools-ConsensusGraph] error: "
                   << (GenomeFile+".fai") << " not existed" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // BAM file
    if (optind<argc)
    {
        HasBamAlnFile = true;
        BamAlnFile = argv[optind++];

        // check the existence of the bam file
        BamReader bamReader;
        if (!bamReader.Open(BamAlnFile))
        {
            cerr << "[PyroTools-ConsensusGraph] error: "
                 << BamAlnFile << " not existed or invalid" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // Region
    if (optind<argc)
    {
        HasRegion = true;
        RegionOfInterest = argv[optind++];
    }
    // check the existence of the region
    if (RegionOfInterest.empty())
    {
        cerr << "[PyroTools-ConsensusGraph] error: "
             << "the region is not specified" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

int ConsensusGraphTool::Run(int argc, char *argv[])
{
    // print help message if no arguments provided
    if (argc==2)
    {
        Help();
        exit(EXIT_SUCCESS);
    }

    // parse the command line
    commandOptionsParser(argc-1, argv+1);

    // executation
    if (ConsensusSequence()==EXIT_FAILURE)
        return EXIT_FAILURE;

    return 0;
}


// temporary struct
struct sequence_t
{
    string          m_ID;
    int             m_startPositionShift;
    int             m_endPositionShift;

    string          m_sequence;
    Cigar           m_cigar;
    int             m_type;

    void key(string& s)
    {
        stringstream ss;
        ss << m_sequence;

        for (auto ptr=m_cigar.begin(); ptr!=m_cigar.end(); ptr++)
        {
            ss << ptr->Length << ptr->Type;
        }

        ss << m_startPositionShift;

        s = ss.str();
    }
};

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

int ConsensusGraphTool::ConsensusSequence()
{
    // open BAM file
    BamReader bamReader;
    bamReader.Open(BamAlnFile);

    // genome dictionary
    RefVector genomeDict;
    genomeDict = bamReader.GetReferenceData();

    // parse the region
    string roiChrName;
    int roiChrLeftPos, roiChrRightPos;
    GenericRegionTools::parse(RegionOfInterest, roiChrName, roiChrLeftPos, roiChrRightPos);
    if (roiChrRightPos<0)
    {
        for (auto g : genomeDict)
        {
            if (roiChrName==g.RefName)
                roiChrRightPos=g.RefLength;
        }
    }

    // retrieve fasta sequence
    Fasta genomeFasta;
    genomeFasta.Open(GenomeFile, GenomeFile+".fai");

    int roiChrID=0;
    for (auto genomeItem=genomeDict.begin(); genomeItem!=genomeDict.end(); genomeItem++, roiChrID++)
    {
        if (genomeItem->RefName == roiChrName)
        {
            break;
        }
    }

    string genomeSequence;
    genomeFasta.GetSequence(roiChrID, roiChrLeftPos, roiChrRightPos-1, genomeSequence);


    // set BAM region
    bamReader.SetRegion(roiChrID, roiChrLeftPos, roiChrID, roiChrRightPos);

    if (verbose>=1)
        Verbose("retrieve the alignments");

    // retrieve BAM alignment
    int numTotAln = 0;
    int numAln = 0;
    BamAlignment alnObj;
    unordered_set<string> alnObjPool;
    vector<sequence_t> seqObjPool;
    while(bamReader.GetNextAlignment(alnObj))
    {
        // skip alignment out of ROI
        if (alnObj.Position>=roiChrRightPos || alnObj.GetEndPosition()<=roiChrLeftPos)
            continue;

        numTotAln += 1;

        if (!isValidAlignment(alnObj, minReadLength, minMapQuality, filterFlag))
            continue;

        // alignment name
        string alnName = GenericBamAlignmentTools::getBamAlignmentName(alnObj);
        auto ptr = skipReadList.find(alnName);
        if (ptr!=skipReadList.end())
            continue;

        // alignment id
        string alnID = GenericBamAlignmentTools::getBamAlignmentName(alnObj);

        // continue if read already meet
        if (alnObjPool.find(alnID)!=alnObjPool.end())
            continue;
        alnObjPool.emplace(alnID);

        // extract local sequences
        string localRead, localGenome;
        Cigar  cigar;
        int numDiff = GenericBamAlignmentTools::getLocalAlignment(alnObj, roiChrLeftPos, roiChrRightPos-roiChrLeftPos,
                                                                  localRead, localGenome, cigar);

        if (localRead.empty() || localGenome.empty())
            continue;

        // make sequence object
        sequence_t seqObj;
        seqObj.m_ID           = GenericBamAlignmentTools::getBamAlignmentName(alnObj);
        seqObj.m_sequence     = localRead;
        seqObj.m_cigar        = cigar;

        // skip if sequence is too short
        if (seqObj.m_sequence.length()*10<(roiChrRightPos-roiChrLeftPos) &&
            seqObj.m_sequence.length()*10<GenericBamAlignmentTools::getBamAlignmentReadLength(alnObj))
            continue;

        if (alnObj.Position>roiChrLeftPos)
            seqObj.m_startPositionShift = alnObj.Position-roiChrLeftPos;
        else
            seqObj.m_startPositionShift = 0;

        if (alnObj.GetEndPosition()<roiChrRightPos)
            seqObj.m_endPositionShift = roiChrRightPos-alnObj.GetEndPosition();
        else
            seqObj.m_endPositionShift = 0;

        seqObj.m_type = READ_CROSS_REGION;
        if (seqObj.m_startPositionShift>0)
            seqObj.m_type = READ_BEGIN_REGION;
        if (seqObj.m_endPositionShift>0)
            seqObj.m_type = READ_END_REGION;

        if (numDiff==0 && useDiffRead)
            continue;

        seqObjPool.emplace_back(seqObj);

        // read counting
        numAln ++;
    }
    // sort by sequence length
    sort(seqObjPool.begin(), seqObjPool.end(), [](const sequence_t& a, const sequence_t &b)->bool{
       return a.m_sequence.length()>b.m_sequence.length();
    });

    if (verbose>=1)
        Verbose("make the unique sequences");

    // make unique sequences
    vector<string> uniqueSeqObjSeqPool, uniqueSeqObjSeqPool2;
    vector<Cigar> uniqueSeqObjCigarPool, uniqueSeqObjCigarPool2;
    vector<int> uniqueSeqObjStartPosPool, uniqueSeqObjStartPosPool2;
    vector<int> uniqueSeqObjHitCountPool, uniqueSeqObjHitCountPool2;
    unordered_map<string,int> uniqueSeqObjIndex;
    int idx=0;
    for (auto ptr=seqObjPool.begin(); ptr!=seqObjPool.end(); ptr++)
    {
        bool isduplicate = false;
        string key;

        ptr->key(key);
        auto ptr2 = uniqueSeqObjIndex.find(key);
        if (ptr2!=uniqueSeqObjIndex.end())
        {
            uniqueSeqObjHitCountPool[ptr2->second] += 1;
            isduplicate = true;
        }

        if (isduplicate)
            continue;

        uniqueSeqObjSeqPool.emplace_back(ptr->m_sequence);
        uniqueSeqObjCigarPool.emplace_back(ptr->m_cigar);
        uniqueSeqObjStartPosPool.emplace_back(ptr->m_startPositionShift);
        uniqueSeqObjIndex.emplace(make_pair(key,idx));
        uniqueSeqObjHitCountPool.emplace_back(1);

        idx++;
    }

    numAln = 0;
    for (int i=0; i<uniqueSeqObjHitCountPool.size(); i++)
    {
        if (uniqueSeqObjHitCountPool[i]>=minUniqHit)
        {
            uniqueSeqObjSeqPool2.emplace_back(uniqueSeqObjSeqPool[i]);
            uniqueSeqObjCigarPool2.emplace_back(uniqueSeqObjCigarPool[i]);
            uniqueSeqObjStartPosPool2.emplace_back(uniqueSeqObjStartPosPool[i]);
            uniqueSeqObjHitCountPool2.emplace_back(uniqueSeqObjHitCountPool[i]);
            numAln += uniqueSeqObjHitCountPool[i];
        }
    }

    vector<tuple<string,Cigar,int,int>> uniqueSeqPool;
    for (int i=0; i<uniqueSeqObjSeqPool2.size(); i++)
    {
        uniqueSeqPool.emplace_back(tuple<string,Cigar,int,int>(uniqueSeqObjSeqPool2[i],uniqueSeqObjCigarPool2[i],uniqueSeqObjStartPosPool2[i],uniqueSeqObjHitCountPool2[i]));
    }
    sort(uniqueSeqPool.begin(), uniqueSeqPool.end(), [](const tuple<string,Cigar,int,int> &a, const tuple<string,Cigar,int,int> &b)->bool{return get<2>(a)<get<2>(b);});
    for (int i=0; i<uniqueSeqPool.size(); i++)
    {
        uniqueSeqObjSeqPool2[i] = get<0>(uniqueSeqPool[i]);
        uniqueSeqObjCigarPool2[i] = get<1>(uniqueSeqPool[i]);
        uniqueSeqObjStartPosPool2[i] = get<2>(uniqueSeqPool[i]);
        uniqueSeqObjHitCountPool2[i] = get<3>(uniqueSeqPool[i]);
    }

    if (verbose>=1)
        Verbose("construct the graph");

    clock_t t1 = clock();
    // construct the graph
    GenericDagGraph consensusGraph;
    consensusGraph.buildDagGraph(genomeSequence, uniqueSeqObjSeqPool2, uniqueSeqObjCigarPool2, uniqueSeqObjStartPosPool2, uniqueSeqObjHitCountPool2);
    if (edgePruneLevel>0 || edgePruneFrac>0)
        consensusGraph.edgePruningMix(edgePruneLevel, edgePruneFrac, true);

    if (OutputFormat==DOT)
    {
        if (verbose>=1)
            Verbose("output the graph");

        cout << consensusGraph << endl;
    }else
    {
        if (verbose>=1)
            Verbose(string("search the top-" + to_string(topK) + " paths").c_str());

        vector<string>       topRankGraphPaths;
        vector<list<Vertex>> topRankGraphPathVertexs;
        vector<double>       topRankGraphPathWeights;
        consensusGraph.topRankPathsByErrMdl(topK, topK2, scoreMethod, useMarkovian,topRankGraphPaths, topRankGraphPathVertexs, topRankGraphPathWeights);

        int numPath = topK;
        if (numPath > topRankGraphPaths.size())
            numPath = topRankGraphPaths.size();

        if (verbose>=1)
            Verbose(string("output the top-" + to_string(numPath) + " paths").c_str());

        for (int i=0; i<numPath; i++)
        {
            Cigar pathCigar;
            consensusGraph.pathCigar(topRankGraphPathVertexs[i], pathCigar);

            string pathCigarString;
            GenericBamAlignmentTools::convertCigarToString(pathCigar, pathCigarString);

            if (OutputFormat==FASTA)
            {
                cout << ">GraphCns" << i << "\t"
                     << topRankGraphPathWeights[i] << endl;
                cout << topRankGraphPaths[i] << endl;
            }else
            {
                cout << roiChrName << "\t" << roiChrLeftPos << "\t" << roiChrRightPos << "\t"
                     << "GraphCns" << i << "\t"
                     << topRankGraphPaths[i] << "\t" << pathCigarString << "\t" << topRankGraphPathWeights[i] << endl;
            }
        }

        // output vcf file
        if (!vcfOut.empty()){
            vector<string> pathname;
            // collect variants
            map<int, GenericVariant> variants;
            for (int i=0; i<numPath; i++){
                string p = "GraphCns" + to_string(i);
                pathname.emplace_back(p);
                auto lvs = topRankGraphPathVertexs[i];
                map<int, GenericAllele> pathVariant;
                consensusGraph.pathVariantCollect(lvs, pathVariant);

                for (auto v : pathVariant){
                    // skip indel
                    if (skipIndel && v.second.m_reference.length()<v.second.m_allele.length())
                        continue;
                    if (skipIndel && v.second.m_allele=="-")
                        continue;

                    // new variant
                    if (variants.count(v.first)==0){
                        GenericVariant variant;
                        variant.m_chrID = roiChrID;
                        variant.m_chrPosition = roiChrLeftPos + v.second.m_chrPosition;
                        variant.m_reference = v.second.m_reference;
                        variant.AlleleInsert(v.second);
                        int idx = variant.AlleleIndex(v.second.m_allele);
                        variant.m_samples[p] = GenericVariant::PhasedAlleles(1,idx);
                        variants[v.first] = variant;
                    }
                    // existing variant
                    else{
                        int idx = variants[v.first].AlleleIndex(v.second.m_allele);
                        if (idx==-1) variants[v.first].AlleleInsert(v.second);
                        idx = variants[v.first].AlleleIndex(v.second.m_allele);
                        variants[v.first].m_samples[p] = GenericVariant::PhasedAlleles(1,idx);
                    }
                }
            }

            // output
            GenericVcfTools::write(genomeDict, vcfOut, pathname, variants);
        }
    }

    clock_t t2 = clock();


    // print summary message
    if (verbose>=1)
        Verbose("process " + to_string(numAln) + " of " + to_string(numTotAln) + " reads in " + RegionOfInterest);

    if (verbose>=1)
        Verbose("time elapsed " + to_string(double(t2-t1)/CLOCKS_PER_SEC) + " in graph analysis");

    // close genome file and bam file
    bamReader.Close();
    genomeFasta.Close();

    return 0;
}
