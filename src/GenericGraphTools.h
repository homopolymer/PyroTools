#ifndef GENERICGRAPHTOOLS_H
#define GENERICGRAPHTOOLS_H

#include <map>
#include <list>
#include <tuple>
#include <vector>
#include <iostream>
#include <iomanip>
#include <unordered_set>
using namespace std;

#include "GenericVariant.h"
#include "GenericSequenceGlobal.h"
#include "GenericBamAlignmentTools.h"
#include "GenericProbabilisticAlignment.h"
using namespace GenericSequenceTools;

namespace GenericSequenceTools
{

#define READ_CROSS_REGION 0
#define READ_BEGIN_REGION 1
#define READ_END_REGION   2

typedef int Vertex;

struct EdgeVertex
{
    int    m_id;
    double m_weight;
    int    m_count;

    bool operator==(EdgeVertex& b)
    {
        return (m_id==b.m_id);
    }
};

#define DAG_BEGIN_VERTEX "begin"
#define DAG_END_VERTEX   "end"

class GenericGraphRead{
    // ctor & dtor
    public:
        GenericGraphRead()
            :startPos(-1)
            ,count(0)
        {}

    // member
    public:
        string         name;        // read name
        int            startPos;    // the start position on the graph backbone
        list<Vertex>   vertexList;  // the graph representation of the read
        int            count;       // the occurrence number of the read in the data
};

class GenericDagGraph
{
    // ctor & dtor
    public:
        GenericDagGraph()
        {
            m_numVertexs  = 0;
            m_numEdges    = 0;
            m_numRead     = 0;
            m_numMismatch = 0;
        }

        ~GenericDagGraph()
        {
            ;
        }

    // function member
    public:
        //-------------------------------------------------------------------//
        // utilities for the graph encoding of a multiple sequence alignment //
        //-------------------------------------------------------------------//
        // construct the graph by giving a reference and a set of alignments
        void buildDagGraph(string& genome, vector<string>& reads, vector<Cigar>& cigars, vector<int>& starts, vector<int>& counts=*static_cast<vector<int>*>(0));

        // build the graph backbone from a given reference sequence
        void buildDagGraphBackbone(string& genome);

        // update the graph by giving an alignment
        void updateDagGraphByRead(string& read, Cigar& cigar, int& start, int& count=*static_cast<int*>(0), bool sorting=false);

        // find the equivalent insert vertices and merge them
        // as traversing the graph in backward direction
        void backwardRefineDagGraph();

        // insert edge
        void insertEdge(Vertex& prev, Vertex& curr, int& count=*static_cast<int*>(0));

        // insert vertex
        void insertVertex(Vertex& curr, string& label, bool ong, Vertex& prev);

        // topological sorting
        void topologicalSorting();

        // compute succede homopolymer length
        void calculateSuccedeHomopolymer();


        //------------------------------//
        // utilities for graph decoding //
        //------------------------------//


        //-------------------------------//
        // utilities for graph searching //
        //-------------------------------//
        // search top k path
        void topRankPaths(int k, vector<string>& pathLabels,
                          vector<list<Vertex>>& pathVertexs=*static_cast<vector<list<Vertex>>*>(0),
                          vector<double>& pathWeights=*static_cast<vector<double>*>(0));

        void topRankPathsExcludeGenome(int k, vector<string>& pathLabels,
                                       vector<list<Vertex>>& pathVertexs=*static_cast<vector<list<Vertex>>*>(0),
                                       vector<double>& pathWeights=*static_cast<vector<double>*>(0));

        // new version of ranking the paths, scoring the path by the average alignment identity
        void topRankPathsByIden(int k, double quantile, vector<string>& pathLabels,
                                vector<list<Vertex>>& pathVertexs=*static_cast<vector<list<Vertex>>*>(0),
                                vector<double>& pathWeights=*static_cast<vector<double>*>(0));
        // new version of ranking the paths, scoring the path by the average alignment identity
        void topRankPathsByErrMdl(int k, int k2, string scoring, bool markovian, vector<string>& pathLabels,
                                vector<list<Vertex>>& pathVertexs=*static_cast<vector<list<Vertex>>*>(0),
                                vector<double>& pathWeights=*static_cast<vector<double>*>(0));

        // consistency evaluation
        bool isConsistentGraph(int k=30);


        //----------------------------------//
        // utilities for graph manipulation //
        //----------------------------------//
        // convert path to cigar
        void pathCigar(list<Vertex>& pathVertex, Cigar& cigar);

        // graph pruning
        void edgePruning(int graphEdgeLevel=1, bool calEdgeWeight=false);
        void edgePruningMix(int graphEdgeLevel=1, double graphEdgeFrac=0.1, bool calEdgeWeight=false);
        void edgePruning2(int graphVertexLevel=1, int graphEdgeLevel=1, bool calEdgeWeight=false);

        // edge weight normalization
        void edgeWeightNormByOut();

        // edge weight by integrating in and out degrees
        void edgeProbWeight();

        // mark the difference on the graph
        void markDifference();

        // build the vertex->reads
        void buildVertexReadArray();

        // build the first order markovian matrix
        void buildMarkov1Matrix();

        //----------------------------------//
        // utilities for graph variants     //
        //----------------------------------//
        // collect the variants in a path
        void pathVariantCollect(list<Vertex>& pathVertices, map<int, GenericAllele>& pathVariants);

    // alignment function
    public:
        // use probabilistic aligner
        double alignReadToGraph(string& read, GenericProbabilisticAlignment& aligner,
                                int readType = READ_CROSS_REGION,
                                double band = 1.0,
                                bool processLargeHomopolymer = true,
                                Cigar& cigar = *static_cast<Cigar*>(0),
                                int& startPosition = *static_cast<int*>(0),
                                string& alnTemplate = *static_cast<string*>(0),
                                string& alnRead = *static_cast<string*>(0));

        // use trivial scoring function
        double alignReadToGraph(string& read, double matchReward=2., double mismatchPenalty=-2.,
                                double gapOpenPenalty=-3., double gapExtendPenalty=-1.,
                                Cigar& cigar = *static_cast<Cigar*>(0),
                                int readType=READ_CROSS_REGION,
                                double band = 1.0,
                                bool processLargeHomopolymer = true,
                                string& alnTemplate = *static_cast<string*>(0),
                                string& alnRead = *static_cast<string*>(0));

        // calculate CIGAR
        void alignReadToGraphCigar(list<Vertex>& alnVertex,
                                   map<Vertex,char>& alnVertexChar,
                                   map<Vertex,list<char>>& alnVertexInsertChars,
                                   Cigar& alnCigar,
                                   int readType);

    // variable member
    public:
        int m_numVertexs;   // number of vertexs
        int m_numEdges;     // number of edges

        Vertex m_begin;
        Vertex m_end;

        map<Vertex, list<EdgeVertex>> m_inEdges;      // list of in edges to a vertex
        map<Vertex, list<EdgeVertex>> m_outEdges;     // list of out edges to a vertex
        map<Vertex, string> m_labels;                 // vertex label
        map<Vertex, int> m_precedeHomopolymerLengths; // homopolymer length before a vertex
        map<Vertex, int> m_succedeHomopolymerLengths; // homopolymer length after a vertex
        map<Vertex, bool> m_isOnGenome;               // vertex is on the backbone
        set<int> m_mismatches;                        // mark the position of mismatch on the backbone
        set<int> m_inserts;                           // mark the position of insert on the backbone
        map<int,int> m_insertSize;                    // mark the maximal size of insert
        set<int> m_deletes;                           // mark the position of delete on the backbone
        int m_numMismatch;                            // the number of mismatch loci
        map<Vertex, bool> m_isMismatch;               // vertex is aligned to another vertex on the backbone
        map<Vertex, bool> m_isInsert;                 // vertex is inserted
        map<Vertex, int> m_genomePosition;            // the position on the backbone
        map<Vertex, tuple<int,int>> m_insertPosition; // the position of insert on the backbone
        map<Vertex, list<Vertex>> m_genomeSiblings;   // the list of vertices that aligned to a vertex on the backbone
        map<Vertex, bool> m_skip;                     // the vertex is removed from the graph
        map<tuple<Vertex,Vertex>, double> m_edgeProbWeight; // the probability weight of an edge
        map<tuple<Vertex,Vertex>, double> m_edgeWeight;     // the absolute weight of an edge
        map<int,map<string,map<string,double>>> m_markov1;  // first order markovian matrix for MD
        map<int,map<string,map<string,double>>> m_markov1Prob; // the normalization of the first order markovian matrix

        vector<Vertex> m_topoSortVertexs;   // topologically sorted vertexs

        // the pool of reads
        map<string, GenericGraphRead> m_readPool;           // the set of sequencing/unique reads
        vector<string> m_reads;                             // save the reads
        map<tuple<string,int>, Vertex> m_readPosVertex;     // save the read seq excluding insert
        map<tuple<string,int,int>, Vertex> m_readPosInsert; // save the read insert
        map<tuple<string,int>, Vertex> m_readIndex2Vertex;  // index to vertex
        map<tuple<string,Vertex>, int> m_readVertex2Index;  // vertex to index
        map<string,int> m_readStartPos;                     // the start position of the read on the backbone
        map<string,int> m_readEndPos;                       // the end position of the read on the backbone
        map<string,int> m_readCount;                        // the number of read occurrence
        map<string,int> m_readIndex;                        // the index of read
        map<string,list<Vertex>> m_readVertices;            // the vertices of read
        map<tuple<string,int>, Vertex> m_readMismatches;    // mismatches on a read
        map<tuple<string,int,int>, Vertex> m_readInserts;   // insertions on a read
        set<tuple<string,int>> m_readDeletes;               // deletions on a read
        map<string,map<int,string>> m_readMd;               // matches, mismatches and deletes on a read
        int m_numRead;                                      // the number of sequencing/unique reads

        map<Vertex, vector<string>> m_vertexReadArray;

        map<Vertex, Vertex> m_insertVertexUID;         // the unique id of a set of the equivalent insert vertices

        map<int,map<Vertex,vector<int>>> m_readMatIdxByPos;              // pos->(vertex->(read))
        map<int,map<Vertex,vector<int>>> m_readMisIdxByPos;              // pos->(vertex->(read))
        map<tuple<int,int>,map<Vertex,vector<int>>> m_readInsIdxByPos;   // (pos,n)->(vertex->(read))
        map<int,vector<int>> m_readDelIdxByPos;                          // pos->(read)

    public:
        friend ostream& operator<<(ostream& out, GenericDagGraph& graph)
        {
            // sort the graph vertexs
            if (graph.m_topoSortVertexs.empty())
                graph.topologicalSorting();

            // output
            out << "digraph";
            out << "{" << endl;

            // output
            out << "\t" << "rankdir=LR;" << endl;

            // output
            for (int i=0; i<graph.m_topoSortVertexs.size(); ++i)
            {
                Vertex v = graph.m_topoSortVertexs[i];

                if (graph.m_skip[v])
                    continue;

                // output
                out << "\t";
                out << "node" << v;
                out << "[";
                out << "label=" << "\"" << graph.m_labels[v] << ":" << v << ":" << graph.m_genomePosition[v] << ":" << ((graph.m_isMismatch[v]?"1":"0")) << "\"";
                if (graph.m_isOnGenome[v])
                    out << "," << "color=\"red\"";
                out << "," << "shape=\"circle\"";
                out << "]" ;
                out << ";" << endl;

            }

            // output
            for (int i=0; i<graph.m_topoSortVertexs.size(); ++i)
            {
                Vertex u = graph.m_topoSortVertexs[i];

                if (graph.m_skip[u])
                    continue;

                // output
                for (list<EdgeVertex>::iterator iter=graph.m_outEdges[u].begin(); iter!=graph.m_outEdges[u].end(); ++iter)
                {
                    if (graph.m_skip[iter->m_id])
                        continue;

                    Vertex v = iter->m_id;
                    double w = graph.m_edgeWeight[tuple<Vertex,Vertex>(u,v)];

                    out << "\t";
                    out << "node" << u << " -> " << "node" << v;
                    out << " " << "[" << "label=" << w << "]";
                    out << ";" << endl;
                }
            }

            // output
            out << "}" << endl;


            return out;
        }
};



class ConsensusGraphTool : public GenericAbstractTool
{
public:
    // ctor
    ConsensusGraphTool();

public:
    int Run(int argc, char *argv[]);
    int Help();

public:
    int commandOptionsParser(int argc, char *argv[]);

public:
    int ConsensusSequence();

public:
    bool   HasGenomeFile;
    string GenomeFile;

    bool   HasBamAlnFile;
    string BamAlnFile;

    bool   HasRegion;
    string RegionOfInterest;

    unordered_set<string> skipReadList;

    int    topK;
    int    topK2;
    int    edgePruneLevel;
    int    minReadLength;
    int    minMapQuality;
    int    filterFlag;
    double edgePruneFrac;
    int    minUniqHit;
    bool   useDiffRead;
    bool   useMarkovian;
    string scoreMethod;
    string vcfOut;
    bool   skipIndel;

    int    verbose;

    int    OutputFormat;
};

// to call strains in a local region
class LocalStrainCallTool : public GenericAbstractTool{
public:
    // constructor
    LocalStrainCallTool();
    // destructor
    ~LocalStrainCallTool(){}

public:
    int Run(int argc, char *argv[]);
    int Help();

public:
    int commandOptionParser(int argc, char *argv[]);

    int LocalStrainInference();

public:
    bool   HasGenomeFile;
    string GenomeFile;

    bool   HasBamAlnFile;
    string BamAlnFile;

    bool   HasRegion;
    string RegionOfInterest;

    // graph searching options
    int topk;
    int topK;
    // graph construction options
    int edgePruneLevel;
    int minUniqHit;
    int minReadLength;
    int minMapQuality;
    int filterFlag;
    bool useDiffRead;
    // vcf output options
    bool vcfOut;
    bool skipIndel;
    // strain clustering options
    double dissim;
    double lambda;
    // strain inference options
    bool useGenome;
    // output
    string prefix;
    int verbose;

};

// to bin long reads or assembled contigs/scaffolds into strains
class StrainBinning : public GenericAbstractTool
{
public:
    // constructor
    StrainBinning();
    // deconstructor
    ~StrainBinning(){}
public:
    int Run(int argc, char *argv[]);
    int Help();

public:
    int commandOptionParser(int argc, char *argv[]);

    int Binning();

public:
    // genome file
    bool   HasGenomeFile;
    string GenomeFile;

    RefVector GenomeDict;

    // short sequence file
    bool   HasShortFile;
    string ShortSeqFile;

    // long sequence file
    bool   HasLongFile;
    string LongSeqFile;

    // options for computational resource
    int    CR_CORES;
    string CR_BT2LONG;

    // options for regional specification
    int    RS_width;
    int    RS_space;
    int    RS_diff;
    string RS_roi;

    // options for LocalStrainCall
    int    LSC_topk;
    int    LSC_edgePruneLevel;
    int    LSC_minUniqHit;
    int    LSC_minReadLength;
    int    LSC_minMapQuality;
    int    LSC_filterFlag;
    bool   LSC_useDiffRead;
    bool   LSC_useGenome;
    double LSC_dissim;

    // options for SimpleVarCall
    int    SVC_minReadLength;
    int    SVC_minMapQuality;
    int    SVC_filterFlag;
    double SVC_minIdentity;
    int    SVC_minBaseQuality;
    int    SVC_minSnpHit;
    double SVC_minSnpRatio;

    // options for graph clustering
    int    GC_editdist;

    // options for message
    int    MS_verbose;

};
}   // namespace

#endif // GENERICGRAPHTOOLS_H
