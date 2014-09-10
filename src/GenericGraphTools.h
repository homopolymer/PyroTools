#ifndef GENERICGRAPHTOOLS_H
#define GENERICGRAPHTOOLS_H

#include <map>
#include <list>
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;

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

class GenericDagGraph
{
    // ctor & dtor
    public:
        GenericDagGraph()
        {
            m_numVertexs = 0;
            m_numEdges = 0;
        }

        ~GenericDagGraph()
        {
            ;
        }

    // function member
    public:
        // construct the graph by giving a reference and a set of alignments
        void buildDagGraph(string& genome, vector<string>& reads, vector<Cigar>& cigars, vector<int>& starts, vector<int>& counts=*static_cast<vector<int>*>(0));

        // build the graph backbone from a given reference sequence
        void buildDagGraphBackbone(string& genome);

        // update the graph by giving an alignment
        void updateDagGraphByRead(string& read, Cigar& cigar, int& start, int& count=*static_cast<int*>(0), bool sorting=false);

        // insert edge
        void insertEdge(Vertex& prev, Vertex& curr, int& count=*static_cast<int*>(0));

        // insert vertex
        void insertVertex(Vertex& curr, string& label, bool ong, Vertex& prev);

        // topological sorting
        void topologicalSorting();

        // compute succede homopolymer length
        void calculateSuccedeHomopolymer();

        // search top k path
        void topRankPaths(int k, vector<string>& pathLabels,
                          vector<list<Vertex>>& pathVertexs=*static_cast<vector<list<Vertex>>*>(0),
                          vector<double>& pathWeights=*static_cast<vector<double>*>(0));

        void topRankPathsExcludeGenome(int k, vector<string>& pathLabels,
                                       vector<list<Vertex>>& pathVertexs=*static_cast<vector<list<Vertex>>*>(0),
                                       vector<double>& pathWeights=*static_cast<vector<double>*>(0));

        // consistency evaluation
        bool isConsistentGraph(int k=30);

        // convert path to cigar
        void pathCigar(list<Vertex>& pathVertex, Cigar& cigar);

        // graph pruning
        void edgePruning(int graphEdgeLevel=1, bool calEdgeWeight=false);
        void edgePruning2(int graphVertexLevel=1, int graphEdgeLevel=1, bool calEdgeWeight=false);

        // edge weight normalization
        void edgeWeightNormByOut();

        // edge weight by integrating in and out degrees
        void edgeProbWeight();

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
        map<Vertex, bool> m_isOnGenome;
        map<Vertex, bool> m_isMismatch;
        map<Vertex, int> m_genomePosition;
        map<Vertex, list<Vertex>> m_genomeSiblings;
        map<Vertex, bool> m_skip;
        map<tuple<Vertex,Vertex>, double> m_edgeProbWeight;

        vector<Vertex> m_topoSortVertexs;   // topologically sorted vertexs


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
                out << "label=" << "\"" << graph.m_labels[v] << ":" << graph.m_genomePosition[v] << ":" << ((graph.m_isMismatch[v]?"1":"0")) << "\"";
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
//                    double w = graph.m_edgeProbWeight[tuple<Vertex,Vertex>(u,v)];

                    out << "\t";
                    out << "node" << u << " -> " << "node" << v;
//                    out << " " << "[" << "label=" << setprecision(3) << w << "]";
                    out << ";" << endl;
                }
            }

            // output
            out << "}" << endl;


            return out;
        }
};



class GenericGraphTools
{
public:
    GenericGraphTools()
    {
        ;
    }
};

}   // namespace

#endif // GENERICGRAPHTOOLS_H
