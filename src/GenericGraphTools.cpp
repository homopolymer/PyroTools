#include "GenericGraphTools.h"
#include "GenericFastaTools.h"
using namespace GenericSequenceTools;

#include <cmath>
using namespace std;

// build the graph from a given reference genome and a set of alignments
void GenericDagGraph::buildDagGraph(string &genome, vector<string> &reads, vector<Cigar> &cigars, vector<int> &starts, vector<int> &counts)
{
    // first, construct the backbone of graph
    buildDagGraphBackbone(genome);

    // second, update the graph by iterately giving a new alignment
    for (int i=0; i<reads.size(); ++i)
    {
        if (&counts!=0)
        {
            updateDagGraphByRead(reads[i], cigars[i], starts[i], counts[i]);
        }else
        {
            updateDagGraphByRead(reads[i], cigars[i], starts[i]);
        }
    }

    // edge weight normalization
    edgeProbWeight();
}

// construct the backbone of graph from the given reference genome
void GenericDagGraph::buildDagGraphBackbone(string &genome)
{
    Vertex t_vertex = 0;    // current vertex

    // preceding homopolymer length
    vector<int> t_precedeHomopolymerLength;
    GenericFastaTools::markPrecedeHomopolymer(genome, t_precedeHomopolymerLength);

    // succeding homopolymer length
    vector<int> t_succedeHomopolymerLength;
    GenericFastaTools::markSuccedeHomopolymer(genome, t_succedeHomopolymerLength);

    // insert a dummy begin node to the graph
    Vertex t_begin = t_vertex++;
    m_inEdges[t_begin]  = list<EdgeVertex>();
    m_outEdges[t_begin] = list<EdgeVertex>();
    m_labels[t_begin]   = DAG_BEGIN_VERTEX;
    m_precedeHomopolymerLengths[t_begin] = 0;
    m_succedeHomopolymerLengths[t_begin] = 0;
    m_isOnGenome[t_begin] = false;
    m_isMismatch[t_begin] = false;
    m_begin = t_begin;
    m_skip[t_begin] = false;

    // iterately insert a base in genome
    for (int i=0; i<genome.length(); ++i)
    {
        Vertex t_curr = t_vertex++;

        // in edge of current vertex
        EdgeVertex inEdge;
        inEdge.m_id     = t_curr-1;     // parent vertex
        inEdge.m_weight = 1;            // edge weight
        inEdge.m_count  = 1;            // edge count
        m_inEdges[t_curr] = list<EdgeVertex>(1, inEdge);

        // out edge of current vertex
        m_outEdges[t_curr] = list<EdgeVertex>();

        // out edge of parent vertex
        EdgeVertex outEdge;
        outEdge.m_id     = t_curr;      // current vertex
        outEdge.m_weight = 1;           // edge weight
        outEdge.m_count  = 1;           // edge count
        m_outEdges[t_curr-1].emplace(m_outEdges[t_curr-1].end(), outEdge);

        // vertex label
        m_labels[t_curr] = genome[i];

        // preceding homopolymer length
        m_precedeHomopolymerLengths[t_curr] = t_precedeHomopolymerLength[i];
        // succeding homopolymer length
        m_succedeHomopolymerLengths[t_curr] = t_succedeHomopolymerLength[i];

        // on genome
        m_isOnGenome[t_curr] = true;

        // not mismatch
        m_isMismatch[t_curr] = false;

        // genomic position
        m_genomePosition[t_curr] = t_curr;

        // siblings
        m_genomeSiblings[t_curr] = list<Vertex>();

        // skip state
        m_skip[t_curr] = false;
    }

    // insert the dummy end vertex
    Vertex t_end = t_vertex++;
    // in edge of end vertex
    EdgeVertex inEdge;
    inEdge.m_id     = t_end-1;
    inEdge.m_weight = 1;
    inEdge.m_count  = 1;
    m_inEdges[t_end] = list<EdgeVertex>(1, inEdge);
    // out edge of parent vertex
    EdgeVertex outEdge;
    outEdge.m_id     = t_end;
    outEdge.m_weight = 1;
    outEdge.m_count  = 1;
    m_outEdges[t_end-1].emplace(m_outEdges[t_end-1].end(), outEdge);
    // vertex label
    m_labels[t_end] = DAG_END_VERTEX;
    // preceding homopolymer length
    m_precedeHomopolymerLengths[t_end] = 0;
    // succeding homopolymer length
    m_succedeHomopolymerLengths[t_end] = 0;

    // the number of vertexs at current stage
    m_numVertexs = t_vertex;
    // the number of edges at current stage
    m_numEdges   = m_numVertexs-1;

    // on genome
    m_isOnGenome[t_end] = false;

    // not mismatch
    m_isMismatch[t_end] = false;

    // skip state
    m_skip[t_end] = false;

    m_end = t_end;
}

// update the graph by giving an alignment
void GenericDagGraph::updateDagGraphByRead(string &read, Cigar &cigar, int &start, int& count, bool sorting)
{
    Vertex prev = start;      // previous vertex
    Vertex curr;              // current vertex

    Vertex g = start+1;       // genome vertex
    int    j = 0;             // pointer in read

    curr=prev;
    for (Cigar::iterator iter=cigar.begin(); iter!=cigar.end(); ++iter)
    {
        // is match or mismatch
        if (GenericBamAlignmentTools::isMatch(iter->Type) || GenericBamAlignmentTools::isMismatch(iter->Type))
        {
            for (int i=0; i<iter->Length; i++, j++, g++)
            {
                // if read base is 'N'
                if (read[j]==Amb)
                {
                    curr = g;
                    insertEdge(prev, curr, count);
                    prev = curr;
                    continue;
                }

                // if read base is the same as reference
                if (read.substr(j,1)==m_labels[g])
                {
                    curr = g;
                    insertEdge(prev, curr, count);
                    prev = curr;
                    continue;
                }

                // if read base is different from reference
                list<Vertex>::iterator siblingIter = m_genomeSiblings[g].begin();
                for (; siblingIter!=m_genomeSiblings[g].end(); ++siblingIter)
                {
                    // it is not a new vertex
                    Vertex y = *siblingIter;
                    if (read.substr(j,1)==m_labels[y].c_str() && m_isOnGenome[y]==false && m_isMismatch[y]==true)
                    {
                        curr = y;
                        insertEdge(prev, curr, count);
                        break;
                    }
                }
                // it is a new vertex
                if (siblingIter==m_genomeSiblings[g].end())
                {
                    curr = m_numVertexs++;
                    string t_label = read.substr(j,1);
                    insertVertex(curr, t_label, false, prev);
                    insertEdge(prev, curr, count);

                    // genomic position
                    m_genomePosition[curr] = g;

                    // mismatch
                    m_isMismatch[curr] = true;

                    // genome sibling
                    m_genomeSiblings[g].emplace_back(curr);
                }

                // update previous vertex
                prev = curr;
            }
        }

        // it is delete
        if (GenericBamAlignmentTools::isDelete(iter->Type))
        {
            g += iter->Length;
            curr = g;
//            insertEdge(prev, curr);
//            prev = curr;
        }

        // it is insert
        if (GenericBamAlignmentTools::isInsert(iter->Type))
        {
            for (int i=0; i<iter->Length; i++, j++)
            {
                list<EdgeVertex>::iterator outEdgeIter = m_outEdges[prev].begin();
                for (; outEdgeIter!=m_outEdges[prev].end(); ++outEdgeIter)
                {
                    int y = outEdgeIter->m_id;
                    if (read.substr(j,1)==m_labels[y] && m_isOnGenome[y]==false && m_isMismatch[y]==false)
                    {
                        // it is not a new vertex
                        curr = y;
                        insertEdge(prev, curr, count);
                        break;
                    }
                }
                // it is a new vertex
                if (outEdgeIter==m_outEdges[prev].end())
                {
                    string t_label = read.substr(j,1);
                    if (t_label!="N")
                    {
                        curr = m_numVertexs++;
                        insertVertex(curr, t_label, false, prev);
                        insertEdge(prev, curr, count);
                    }
                }
                // update previous vertex
                prev = curr;
            }
        }
    }

//    insertEdge(prev, m_end);
    insertEdge(prev, g, count);

    // topological sorting
    if (sorting)
    {
        topologicalSorting();
    }
}

// insert edge
void GenericDagGraph::insertEdge(Vertex &prev, Vertex &curr, int &count)
{
    if (prev==curr)
        return ;

    // in edge of curr
    list<EdgeVertex>::iterator inEdgeIter = m_inEdges[curr].begin();
    for (; inEdgeIter!=m_inEdges[curr].end(); ++inEdgeIter)
    {
        if (inEdgeIter->m_id==prev)
        {
            if (&count!=0)
            {
                inEdgeIter->m_weight += count;
                inEdgeIter->m_count  += count;
            }else
            {
                inEdgeIter->m_weight += 1;
                inEdgeIter->m_count  += 1;
            }
            break;
        }
    }
    if (inEdgeIter==m_inEdges[curr].end())
    {
        EdgeVertex v;
        v.m_id = prev;
        if (&count!=0)
        {
            v.m_weight = count;
            v.m_count  = count;
        }else
        {
            v.m_weight = 1;
            v.m_count  = 1;
        }
        m_inEdges[curr].emplace(m_inEdges[curr].end(), v);
    }

    // out edge of prev
    list<EdgeVertex>::iterator outEdgeIter = m_outEdges[prev].begin();
    for (; outEdgeIter!=m_outEdges[prev].end(); ++outEdgeIter)
    {
        if (outEdgeIter->m_id==curr)
        {
            if (&count!=0)
            {
                outEdgeIter->m_weight += count;
                outEdgeIter->m_count  += count;
            }else
            {
                outEdgeIter->m_weight += 1;
                outEdgeIter->m_count  += 1;
            }
            break;
        }
    }
    if (outEdgeIter==m_outEdges[prev].end())
    {
        EdgeVertex v;
        v.m_id = curr;
        if (&count!=0)
        {
            v.m_weight = count;
            v.m_count  = count;
        }else
        {
            v.m_weight = 1;
            v.m_count  = 1;
        }
        m_outEdges[prev].emplace(m_outEdges[prev].end(), v);
    }
}

void GenericDagGraph::insertVertex(Vertex &curr, string &label, bool ong, Vertex &prev)
{
    // in edge
    m_inEdges[curr] = list<EdgeVertex>();

    // out edge
    m_outEdges[curr] = list<EdgeVertex>();

    // label
    m_labels[curr] = label;

    // on genome
    m_isOnGenome[curr] = ong;

    // not mismatch
    m_isMismatch[curr] = false;

    // skip state
    m_skip[curr] = false;

    // preceding homopolymer length
    if (m_labels[curr]==m_labels[prev])
        m_precedeHomopolymerLengths[curr] = m_precedeHomopolymerLengths[prev]+1;
    else
        m_precedeHomopolymerLengths[curr] = 0;

    // succeding homopolymer length
    m_succedeHomopolymerLengths[curr] = 0;
}

// topological sort
void dfs(GenericDagGraph& graph, Vertex& x, map<Vertex, bool>& discovered, vector<Vertex>& sorted)
{
    if (graph.m_skip[x])
        return;

    discovered[x] = true;   // change the state

    // depth first search
    list<EdgeVertex>::iterator iter=graph.m_outEdges[x].begin();
    while (iter!=graph.m_outEdges[x].end()) {

        if (!discovered[iter->m_id] && !graph.m_skip[iter->m_id])
        {
            dfs(graph, iter->m_id, discovered, sorted);
        }
        iter++;
    }

    // push into array
    sorted.push_back(x);
}

// The Algorithm Design Manual, page 181
void GenericDagGraph::topologicalSorting()
{
    map<Vertex, bool> t_discovered;
    vector<Vertex> t_stack;

    for (int i=0; i<m_numVertexs; ++i)
    {
        t_discovered[i] = false;
    }

    for (int i=0; i<m_numVertexs; ++i)
    {
        if (m_skip[i])
            continue;

        if (!t_discovered[i])
            dfs(*this, i, t_discovered, t_stack);
    }

    m_topoSortVertexs.assign(t_stack.rbegin(), t_stack.rend());
}


void getFeature(GenericDagGraph& graph, string& read,
                Vertex& v, int& j, ProbalnState& currState,
                Vertex& pv, int& pj, ProbalnState& prevState,
                Matrix3Double& viterbiReadPrecedeHomopolymerLength,
                Matrix3Double& viterbiGenomePrecedeHomopolymerLength,
                vector<int>& readPrecedeHomopolymerLength,
                ProbalnTransitionFeature& trFt,
                ProbalnReadHomopolymerFeature& rhFt,
                ProbalnGenomeHomopolymerFeature& ghFt,
                ProbalnEmissionFeature& emFt)
{
    // state transition feature
    trFt = ProbalnTransitionFeature(prevState, currState);

    // emission feature
    emFt = ProbalnEmissionFeature();
    // current state is match
    if (currState==PROBALN_MATCH)
    {
        if (read.substr(j,1)==graph.m_labels[v])
        {
            emFt = ProbalnEmissionFeature(currState, Alpha2Base[graph.m_labels[v][0]], Alpha2Base[read[j]]);
        }else if (read.substr(j,1)=="N")
        {
            emFt = ProbalnEmissionFeature(currState, Alpha2Base[graph.m_labels[v][0]], Alpha2Base[graph.m_labels[v][0]]);
        }else if (graph.m_labels[v]=="N")
        {
            emFt = ProbalnEmissionFeature(currState, Alpha2Base[read[j]], Alpha2Base[read[j]]);
        }
    }
    // current state is mismatch
    if (currState==PROBALN_MISMATCH)
    {
        if (read.substr(j,1)!=graph.m_labels[v] && read.substr(j,1)!="N" && graph.m_labels[v]!="N")
        {
            emFt = ProbalnEmissionFeature(currState, Alpha2Base[graph.m_labels[v][0]], Alpha2Base[read[j]]);
        }
    }
    // current state is insert
    if (currState==PROBALN_INSERT)
    {
        if (read.substr(j,1)!="N")
        {
            emFt = ProbalnEmissionFeature(currState, SPACE, Alpha2Base[read[j]]);
        }else
        {
            emFt = ProbalnEmissionFeature(currState, SPACE, ADENINE);
        }
    }
    // current state is delete
    if (currState==PROBALN_DELETE)
    {
        if (graph.m_labels[v]!="N")
        {
            emFt = ProbalnEmissionFeature(currState, Alpha2Base[graph.m_labels[v][0]], SPACE);
        }else
        {
            emFt = ProbalnEmissionFeature(currState, ADENINE, SPACE);
        }
    }

    // read homopolymer feature
    rhFt = ProbalnReadHomopolymerFeature();
    if (currState==PROBALN_INSERT && (prevState==PROBALN_MATCH || prevState==PROBALN_INSERT))
    {
        if (read.substr(j,1)==graph.m_labels[v])
        {
            int rhl, thl, delta;
            rhl = viterbiReadPrecedeHomopolymerLength[pv][pj][prevState]+1+1;
            thl = viterbiGenomePrecedeHomopolymerLength[pv][pj][prevState]+1;
//            int rhl = readPrecedeHomopolymerLength[j]+1;
//            int thl = graph.m_precedeHomopolymerLengths[v]+1;

            delta = rhl-thl;

            if (rhl>maxHomopolymerSize)
                rhl = maxHomopolymerSize;
            if (delta>5)
                delta = 5;

            if (rhl>=midHomopolymerSize && rhl<=maxHomopolymerSize && delta>=1 && delta<=5)
                rhFt = ProbalnReadHomopolymerFeature(currState, rhl, delta);
        }
    }

    // genome homopolymer length
    ghFt = ProbalnGenomeHomopolymerFeature();
    if (currState==PROBALN_DELETE && (prevState==PROBALN_MATCH || prevState==PROBALN_DELETE))
    {
        if (graph.m_labels[v]==read.substr(j,1))
        {
            int rhl, thl, delta;
            rhl = viterbiReadPrecedeHomopolymerLength[pv][pj][prevState]+1;
            thl = viterbiGenomePrecedeHomopolymerLength[pv][pj][prevState]+1+1;
//            int rhl = readPrecedeHomopolymerLength[j]+1;
//            int thl = graph.m_precedeHomopolymerLengths[v]+1;

            delta = thl-rhl;

            if (thl>maxHomopolymerSize)
                thl = maxHomopolymerSize;
            if (delta>5)
                delta = 5;

            if (thl>=midHomopolymerSize && thl<=maxHomopolymerSize && delta>=1 && delta<=5)
                ghFt = ProbalnGenomeHomopolymerFeature(currState, thl, delta);
        }
    }
}

void updateViterbiReadPrecedeHomopolymerLength(GenericDagGraph& graph, string& read,
                                               Vertex currV, int currJ, int currState,
                                               Vertex prevV, int prevJ, int prevState,
                                               Matrix3Double& vrhl)
{
    // current state is match
    if (currState==PROBALN_MATCH)
    {
        // previous state is match
        if (prevState==PROBALN_MATCH)
        {
            if (read[currJ]==read[prevJ])
            {
                vrhl[currV][currJ][currState] = vrhl[prevV][prevJ][prevState]+1;
            }
        }
    }

    // current state is insert
    if (currState==PROBALN_INSERT)
    {
        // previous state is match
        if (prevState==PROBALN_MATCH)
        {
            if (read[currJ]==read[prevJ])
            {
                vrhl[currV][currJ][currState] = vrhl[prevV][prevJ][prevState]+1;
            }
        }

        // previous state is insert
        if (prevState==PROBALN_INSERT)
        {
            if (read[currJ]==read[prevJ])
            {
                vrhl[currV][currJ][currState] = vrhl[prevV][prevJ][prevState]+1;
            }
        }
    }
}

void updateViterbiGenomePrecedeHomopolymerLength(GenericDagGraph& graph, string& read,
                                                 Vertex currV, int currJ, int currState,
                                                 Vertex prevV, int prevJ, int prevState,
                                                 Matrix3Double& vghl)
{
    // current state is match
    if (currState==PROBALN_MATCH)
    {
        // previous state is match
        if (prevState==PROBALN_MATCH)
        {
            if (graph.m_labels[currV]==graph.m_labels[prevV])
            {
                vghl[currV][currJ][currState] = vghl[prevV][prevJ][prevState]+1;
            }
        }
    }

    // current state is delete
    if (currState==PROBALN_DELETE)
    {
        // previous state is match
        if (prevState==PROBALN_MATCH)
        {
            if (graph.m_labels[currV]==graph.m_labels[prevV])
            {
                vghl[currV][currJ][currState] = vghl[prevV][prevJ][prevState]+1;
            }
        }
        // previous state is delete
        if (prevState==PROBALN_DELETE)
        {
            if (graph.m_labels[currV]==graph.m_labels[prevV])
            {
                vghl[currV][currJ][currState] = vghl[prevV][prevJ][prevState]+1;
            }
        }
    }
}

void removeSpaceToSpace(string& t_alnRead, string& t_alnTemplate)
{
    // remove -/-
    vector<int> t_remove;
    for (int i=0; i<t_alnRead.length(); i++)
    {
        if (t_alnRead[i]=='-' && t_alnTemplate[i]=='-')
            t_remove.push_back(i);
    }
    for (vector<int>::reverse_iterator rj=t_remove.rbegin(); rj!=t_remove.rend(); ++rj)
    {
        t_alnRead.erase(t_alnRead.begin()+(*rj));
        t_alnTemplate.erase(t_alnTemplate.begin()+(*rj));
    }
}

void smallGapRefine(string& t_alnRead, string& t_alnTemplate)
{
    vector<int> rphl;
    GenericFastaTools::markPrecedeHomopolymer(t_alnRead, rphl);
    vector<int> rshl;
    GenericFastaTools::markSuccedeHomopolymer(t_alnRead, rshl);
    vector<int> gphl;
    GenericFastaTools::markPrecedeHomopolymer(t_alnTemplate, gphl);
    vector<int> gshl;
    GenericFastaTools::markSuccedeHomopolymer(t_alnTemplate, gshl);

    for (int i=1; i<t_alnRead.length(); ++i)
    {

        if (t_alnRead[i-1]=='-' && t_alnRead[i]!='-' && t_alnRead[i]!=t_alnTemplate[i] && t_alnRead[i]!='N')
        {
            if (t_alnTemplate[i-1]==t_alnRead[i])
            {
                t_alnRead[i-1] = t_alnRead[i];
                t_alnRead[i] = '-';
                continue;
            }
        }

        if (t_alnTemplate[i-1]=='-' && t_alnTemplate[i]!='-' && t_alnTemplate[i]!=t_alnRead[i] && t_alnRead[i]!='N')
        {
            if (t_alnRead[i-1]==t_alnTemplate[i])
            {
                t_alnTemplate[i-1] = t_alnTemplate[i];
                t_alnTemplate[i] = '-';
                continue;
            }
        }

        if (t_alnTemplate[i-1]=='-' && t_alnRead[i]=='-' && gshl[i]==0)
        {
            t_alnTemplate[i-1] = t_alnTemplate[i];
            t_alnTemplate[i] = '-';
            continue;
        }

        if (t_alnTemplate[i]=='-' && t_alnRead[i-1]=='-' && rshl[i]==0)
        {
            t_alnRead[i-1] = t_alnRead[i];
            t_alnRead[i] = '-';
            continue;
        }

        if (t_alnTemplate[i]=='-' && t_alnRead[i-1]!=t_alnTemplate[i-1])
        {
            if (t_alnRead[i]==t_alnTemplate[i-1])
            {
                t_alnTemplate[i] = t_alnTemplate[i-1];
                t_alnTemplate[i-1] = '-';
                continue;
            }

            if (t_alnRead[i]==t_alnRead[i-1] && gphl[i-1]==0)
            {
                t_alnTemplate[i] = t_alnTemplate[i-1];
                t_alnTemplate[i-1] = '-';
                continue;
            }
        }

        if (t_alnRead[i]=='-' && t_alnRead[i-1]!=t_alnTemplate[i-1])
        {
            if (t_alnTemplate[i]==t_alnRead[i-1])
            {
                t_alnRead[i] = t_alnRead[i-1];
                t_alnRead[i-1] = '-';
                continue;
            }

            if (t_alnTemplate[i]==t_alnTemplate[i-1] && rphl[i-1]==0)
            {
                t_alnRead[i] = t_alnRead[i-1];
                t_alnRead[i-1] = '-';
                continue;
            }
        }
    }
    GenericBamAlignmentTools::moveInDelRight(t_alnRead, t_alnTemplate);

    // remove -/-
    removeSpaceToSpace(t_alnRead, t_alnTemplate);

    //
    for (int i=0; i<t_alnRead.length(); ++i)
    {
        if (t_alnRead[i]=='-')
        {
            for (int j=i+1; j<t_alnRead.length(); ++j)
            {
                if (t_alnRead[j]!='-' && t_alnRead[j]!=t_alnTemplate[i])
                    break;

                if (t_alnRead[j]==t_alnTemplate[i])
                {
                    t_alnRead[i] = t_alnRead[j];
                    t_alnRead[j] = '-';
                    break;
                }
            }
        }

        if (t_alnTemplate[i]=='-')
        {
            for (int j=i+1; j<t_alnRead.length(); ++j)
            {
                if (t_alnTemplate[j]!='-' && t_alnTemplate[j]!=t_alnRead[i])
                    break;

                if (t_alnTemplate[j]==t_alnRead[i])
                {
                    t_alnTemplate[i] = t_alnTemplate[j];
                    t_alnTemplate[j] = '-';
                    break;
                }
            }
        }
    }

    // remove -/-
    removeSpaceToSpace(t_alnRead, t_alnTemplate);
}

double GenericDagGraph::alignReadToGraph(string &read, GenericProbabilisticAlignment &aligner, int readType, double band, bool processLargeHomopolymer, Cigar &cigar, int &startPosition, string &alnTemplate, string &alnRead)
{
    long double ViterbiScore;

    // make sure graph is sorted
    if (m_topoSortVertexs.empty())
        topologicalSorting();


    // the number of graph vertices excluding end
    int m = m_numVertexs;
    // the number of bases in read
    int n = read.length();

    // precede homopolymer length in read
    VectorInteger readPrecedeHomopolymerLength;
    GenericFastaTools::markPrecedeHomopolymer(read, readPrecedeHomopolymerLength);
    // succede homopolymer length in read
    VectorInteger readSuccedeHomopolymerLength;
    GenericFastaTools::markSuccedeHomopolymer(read, readSuccedeHomopolymerLength);
    // left neightbor homopolymer length in read
    VectorInteger leftNeighborHomopolymerLength;
    GenericFastaTools::markLeftNeighborHomopolymer(read, leftNeighborHomopolymerLength);
    // right neighbor homopolymer length in read
    VectorInteger rightNeighborHomopolymerLength;
    GenericFastaTools::markRightNeighborHomopolymer(read, rightNeighborHomopolymerLength);

    // viterbi matrix
    Matrix3Double ViterbiMatrix;
    setMatrixValue(ViterbiMatrix, m, n, PROBALN_STATE_SPACE_SIZE, DOUBLE_NEGATIVE_INFINITY);
    // previous state
    Matrix3Integer ViterbiPreviousState;
    setMatrixValue(ViterbiPreviousState, m, n, PROBALN_STATE_SPACE_SIZE, PROBALN_UNDEFINE);
    // previous node in graph
    Matrix3Integer ViterbiPreviousVertex;
    setMatrixValue(ViterbiPreviousVertex, m, n, PROBALN_STATE_SPACE_SIZE, -1);
    // previous node in read
    Matrix3Integer ViterbiPreviousJ;
    setMatrixValue(ViterbiPreviousJ, m, n, PROBALN_STATE_SPACE_SIZE, -1);
    // read precede homopolymer length
    Matrix3Double ViterbiReadPrecedeHomopolymer;
    setMatrixValue(ViterbiReadPrecedeHomopolymer, m, n, PROBALN_STATE_SPACE_SIZE, 0);
    // genome precede homopolymer length
    Matrix3Double ViterbiGenomePrecedeHomopolymer;
    setMatrixValue(ViterbiGenomePrecedeHomopolymer, m, n, PROBALN_STATE_SPACE_SIZE, 0);

    // end matrix
    VectorDouble ViterbiEndScore;
    setVectorValue(ViterbiEndScore, m, DOUBLE_NEGATIVE_INFINITY);
    // previous state at the end
    VectorInteger ViterbiEndPreviousState;
    setVectorValue(ViterbiEndPreviousState, m, PROBALN_UNDEFINE);
    // previous node in graph
    VectorInteger ViterbiEndPreviousVertex;
    setVectorValue(ViterbiEndPreviousVertex, m, -1);
    // previous node in read
    VectorInteger ViterbiEndPreviousJ;
    setVectorValue(ViterbiEndPreviousJ, m, n-1);

    // left boundary for band computation
    VectorInteger LeftBoundary;
    setVectorValue(LeftBoundary, m, 0);
    // right boundary for band computation
    VectorInteger RightBoundary;
    setVectorValue(RightBoundary, m, n);

    // utilities for Viterbi computation
    long double t_prevScore;
    long double t_localScore;
    long double t_currScore;

    int         t_prevState;
    int         t_currState;

    int         t_currVertex;
    int         t_currJ;
    int         t_prevVertex;
    int         t_prevJ;

    long double t_maxScore;
    int         t_maxPrevState;
    int         t_maxPrevVertex;
    int         t_maxPrevJ;

    ProbalnTransitionFeature t_transFt;
    ProbalnReadHomopolymerFeature t_rhFt;
    ProbalnGenomeHomopolymerFeature t_ghFt;
    ProbalnEmissionFeature t_emisFt;

    // Viterbi computation
    for (vector<Vertex>::iterator viter=m_topoSortVertexs.begin(); viter!=m_topoSortVertexs.end(); ++viter)
    {
        int i = *viter;
        if (i==m_begin || i==m_end)
            continue;

        if (m_skip[i])
            continue;

        // for band computation
        int leftBound  = LeftBoundary[i];
        int rightBound = RightBoundary[i];

        for (int j=leftBound; j<rightBound; ++j)
        {
            if ((j==0 && readType==READ_BEGIN_REGION) || (i==m_begin+1 && j==0 && readType!=READ_BEGIN_REGION))   // previous must be begin state
            {
                t_prevState = PROBALN_BEGIN;
                t_prevVertex = m_inEdges[i].begin()->m_id;
                t_prevJ = j-1;
                for (t_currState=PROBALN_MATCH; t_currState<=PROBALN_MISMATCH; ++t_currState)
                {
                    // check match
                    if (t_currState==PROBALN_MATCH)
                    {
                        if (m_labels[i]!=read.substr(j,1) && m_labels[i]!="N" && read.substr(j,1)!="N")
                            continue;
                    }

                    // check mismatch
                    if (t_currState==PROBALN_MISMATCH)
                    {
                        if (m_labels[i]==read.substr(j,1) || m_labels[i]=="N" || read.substr(j,1)=="N")
                            continue;
                    }

                    // get features
                    getFeature(*this, read, i, j, t_currState,
                               t_prevVertex, t_prevJ, t_prevState,
                               ViterbiReadPrecedeHomopolymer, ViterbiGenomePrecedeHomopolymer,
                               readPrecedeHomopolymerLength,
                               t_transFt, t_rhFt, t_ghFt, t_emisFt);

                    // previous score alway is 0
                    t_prevScore = 0;
                    // local score
                    t_localScore = 0;
                    if (t_transFt!=ProbalnTransitionFeature())
                        t_localScore += aligner.TransitionFeatureWeight[t_transFt];
                    if (t_rhFt!=ProbalnReadHomopolymerFeature())
                        t_localScore += aligner.ReadHomopolymerFeatureWeight[t_rhFt];
                    if (t_ghFt!=ProbalnGenomeHomopolymerFeature())
                        t_localScore += aligner.GenomeHomopolymerFeatureWeight[t_ghFt];
                    if (t_emisFt!=ProbalnEmissionFeature())
                        t_localScore += aligner.EmissionFeatureWeight[t_emisFt];
                    // current score
                    t_currScore = t_prevScore+t_localScore;

                    // save computation results
                    ViterbiMatrix[i][j][t_currState]         = t_currScore;
                    ViterbiPreviousState[i][j][t_currState]  = t_prevState;
                    ViterbiPreviousVertex[i][j][t_currState] = m_inEdges[i].begin()->m_id;
                    ViterbiPreviousJ[i][j][t_currState]      = j-1;
                }

            }else   // previous can not be begin state
            {
                for (t_currState=PROBALN_MATCH; t_currState<=PROBALN_DELETE; ++t_currState)
                {
                    t_maxScore      = DOUBLE_NEGATIVE_INFINITY;
                    t_maxPrevState  = PROBALN_UNDEFINE;
                    t_maxPrevVertex = -1;
                    t_maxPrevJ      = -1;

                    for (t_prevState=PROBALN_MATCH; t_prevState<=PROBALN_DELETE; ++t_prevState)
                    {

                        // check match
                        if (t_currState==PROBALN_MATCH)
                        {
                            if (m_labels[i]!=read.substr(j,1) && m_labels[i]!="N" && read.substr(j,1)!="N")
                                continue;
                        }

                        // check mismatch
                        if (t_currState==PROBALN_MISMATCH)
                        {
                            if (m_labels[i]==read.substr(j,1) || m_labels[i]=="N" || read.substr(j,1)=="N")
                                continue;
                        }

                        if (t_prevState==PROBALN_DELETE && t_currState==PROBALN_INSERT)
                            continue;
                        if (t_prevState==PROBALN_INSERT && t_currState==PROBALN_DELETE)
                            continue;

                        // match / mismatch / deleta
                        if (t_currState==PROBALN_MATCH || t_currState==PROBALN_MISMATCH || t_currState==PROBALN_DELETE)
                        {
                            if (t_currState==PROBALN_MATCH || t_currState==PROBALN_MISMATCH)
                                t_prevJ = j-1;
                            if (t_currState==PROBALN_DELETE)
                                t_prevJ = j;

                            if (t_prevJ<0)
                                continue;

                            for (list<EdgeVertex>::iterator iter=m_inEdges[i].begin(); iter!=m_inEdges[i].end(); ++iter)
                            {
                                if (m_skip[iter->m_id])
                                    continue;

                                // previous result
                                t_prevVertex = iter->m_id;
                                t_prevScore  = ViterbiMatrix[t_prevVertex][t_prevJ][t_prevState];

                                // get features
                                getFeature(*this, read,
                                           i, j, t_currState,
                                           t_prevVertex, t_prevJ, t_prevState,
                                           ViterbiReadPrecedeHomopolymer, ViterbiGenomePrecedeHomopolymer,
                                           readPrecedeHomopolymerLength,
                                           t_transFt, t_rhFt, t_ghFt, t_emisFt);

                                // local score
                                t_localScore = 0;
                                if (t_transFt!=ProbalnTransitionFeature())
                                    t_localScore += aligner.TransitionFeatureWeight[t_transFt];
                                if (t_rhFt!=ProbalnReadHomopolymerFeature())
                                    t_localScore += aligner.ReadHomopolymerFeatureWeight[t_rhFt];
                                if (t_ghFt!=ProbalnGenomeHomopolymerFeature())
                                {
                                    t_localScore += aligner.GenomeHomopolymerFeatureWeight[t_ghFt];
                                    int t_vrphl = ViterbiReadPrecedeHomopolymer[t_prevVertex][t_prevJ][t_prevState]+2;
                                    int t_vgphl = ViterbiGenomePrecedeHomopolymer[t_prevVertex][t_prevJ][t_prevState]+2;
                                    int t_dphl  = abs(t_vrphl-t_vgphl);
                                    if (t_vgphl>maxHomopolymerSize && processLargeHomopolymer)
                                    {
                                        t_localScore += t_vrphl/maxHomopolymerSize;
                                        t_localScore += t_vgphl/maxHomopolymerSize;
                                        t_localScore += t_dphl/maxHomopolymerSize;
                                    }
                                }
                                if (t_emisFt!=ProbalnEmissionFeature())
                                    t_localScore += aligner.EmissionFeatureWeight[t_emisFt];

                                // current result
                                t_currScore = t_prevScore+t_localScore;

                                // maximal
                                if (t_maxScore<t_currScore)
                                {
                                    t_maxScore      = t_currScore;
                                    t_maxPrevState  = t_prevState;
                                    t_maxPrevVertex = t_prevVertex;
                                    t_maxPrevJ      = t_prevJ;
                                }
                            }
                        }else if (t_currState==PROBALN_INSERT)
                        {
                            // previous result
                            t_prevVertex = i;
                            t_prevJ      = j-1;

                            if (t_prevJ<0)
                                continue;

                            t_prevScore  = ViterbiMatrix[t_prevVertex][t_prevJ][t_prevState];

                            // get features
                            getFeature(*this, read, i, j, t_currState,
                                       t_prevVertex, t_prevJ, t_prevState,
                                       ViterbiReadPrecedeHomopolymer, ViterbiGenomePrecedeHomopolymer,
                                       readPrecedeHomopolymerLength,
                                       t_transFt, t_rhFt, t_ghFt, t_emisFt);

                            // local score
                            t_localScore = 0;
                            if (t_transFt!=ProbalnTransitionFeature())
                                t_localScore += aligner.TransitionFeatureWeight[t_transFt];
                            if (t_rhFt!=ProbalnReadHomopolymerFeature())
                            {
                                t_localScore += aligner.ReadHomopolymerFeatureWeight[t_rhFt];
                                int t_vrphl = ViterbiReadPrecedeHomopolymer[t_prevVertex][t_prevJ][t_prevState]+2;
                                int t_vgphl = ViterbiGenomePrecedeHomopolymer[t_prevVertex][t_prevJ][t_prevState]+2;
                                int t_dphl  = abs(t_vrphl-t_vgphl);
                                if (t_vrphl>maxHomopolymerSize && processLargeHomopolymer)
                                {
                                    t_localScore += t_vrphl/maxHomopolymerSize;
                                    t_localScore += t_vgphl/maxHomopolymerSize;
                                    t_localScore += t_dphl/maxHomopolymerSize;
                                }
                            }
                            if (t_ghFt!=ProbalnGenomeHomopolymerFeature())
                                t_localScore += aligner.GenomeHomopolymerFeatureWeight[t_ghFt];
                            if (t_emisFt!=ProbalnEmissionFeature())
                                t_localScore += aligner.EmissionFeatureWeight[t_emisFt];

                            // neighbor context
                            if (read.substr(j,1)!=m_labels[i] || (m_isOnGenome[i]!=true && m_isOnGenome[m_genomePosition[i]]!=true))
                            {
                                if (readSuccedeHomopolymerLength[j]==0 || readPrecedeHomopolymerLength[j]==0)
                                {
                                    int nl = leftNeighborHomopolymerLength[j];
                                    if (nl<rightNeighborHomopolymerLength[j])
                                        nl = rightNeighborHomopolymerLength[j];

                                    if (nl>8 && processLargeHomopolymer)
                                        t_localScore += nl/maxHomopolymerSize;
                                }
                            }

                            // current result
                            t_currScore = t_prevScore+t_localScore;

                            // maximal
                            if (t_maxScore<t_currScore)
                            {
                                t_maxScore      = t_currScore;
                                t_maxPrevState  = t_prevState;
                                t_maxPrevVertex = t_prevVertex;
                                t_maxPrevJ      = t_prevJ;
                            }
                        }
                    }

                    // save result
                    ViterbiMatrix[i][j][t_currState]         = t_maxScore;
                    ViterbiPreviousState[i][j][t_currState]  = t_maxPrevState;
                    ViterbiPreviousVertex[i][j][t_currState] = t_maxPrevVertex;
                    ViterbiPreviousJ[i][j][t_currState]      = t_maxPrevJ;
                    updateViterbiReadPrecedeHomopolymerLength(*this, read,
                                                              i, j, t_currState,
                                                              t_maxPrevVertex, t_maxPrevJ, t_maxPrevState,
                                                              ViterbiReadPrecedeHomopolymer);
                    updateViterbiGenomePrecedeHomopolymerLength(*this, read,
                                                                i, j, t_currState,
                                                                t_maxPrevVertex, t_maxPrevJ, t_maxPrevState,
                                                                ViterbiGenomePrecedeHomopolymer);
//                    // debug
//                    cout << i << ", " << j << ", ";
//                    cout << ProbalnStateName[t_currState] << ", ";
//                    cout << t_maxPrevVertex << ", " << t_maxPrevJ << ", ";
//                    cout << ProbalnStateName[t_maxPrevState] << ", ";
//                    cout << t_maxScore << endl;
                }
            }
        }
    }

    // run into end vertex
    for (vector<Vertex>::iterator viter=m_topoSortVertexs.begin(); viter!=m_topoSortVertexs.end(); ++viter)
    {
        int i = *viter;
        if (i==m_begin || i==m_end)
            continue;

        if (m_skip[i])
            continue;

        // current state
        t_currState = PROBALN_END;

        // maximal utilities
        t_maxScore      = DOUBLE_NEGATIVE_INFINITY;
        t_maxPrevState  = PROBALN_UNDEFINE;
        t_maxPrevVertex = i;
        t_maxPrevJ      = n-1;

        // previous state
        for (t_prevState=PROBALN_MATCH; t_prevState<=PROBALN_DELETE; ++t_prevState)
        {
            t_prevVertex = i;
            t_prevJ      = n-1;
            // previous result
            t_prevScore = ViterbiMatrix[t_prevVertex][t_prevJ][t_prevState];

            // local score
            t_transFt = ProbalnTransitionFeature(t_prevState, t_currState);
            t_localScore = aligner.TransitionFeatureWeight[t_transFt];

            // current score
            t_currScore = t_prevScore+t_localScore;

            // maximal
            if (t_maxScore<t_currScore)
            {
                t_maxScore     = t_currScore;
                t_maxPrevState = t_prevState;
            }
        }

        // save result
        ViterbiEndScore[i]          = t_maxScore;
        ViterbiEndPreviousState[i]  = t_maxPrevState;
        ViterbiEndPreviousVertex[i] = t_maxPrevVertex;
        ViterbiEndPreviousJ[i]      = t_maxPrevJ;
    }

    // viterbi score
    if (readType==READ_CROSS_REGION || readType==READ_BEGIN_REGION)
    {
        t_maxScore = DOUBLE_NEGATIVE_INFINITY;
        for (list<EdgeVertex>::iterator iter=m_inEdges[m_end].begin(); iter!=m_inEdges[m_end].end(); ++iter)
        {
            if (m_skip[iter->m_id])
                continue;

            int y = iter->m_id;
            if (t_maxScore < ViterbiEndScore[y])
            {
                t_maxScore   = ViterbiEndScore[y];
                t_currState  = ViterbiEndPreviousState[y];
                t_currVertex = ViterbiEndPreviousVertex[y];
                t_currJ      = ViterbiEndPreviousJ[y];
            }
        }
        ViterbiScore = t_maxScore;
    }else
    {
        t_maxScore = DOUBLE_NEGATIVE_INFINITY;
        for (int i=1; i<m; ++i)
        {
            if (i==m_end)
                continue;
            if (t_maxScore<ViterbiEndScore[i])
            {
                t_maxScore   = ViterbiEndScore[i];
                t_currState  = ViterbiEndPreviousState[i];
                t_currVertex = ViterbiEndPreviousVertex[i];
                t_currJ      = ViterbiEndPreviousJ[i];
            }
        }
        ViterbiScore = t_maxScore;
    }

    // trace back
    list<Vertex>             alignVertex;
    map<Vertex,char>         alignVertexChar;
    map<Vertex,list<char>>   alignVertexInsertChars;
    list<char> insertChars;

    while (t_currState!=PROBALN_BEGIN && t_currState!=PROBALN_UNDEFINE)
    {
//        // debug
//        cout << t_currVertex << ", " << t_currJ << ", "
//             << (m_isOnGenome[t_currVertex]?"on":"off") << ", "
//             << m_genomePosition[t_currVertex] << ", "
//             << ProbalnStateName[t_currState] << ", "
//             << m_labels[t_currVertex] << ", "
//             << read[t_currJ] << endl;

        if (t_currState==PROBALN_MATCH)
        {
            if (m_isOnGenome[t_currVertex])
            {
                // save to aligned vertices
                alignVertex.insert(alignVertex.begin(), t_currVertex);

                // save read base into vertex char
                alignVertexChar[t_currVertex] = read[t_currJ];

                // save insert chars into vertex insert chars
                if (!insertChars.empty())
                    alignVertexInsertChars[t_currVertex] = list<char>(insertChars);

                // clear insert chars
                insertChars.clear();

            }else if (m_isOnGenome[m_genomePosition[t_currVertex]])
            {

                // save to aligned vertices
                alignVertex.insert(alignVertex.begin(), m_genomePosition[t_currVertex]);

                // save read base into vertex char

                alignVertexChar[m_genomePosition[t_currVertex]] = read[t_currJ];

                // save insert chars into vertex insert chars
                if (!insertChars.empty())
                    alignVertexInsertChars[m_genomePosition[t_currVertex]] = list<char>(insertChars);

                // clear insert chars
                insertChars.clear();
            }else
            {
                insertChars.insert(insertChars.begin(), read[t_currJ]);
            }
        }

        if (t_currState==PROBALN_MISMATCH)
        {
            if (m_isOnGenome[t_currVertex])
            {
                // save to aligned vertices
                alignVertex.insert(alignVertex.begin(), t_currVertex);

                // save read base into vertex char
                alignVertexChar[t_currVertex] = read[t_currJ];

                // save insert chars into vertex insert chars
                if (!insertChars.empty())
                    alignVertexInsertChars[t_currVertex] = list<char>(insertChars);

                // clear insert chars
                insertChars.clear();

            }else
            {
                insertChars.insert(insertChars.begin(), read[t_currJ]);
            }
        }

        if (t_currState==PROBALN_INSERT)
        {
            // save insert char to insert chars
            insertChars.insert(insertChars.begin(), read[t_currJ]);
        }

        // update state variable
        t_prevState  = ViterbiPreviousState[t_currVertex][t_currJ][t_currState];
        t_prevVertex = ViterbiPreviousVertex[t_currVertex][t_currJ][t_currState];
        t_prevJ      = ViterbiPreviousJ[t_currVertex][t_currJ][t_currState];

        t_currState  = t_prevState;
        t_currVertex = t_prevVertex;
        t_currJ      = t_prevJ;
    }

    // start position on genome
    startPosition = (*alignVertex.begin())-1;

    // aligned read and aligned genome
    string t_alnRead     = "";
    string t_alnTemplate = "";

    Vertex g,h;
    list<Vertex>::iterator v = alignVertex.begin();

    g = *v;

    // if read begin in region
    if (g!=m_begin+1 && readType!=READ_BEGIN_REGION)
    {
        for (int x=m_begin+1; x<g; x++)
        {
            t_alnRead     += "-";
            t_alnTemplate += m_labels[x];
        }
    }

    // next
    h = g;
    for (; v!=alignVertex.end(); ++v)
    {
        g = *v;

        // delete
        if (g>h+1)
        {
            h += 1;
            for (; h<g; h++)
            {
                t_alnRead     += "-";
                t_alnTemplate += m_labels[h];
            }
        }

        // match / mismatch
        t_alnRead     += alignVertexChar[g];
        t_alnTemplate += m_labels[g];

        // insert
        map<Vertex, list<char>>::iterator iter = alignVertexInsertChars.find(g);
        if (iter!=alignVertexInsertChars.end())
        {
            list<char>::iterator iter2 = iter->second.begin();
            for (; iter2!=iter->second.end(); ++iter2)
            {
                t_alnRead     += (*iter2);
                t_alnTemplate += "-";
            }
        }

        // previous vertex
        h = g;
    }

    // end
    if (h+1<m_end && readType!=READ_END_REGION)
    {
        h += 1;
        for (; h<m_end; h++)
        {
            t_alnRead     += "-";
            t_alnTemplate += m_labels[h];
        }
    }

    // gap tuning
    GenericBamAlignmentTools::moveInDelRight(t_alnRead, t_alnTemplate);
    smallGapRefine(t_alnRead, t_alnTemplate);

    if (&alnRead!=0 && &alnTemplate!=0)
    {
        alnRead     = t_alnRead;
        alnTemplate = t_alnTemplate;
    }

    // cigar
    if (&cigar!=0)
    {
        GenericBamAlignmentTools::calculateCigar(t_alnRead, t_alnTemplate, cigar);
    }


    // return value
    return ViterbiScore;
}


void GenericDagGraph::calculateSuccedeHomopolymer()
{
    // make sure it is sorted
    if (m_topoSortVertexs.empty())
        topologicalSorting();

    // visit every vertex
    vector<Vertex>::reverse_iterator riter = next(m_topoSortVertexs.rbegin());
    for (; riter!=m_topoSortVertexs.rend(); ++riter)
    {
        Vertex x = *riter;

        if (m_isOnGenome[x])
            continue;

        if (m_skip[x])
            continue;

        int m=0;
        list<EdgeVertex>::iterator v = m_outEdges[x].begin();
        for (; v!=m_outEdges[x].end(); ++v)
        {
            Vertex y = v->m_id;

            if (m_labels[x]==m_labels[y])
            {
                if (m<m_succedeHomopolymerLengths[y]+1)
                {
                    m = m_succedeHomopolymerLengths[y]+1;
                }
            }
        }
        m_succedeHomopolymerLengths[x] = m;
    }
}

struct SortPathByWeightDescent
{
    bool operator()(const tuple<string, double, double, list<Vertex>>& a, const tuple<string, double, double, list<Vertex>>& b)
    {
        return (get<2>(a)>get<2>(b));
    }
};

void GenericDagGraph::topRankPaths(int k, vector<string> &pathLabels, vector<list<Vertex>>& pathVertexs, vector<double> &pathWeights)
{
    // make sure graph is sorted
    if (m_topoSortVertexs.empty())
        topologicalSorting();

    // from begin to end
    map<Vertex, vector<tuple<string, double, double, list<Vertex>>>> pathToVertex;

    for (vector<Vertex>::iterator iter=m_topoSortVertexs.begin(); iter!=m_topoSortVertexs.end(); ++iter)
    {
        Vertex v = *iter;

        if (m_skip[v])
            continue;

        // it is begin node
        if (v==m_begin)
        {
            pathToVertex[v] = vector<tuple<string, double, double, list<Vertex>>>(1, tuple<string, double, double, list<Vertex>>("", 0., 0., list<Vertex>(1, m_begin)));
            continue;
        }

        // it is middle or end node
        if (v!=m_begin)
        {
            vector<tuple<string, double, double, list<Vertex>>> pathBuffer;
            for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
            {
                Vertex u = iter2->m_id;

                if (m_skip[u])
                    continue;

                for (int i=0; i<pathToVertex[u].size(); ++i)
                {
                    string pl       = get<0>(pathToVertex[u][i]);
                    double pw       = get<1>(pathToVertex[u][i]);   // sum of path weight
                    double apw      = get<2>(pathToVertex[u][i]);   // geometric average of path weight
                    list<Vertex> pv = get<3>(pathToVertex[u][i]);

                    pv.emplace_back(v);

                    pw += log(m_edgeProbWeight[tuple<Vertex,Vertex>(u,v)]);
                    apw = pw/pv.size();

//                    pw = pw*(pv.size()-1)+log(m_edgeProbWeight[tuple<Vertex,Vertex>(u,v)]);
//                    pw /= pv.size();
//                    pw *= m_edgeProbWeight[tuple<Vertex,Vertex>(u,v)];
//                    apw = pw;

                    // push into vector
                    if (v!=m_end)
                        pathBuffer.push_back(tuple<string, double, double, list<Vertex>>(pl+m_labels[v], pw, apw, pv));
                    else
                        pathBuffer.push_back(tuple<string, double, double, list<Vertex>>(pl, pw, apw, pv));
                }
            }
            // sort path buffer
            sort(pathBuffer.begin(), pathBuffer.end(), SortPathByWeightDescent());

            // erase items more than k
            if (pathBuffer.size()>k)
            {
                pathBuffer.erase(pathBuffer.begin()+k, pathBuffer.end());
            }

            // save it
            pathToVertex[v] = pathBuffer;
        }
    }

    // save results
    for (vector<tuple<string,double,double,list<Vertex>>>::iterator iter=pathToVertex[m_end].begin();
         iter!=pathToVertex[m_end].end(); ++iter)
    {
        pathLabels.push_back(get<0>(*iter));
        if (&pathWeights!=0)
            pathWeights.push_back(get<2>(*iter));
        if (&pathVertexs!=0)
            pathVertexs.push_back(get<3>(*iter));
    }
}


bool GenericDagGraph::isConsistentGraph(int k)
{
    return false;

    // top rank path
    vector<string>       pathLabels;
    vector<double>       pathWeights;
    vector<list<Vertex>> pathVertexs;
    topRankPaths(k, pathLabels, pathVertexs, pathWeights);

    // path counting
    map<string,int> pathCount;
    int i=0;
    for (vector<string>::iterator iter=pathLabels.begin(); iter!=pathLabels.end(); ++iter, ++i)
    {
        string pl0 = *iter;
        string pl  = GenericFastaTools::sequenceSignature(pl0);

        if (pathCount.find(pl)==pathCount.end())
        {
            pathCount[pl] = 1;
        }else
        {
            pathCount[pl] += 1;
        }
    }

    // search for path of multiple time
    bool consistent = true;
    for (map<string,int>::iterator iter=pathCount.begin(); iter!=pathCount.end(); ++iter)
    {

        if (iter->second>1)
        {
            consistent = false;
        }
    }

    return consistent;

}

void GenericDagGraph::alignReadToGraphCigar(list<Vertex> &alnVertex, map<Vertex, char> &alnVertexChar, map<Vertex, list<char> > &alnVertexInsertChars, Cigar &alnCigar, int readType)
{
    // first vertex
    list<Vertex>::iterator iter = alnVertex.begin();
    Vertex v = *iter++;

    // delete at the begin
    if (v!=m_begin+1 && readType!=READ_BEGIN_REGION)
    {
        CigarOp op;
        op.Type   = 'D';
        op.Length = v-m_begin-1;
        alnCigar.push_back(op);
    }

    // iterate over vertex
    Vertex u = v;               // previous vertex
    for (; iter!=alnVertex.end(); ++iter)
    {
        // current vertex
        v = *iter;

        // delete
        if (v>u+1)
        {
            CigarOp op;
            op.Type   = 'D';
            op.Length = v-u-1;
            alnCigar.push_back(op);
        }

        // match / mismatch
        Cigar::reverse_iterator riter = alnCigar.rbegin();
        if (riter!=alnCigar.rend() && (GenericBamAlignmentTools::isMatch(riter->Type) || GenericBamAlignmentTools::isMismatch(riter->Type)))
        {
            riter->Length += 1;
        }else
        {
            CigarOp op;
            op.Type   = 'M';
            op.Length = 1;
            alnCigar.push_back(op);
        }

        // insert
        map<Vertex,list<char>>::iterator iter2 = alnVertexInsertChars.find(v);
        if (iter2!=alnVertexInsertChars.end())
        {
            CigarOp op;
            op.Type   = 'I';
            op.Length = iter2->second.size();
            alnCigar.push_back(op);
        }

        // update previous vertex
        u = v;
    }

    // delete at the end
    if (m_end>u+1 && readType!=READ_END_REGION)
    {
        CigarOp op;
        op.Type   = 'D';
        op.Length = m_end-u-1;
        alnCigar.push_back(op);
    }
}


void GenericDagGraph::pathCigar(list<Vertex> &pathVertex, Cigar &cigar)
{
    Vertex u,v;

    list<Vertex>::iterator iter=pathVertex.begin();
    for (; iter!=pathVertex.end(); ++iter)
    {
        v = *iter;
        if (v==m_begin)
        {
            ;
        }else if (v==m_end)
        {
            ;
        }else
        {
            if (m_isOnGenome[v]==true || m_isMismatch[v]==true)
            {
                Vertex x=u, y=v;
                if (m_isMismatch[u]==true)
                    x = m_genomePosition[u];
                if (m_isMismatch[v]==true)
                    y = m_genomePosition[v];

                // delete
                if (y>x+1)
                {
                    CigarOp op;
                    op.Type   = 'D';
                    op.Length = y-x-1;
                    cigar.push_back(op);
                }


                Cigar::reverse_iterator riter = cigar.rbegin();
                if (riter!=cigar.rend() && (GenericBamAlignmentTools::isMatch(riter->Type) || GenericBamAlignmentTools::isMismatch(riter->Type)))
                {
                    riter->Length++;
                }else
                {
                    CigarOp op;
                    op.Type   = 'M';
                    op.Length = 1;
                    cigar.push_back(op);
                }
            }else
            {
                Cigar::reverse_iterator riter = cigar.rbegin();
                if (riter!=cigar.rend() && GenericBamAlignmentTools::isInsert(riter->Type))
                {
                    riter->Length++;
                }else
                {
                    CigarOp op;
                    op.Type   = 'I';
                    op.Length = 1;
                    cigar.push_back(op);
                }
            }
        }

        // previous
        u = v;
    }
}

void GenericDagGraph::topRankPathsExcludeGenome(int k, vector<string> &pathLabels, vector<list<Vertex> > &pathVertexs, vector<double> &pathWeights)
{
    // search top k paths
    topRankPaths(k, pathLabels, pathVertexs, pathWeights);

    // genome
    string genome = "";
    for (Vertex v=m_begin+1; v<m_end; ++v)
    {
        genome += m_labels[v];
    }

    // find genome in paths
    int j = -1;
    for (j=0; j<pathLabels.size(); ++j)
    {
        if (pathLabels[j]==genome)
            break;
    }

    // remove genome
    if (j>=0 && j<pathLabels.size())
    {
        pathLabels.erase(pathLabels.begin()+j);
        pathVertexs.erase(pathVertexs.begin()+j);
        pathWeights.erase(pathWeights.begin()+j);
    }
}


void GenericDagGraph::edgePruning(int graphEdgeLevel, bool calEdgeWeight)
{

    // make sure graph is sorted
    if (m_topoSortVertexs.empty())
        topologicalSorting();

    vector<tuple<Vertex,Vertex>> edgeToRemove;
    // visit vertexs
    for (vector<Vertex>::iterator iter=m_topoSortVertexs.begin(); iter!=m_topoSortVertexs.end(); ++iter)
    {
        Vertex v = *iter;

        if (v==m_begin || v==m_end)
            continue;

        for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
        {
            Vertex u = iter2->m_id;

            if (u+1==v)
                continue;

            // edge count below threshold
            if (iter2->m_count<=graphEdgeLevel)
                edgeToRemove.push_back(tuple<Vertex,Vertex>(u,v));
        }

    }

    // remove edge
    for (int i=0; i<edgeToRemove.size(); ++i)
    {
        Vertex u = get<0>(edgeToRemove[i]);
        Vertex v = get<1>(edgeToRemove[i]);

        list<EdgeVertex>::iterator iter;
        // erase out edge of u
        iter=m_outEdges[u].begin();
        for (; iter!=m_outEdges[u].end(); ++iter)
        {
            if (iter->m_id==v)
            {
                m_outEdges[u].erase(iter);
                break;
            }
        }

        // erase in edge of v
        iter=m_inEdges[v].begin();
        for (; iter!=m_inEdges[v].end(); ++iter)
        {
            if (iter->m_id==u)
            {
                m_inEdges[v].erase(iter);
                break;
            }
        }
    }

    // find isolated vertex
    for (int i=0; i<m_topoSortVertexs.size(); ++i)
    {
        Vertex v = m_topoSortVertexs[i];

        if (v==m_begin || v==m_end)
            continue;

        if (m_outEdges[v].empty() || m_inEdges[v].empty())
            m_skip[v] = true;
    }

    // sorting
    topologicalSorting();

    // re-calculate edge weight
    if (calEdgeWeight)
        edgeProbWeight();
}

void GenericDagGraph::edgePruning2(int graphVertexLevel, int graphEdgeLevel, bool calEdgeWeight)
{
    int threshMismatch = graphVertexLevel;
    int threshInDel    = graphEdgeLevel;

    // make sure graph is sorted
    if (m_topoSortVertexs.empty())
        topologicalSorting();

    vector<tuple<Vertex,Vertex>> edgeToRemove;
    // visit vertexs
    for (vector<Vertex>::iterator iter=m_topoSortVertexs.begin(); iter!=m_topoSortVertexs.end(); ++iter)
    {
        Vertex v = *iter;

        if (v==m_begin || v==m_end)
            continue;

        if (m_isOnGenome[v])
        {
            for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
            {
                Vertex u = iter2->m_id;

                if (!m_isOnGenome[u])
                    continue;

                // delete
                if (iter2->m_count<=threshInDel)
                    edgeToRemove.push_back(tuple<Vertex,Vertex>(u,v));
            }
            continue;
        }

        if (m_isMismatch[v]==true)  // mismatch
        {
            double s_weight = 0;
            for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
            {
                s_weight += iter2->m_count;
            }

            // threshold
            if (s_weight<=threshMismatch)
            {
                for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
                {
                    Vertex u = iter2->m_id;
                    edgeToRemove.push_back(tuple<Vertex,Vertex>(u,v));
                }
            }
        }else   // insert
        {
            double s_weight = 0;
            for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
            {
                s_weight += iter2->m_count;
            }

            // threshold
            if (s_weight<=threshInDel)
            {
                for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
                {
                    Vertex u = iter2->m_id;
                    edgeToRemove.push_back(tuple<Vertex,Vertex>(u,v));
                }
            }
        }
    }

    // remove edge
    for (int i=0; i<edgeToRemove.size(); ++i)
    {
        Vertex u = get<0>(edgeToRemove[i]);
        Vertex v = get<1>(edgeToRemove[i]);

        list<EdgeVertex>::iterator iter;
        // erase out edge of u
        iter=m_outEdges[u].begin();
        for (; iter!=m_outEdges[u].end(); ++iter)
        {
            if (iter->m_id==v)
            {
                m_outEdges[u].erase(iter);
                break;
            }
        }

        // erase in edge of v
        iter=m_inEdges[v].begin();
        for (; iter!=m_inEdges[v].end(); ++iter)
        {
            if (iter->m_id==u)
            {
                m_inEdges[v].erase(iter);
                break;
            }
        }
    }

    // find isolated vertex
    for (int i=0; i<m_topoSortVertexs.size(); ++i)
    {
        Vertex v = m_topoSortVertexs[i];

        if (v==m_begin || v==m_end)
            continue;

        if (m_outEdges[v].empty() || m_inEdges[v].empty())
            m_skip[v] = true;
    }

    // sorting
    topologicalSorting();

    // re-calculate edge weight
    if (calEdgeWeight)
        edgeProbWeight();
}


void GenericDagGraph::edgeWeightNormByOut()
{
    // normalize out edge
    for (vector<Vertex>::iterator iter=m_topoSortVertexs.begin(); iter!=m_topoSortVertexs.end(); ++iter)
    {
        Vertex v = *iter;
        double z = 0;
        for (list<EdgeVertex>::iterator iter2=m_outEdges[v].begin(); iter2!=m_outEdges[v].end(); ++iter2)
        {
            z += iter2->m_count;
        }
        for (list<EdgeVertex>::iterator iter2=m_outEdges[v].begin(); iter2!=m_outEdges[v].end(); ++iter2)
        {
            iter2->m_weight /= z;
        }
    }
}


void GenericDagGraph::edgeProbWeight()
{
    // make sure graph is topological sorted
    if (m_topoSortVertexs.empty())
        topologicalSorting();

    map<Vertex, double> t_inEdgeWeight;

    // total weight of in_edge of a vertex
    for (vector<Vertex>::iterator iter=m_topoSortVertexs.begin(); iter!=m_topoSortVertexs.end(); ++iter)
    {
        Vertex v = *iter;
        double w = 0;
        if (v==m_begin)
        {
            for (list<EdgeVertex>::iterator iter2=m_outEdges[v].begin(); iter2!=m_outEdges[v].end(); ++iter2)
            {
                w += iter2->m_count;
            }
        }else
        {
            for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
            {
                w += iter2->m_count;
            }
        }
        t_inEdgeWeight[v] = w;

    }

    // compute edge weight
    for (vector<Vertex>::iterator iter=m_topoSortVertexs.begin(); iter!=m_topoSortVertexs.end(); ++iter)
    {
        Vertex u = *iter;
        for (list<EdgeVertex>::iterator iter2=m_outEdges[u].begin(); iter2!=m_outEdges[u].end(); ++iter2)
        {
            Vertex v = iter2->m_id;
            double w = t_inEdgeWeight[u];
            m_edgeProbWeight[tuple<Vertex,Vertex>(u,v)] = iter2->m_count/w;
        }
    }
}
