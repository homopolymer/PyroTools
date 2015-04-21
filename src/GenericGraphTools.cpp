#include "GenericGraphTools.h"
#include "GenericFastaTools.h"
using namespace GenericSequenceTools;

#include <list>
#include <vector>
#include <tuple>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <parallel/algorithm>
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

    // third, merge the identical vertices by traversing in backward direction
//    backwardRefineDagGraph();

    // mark the difference
    markDifference();

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
    m_isOnGenome[t_begin]     = false;
    m_isMismatch[t_begin]     = false;
    m_isInsert[t_begin]       = false;
    m_begin                   = t_begin;
    m_skip[t_begin]           = false;
    m_genomePosition[t_begin] = t_begin;

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

        // not insert
        m_isInsert[t_curr] = false;

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
    m_numVertexs            = t_vertex;
    // the number of edges at current stage
    m_numEdges              = m_numVertexs-1;
    // on genome
    m_isOnGenome[t_end]     = false;
    // not mismatch
    m_isMismatch[t_end]     = false;
    // not insert
    m_isInsert[t_end]       = false;
    // skip state
    m_skip[t_end]           = false;
    // genome position
    m_genomePosition[t_end] = t_end;

    m_end = t_end;
}

// update the graph by giving an alignment
void GenericDagGraph::updateDagGraphByRead(string &read, Cigar &cigar, int &start, int& count, bool sorting)
{
    int index = 0;            // the index of vertex
    Vertex prev = start;      // previous vertex
    Vertex curr;              // current vertex

    Vertex g = start+1;       // genome vertex
    int    j = 0;             // pointer in read

    list<Vertex> verticesOfRead;

    curr=prev;
    for (Cigar::iterator iter=cigar.begin(); iter!=cigar.end(); ++iter)
    {
        // is match or mismatch
        if (GenericBamAlignmentTools::isMatch(iter->Type) || GenericBamAlignmentTools::isMismatch(iter->Type))
        {
            for (int i=0; i<iter->Length; i++, j++, g++, index++)
            {
                // if read base is 'N'
                if (read[j]==Amb)
                {
                    curr = g;
                    insertEdge(prev, curr, count);
                    prev = curr;

                    // save vertex
                    verticesOfRead.emplace_back(curr);
                    // save vertex
                    m_readPosVertex[tuple<string,int>(read,g)] = curr;
                    // save vertex
                    m_readVertex2Index[tuple<string,Vertex>(read,curr)] = index;
                    // save vertex
                    m_readIndex2Vertex[tuple<string,int>(read,index)] = curr;

                    // index by pos and vertex
                    m_readMatIdxByPos[g][curr].emplace_back(m_numRead);

                    // save read MD
                    m_readMd[read][g] = read[j];

                    continue;
                }

                // if read base is the same as reference
                if (read.substr(j,1)==m_labels[g])
                {
                    curr = g;
                    insertEdge(prev, curr, count);
                    prev = curr;

                    // save vertex
                    verticesOfRead.emplace_back(curr);
                    // save vertex
                    m_readPosVertex[tuple<string,int>(read,g)] = curr;
                    // save vertex
                    m_readVertex2Index[tuple<string,Vertex>(read,curr)] = index;
                    // save vertex
                    m_readIndex2Vertex[tuple<string,int>(read,index)] = curr;

                    // index by pos and vertex
                    m_readMatIdxByPos[g][curr].emplace_back(m_numRead);

                    // save read MD
                    m_readMd[read][g] = read[j];

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

                        // save vertex
                        verticesOfRead.emplace_back(curr);
                        // save vertex
                        m_readPosVertex[tuple<string,int>(read,g)] = curr;
                        // save vertex
                        m_readVertex2Index[tuple<string,Vertex>(read,curr)] = index;
                        // save vertex
                        m_readIndex2Vertex[tuple<string,int>(read,index)] = curr;

                        // save read mismatch
                        m_readMismatches[tuple<string,int>(read,g)] = curr;

                        // index by pos and vertex
                        m_readMisIdxByPos[g][curr].emplace_back(m_numRead);

                        // save read MD
                        m_readMd[read][g] = read[j];

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

                    // count the mismatch
                    if (m_mismatches.find(g)==m_mismatches.end())
                        m_numMismatch++;

                    // mismatch
                    m_isMismatch[curr] = true;

                    // not insert
                    m_isInsert[curr] = false;

                    // not on backbone
                    m_isOnGenome[curr] = false;

                    // not skip
                    m_skip[curr] = false;

                    // genome sibling
                    m_genomeSiblings[g].emplace_back(curr);

                    // save vertex
                    verticesOfRead.emplace_back(curr);
                    // save vertex
                    m_readPosVertex[tuple<string,int>(read,g)] = curr;
                    // save vertex
                    m_readVertex2Index[tuple<string,Vertex>(read,curr)] = index;
                    // save vertex
                    m_readIndex2Vertex[tuple<string,int>(read,index)] = curr;

                    // save read mismatch
                    m_readMismatches[tuple<string,int>(read,g)] = curr;

                    // index by pos and vertex
                    m_readMisIdxByPos[g][curr] = vector<int>(1,m_numRead);

                    // save read md
                    m_readMd[read][g] = read[j];
                }

                // update previous vertex
                prev = curr;
            }
            continue;
        }

        // it is delete
        if (GenericBamAlignmentTools::isDelete(iter->Type))
        {
            int g0 = g;

            g += iter->Length;
            curr = g;

            // save read delete
            for (int gg=g0; gg<g; gg++)
            {
                m_readDeletes.insert(tuple<string,int>(read,gg));
            }

            // index by pos
            for (int gg=g0; gg<g; gg++)
            {
                if (m_readDelIdxByPos.find(gg) == m_readDelIdxByPos.end())
                {
                    m_readDelIdxByPos[gg] = vector<int>(1,m_numRead);
                }else
                {
                    m_readDelIdxByPos[gg].emplace_back(m_numRead);
                }
            }

            // save read MD
            for (int gg=g0; gg<g; gg++)
            {
                m_readMd[read][gg] = "-";
            }

            continue;
        }

        // it is insert
        if (GenericBamAlignmentTools::isInsert(iter->Type))
        {
            for (int i=0; i<iter->Length; i++, j++, index++)
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
                        // save vertex
                        verticesOfRead.emplace_back(curr);
                        // save vertex
                        m_readPosInsert[tuple<string,int,int>(read,g-1,i)] = curr;
                        // save vertex
                        m_readVertex2Index[tuple<string,Vertex>(read,curr)] = index;
                        // save vertex
                        m_readIndex2Vertex[tuple<string,int>(read,index)] = curr;

                        // save read insert
                        m_readInserts[tuple<string,int,int>(read,g-1,i)] = curr;

                        // index by pos and vertex
                        m_readInsIdxByPos[tuple<int,int>(g-1,i)][curr].emplace_back(m_numRead);

                        break;
                    }
                }
                // it is a new vertex
                if (outEdgeIter==m_outEdges[prev].end())
                {
                    string t_label = read.substr(j,1);
                    curr = m_numVertexs++;
                    insertVertex(curr, t_label, false, prev);
                    insertEdge(prev, curr, count);
                    // insert
                    m_isInsert[curr] = true;
                    // not mismatch
                    m_isMismatch[curr] = false;
                    // not on backbone
                    m_isOnGenome[curr] = false;
                    // the position on the backbone
                    m_genomePosition[curr] = g-1;
                    // not skip
                    m_skip[curr] = false;
                    // insert position
                    m_insertPosition[curr] = tuple<int,int>(g-1,i);
                    // save vertex
                    verticesOfRead.emplace_back(curr);
                    // save vertex
                    m_readPosInsert[tuple<string,int,int>(read,g-1,i)] = curr;
                    // save vertex
                    m_readVertex2Index[tuple<string,Vertex>(read,curr)] = index;
                    // save vertex
                    m_readIndex2Vertex[tuple<string,int>(read,index)] = curr;

                    // save read insert
                    m_readInserts[tuple<string,int,int>(read,g-1,i)] = curr;

                    // index by pos and vertex
                    m_readInsIdxByPos[tuple<int,int>(g-1,i)][curr] = vector<int>(1,m_numRead);

                }

                // update previous vertex
                prev = curr;
            }
            continue;
        }
    }

    while (g!=m_end && m_labels[g]==m_labels[prev])
    {
        g++;
    }
    insertEdge(prev, g, count);


    // save read information
    GenericGraphRead graphRead;
    graphRead.name       = read;
    graphRead.startPos   = start+1;
    graphRead.vertexList = verticesOfRead;
    graphRead.count      = count;
    m_readPool[read]     = graphRead;
    m_numRead           += 1;

    // save the read
    m_reads.emplace_back(read);
    // the start position of the read on the backbone
    m_readStartPos[read] = m_genomePosition[*verticesOfRead.begin()];
    // the end position of the read on the backbone
    m_readEndPos[read]   = m_genomePosition[*verticesOfRead.rbegin()];
    // the number of read occurrence
    m_readCount[read]    = count;
    // the index of read
    m_readIndex[read] = m_numRead-1;
    // the vertices of read
    m_readVertices[read] = verticesOfRead;

    // topological sorting
    if (sorting)
    {
        topologicalSorting();
    }

    // debug
    if (m_readIndex[read]==-1)
    {
        string sigar;
        GenericBamAlignmentTools::convertCigarToString(cigar, sigar);

        cerr << read << endl;
        cerr << sigar << endl;
        for (auto v : verticesOfRead)
        {
            cerr << v << " ~ " << m_readVertex2Index[tuple<string,Vertex>(read,v)] << endl;
        }
        cerr << m_end << endl;
        exit(0);
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

typedef long double LDBL;
struct SortPathByWeightDescent
{
    bool operator()(const tuple<string, LDBL, LDBL, list<Vertex>>& a, const tuple<string, LDBL, LDBL, list<Vertex>>& b)
    {
        return (get<2>(a)>get<2>(b));
    }
};

void GenericDagGraph::topRankPaths(int k, vector<string> &pathLabels, vector<list<Vertex>>& pathVertexs, vector<double> &pathWeights)
{
    // make sure graph is sorted
    topologicalSorting();

    // from begin to end
    map<Vertex, vector<tuple<string, LDBL, LDBL, list<Vertex>>>> pathToVertex;

    for (vector<Vertex>::iterator iter=m_topoSortVertexs.begin(); iter!=m_topoSortVertexs.end(); ++iter)
    {
        Vertex v = *iter;

        if (m_skip[v])
            continue;

        // it is begin node
        if (v==m_begin)
        {
            pathToVertex[v] = vector<tuple<string, LDBL, LDBL, list<Vertex>>>(1, tuple<string, LDBL, LDBL, list<Vertex>>("", 1., 1., list<Vertex>(1, m_begin)));
            continue;
        }

        // it is middle or end node
        if (v!=m_begin)
        {
            vector<tuple<string, LDBL, LDBL, list<Vertex>>> pathBuffer;
            for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
            {
                Vertex u = iter2->m_id;

                if (m_skip[u])
                    continue;

                for (int i=0; i<pathToVertex[u].size(); ++i)
                {
                    string pl       = get<0>(pathToVertex[u][i]);
                    LDBL pw         = get<1>(pathToVertex[u][i]);   // sum of path weight
                    LDBL apw        = get<2>(pathToVertex[u][i]);   //
                    list<Vertex> pv = get<3>(pathToVertex[u][i]);

                    pv.emplace_back(v);

                    LDBL ow = 0;
                    for (auto uo : m_outEdges[u])
                        ow += uo.m_count;

                    LDBL iw = 0;
                    for (auto io : m_inEdges[v])
                        iw += io.m_count;

                    LDBL ew = m_edgeWeight[tuple<Vertex,Vertex>(u,v)];
                    pw += 2*log(ew)-log(ow)-log(iw);
                    apw = pw/pv.size();

                    // push into vector
                    if (v!=m_end)
                        pathBuffer.push_back(tuple<string, LDBL, LDBL, list<Vertex>>(pl+m_labels[v], pw, apw, pv));
                    else
                        pathBuffer.push_back(tuple<string, LDBL, LDBL, list<Vertex>>(pl, pw, apw, pv));
                }
            }
            // sort path buffer
            if (k<=100){
                sort(pathBuffer.begin(), pathBuffer.end(), SortPathByWeightDescent());
            }else{
                __gnu_parallel::sort(pathBuffer.begin(), pathBuffer.end(), SortPathByWeightDescent());
            }


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
    for (vector<tuple<string,LDBL,LDBL,list<Vertex>>>::iterator iter=pathToVertex[m_end].begin();
         iter!=pathToVertex[m_end].end(); ++iter)
    {
        pathLabels.push_back(get<0>(*iter));
        if (&pathWeights!=0)
            pathWeights.push_back(get<2>(*iter)*k);
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
    topologicalSorting();

    vector<tuple<Vertex,Vertex>> edgeToRemove;
    // visit vertexs
    for (auto iter=m_topoSortVertexs.rbegin(); iter!=m_topoSortVertexs.rend(); ++iter)
    {
        Vertex v = *iter;

        if (v==m_begin || v==m_end)
            continue;

        if (m_skip[v])
            continue;

        for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
        {
            Vertex u = iter2->m_id;

            if (m_isOnGenome[u] && m_isOnGenome[v] && u+1==v)
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
    for (auto iter=m_topoSortVertexs.rbegin(); iter!=m_topoSortVertexs.rend(); iter++)
    {
        Vertex v = *iter;

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

void GenericDagGraph::edgePruningMix(int graphEdgeLevel, double graphEdgeFrac, bool calEdgeWeight)
{

    // make sure graph is sorted
    topologicalSorting();

    // normalized by out-edge
    edgeWeightNormByOut();

    vector<tuple<Vertex,Vertex>> edgeToRemove;
    // visit vertexs
    for (auto iter=m_topoSortVertexs.rbegin(); iter!=m_topoSortVertexs.rend(); ++iter)
    {
        Vertex v = *iter;

        if (v==m_begin)
            continue;

        if (m_skip[v])
            continue;

        for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
        {
            Vertex u = iter2->m_id;

            if (m_isOnGenome[u] && m_isOnGenome[v] && u+1==v)
                continue;

            // edge count below threshold
            if (iter2->m_count<=graphEdgeLevel ||
                    m_edgeProbWeight[tuple<Vertex,Vertex>(u,v)]<=graphEdgeFrac)
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
    for (auto iter=m_topoSortVertexs.begin(); iter!=m_topoSortVertexs.end(); iter++)
    {
        Vertex v = *iter;

        if (v==m_begin || v==m_end)
            continue;

        if (m_outEdges[v].empty() || m_inEdges[v].empty())
        {
            // in-edges
            for (auto eu : m_inEdges[v])
            {
                for (auto vptr=m_outEdges[eu.m_id].begin(); vptr!=m_outEdges[eu.m_id].end(); vptr++)
                {
                    if (vptr->m_id==v)
                    {
                        m_outEdges[eu.m_id].erase(vptr);
                    }
                }
            }
            // out-edges
            for (auto evv : m_outEdges[v])
            {
                for (auto vptr=m_inEdges[evv.m_id].begin(); vptr!=m_inEdges[evv.m_id].end(); vptr++)
                {
                    if (vptr->m_id==v)
                    {
                        m_inEdges[evv.m_id].erase(vptr);
                    }
                }
            }

            m_skip[v] = true;
        }
    }

    // sorting
    topologicalSorting();

    // re-calculate edge weight
    if (calEdgeWeight)
        edgeProbWeight();

    // mark the differences
    markDifference();
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
    for (auto iter=m_topoSortVertexs.rbegin(); iter!=m_topoSortVertexs.rend(); ++iter)
    {
        Vertex v = *iter;

        if (v==m_begin || v==m_end)
            continue;

        if (m_skip[v])
            continue;

        if (m_isOnGenome[v])
        {
            for (list<EdgeVertex>::iterator iter2=m_inEdges[v].begin(); iter2!=m_inEdges[v].end(); ++iter2)
            {
                Vertex u = iter2->m_id;

                if (m_isOnGenome[u])
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
            m_edgeWeight[tuple<Vertex,Vertex>(u,v)] = iter2->m_count;
        }
    }
}

void GenericDagGraph::markDifference()
{
    // clear the previous marking
    m_mismatches.clear();
    m_inserts.clear();
    m_insertSize.clear();
    m_deletes.clear();

    // make sure topological sort
    topologicalSorting();

    for (auto v : m_topoSortVertexs)
    {
        if (v==m_begin || v==m_end)
            continue;

        if (m_skip[v])
            continue;

        // mismatch and insert
        if (m_isMismatch[v])
        {
            m_mismatches.emplace(m_genomePosition[v]);
        }else if (m_isInsert[v])
        {
            auto ptr = m_inserts.find(m_genomePosition[v]);
            if (ptr == m_inserts.end())
            {
                m_inserts.emplace(m_genomePosition[v]);
                m_insertSize[m_genomePosition[v]] = 1;
            }else
            {
                tuple<int,int> p = m_insertPosition[v];
                if (get<1>(p)>=m_insertSize[m_genomePosition[v]])
                    m_insertSize[m_genomePosition[v]] ++;
            }
        }

        // delete
        for (auto eu : m_inEdges[v])
        {
            auto u = eu.m_id;

            if (m_skip[u])
                continue;

            if (m_genomePosition[u]+1<m_genomePosition[v])
            {
                for (int g=m_genomePosition[u]+1; g<m_genomePosition[v]; g++)
                {
                    m_deletes.emplace(g);
                }
            }
        }
    }
}

void GenericDagGraph::buildVertexReadArray()
{
    // topological sorting
    topologicalSorting();

    for (auto v : m_topoSortVertexs)
    {
        if (v==m_begin || v==m_end)
            continue;

        if (m_skip[v])
            continue;

        int gv = m_genomePosition[v];
        m_vertexReadArray[v] = vector<string>(m_numRead);
        for (int i=0; i<m_numRead; i++)
        {
            string read = m_reads[i];
            if (m_readStartPos[read]>gv)
            {
                m_vertexReadArray[v][i]="^";
            }else if (m_readEndPos[read]<gv)
            {
                m_vertexReadArray[v][i]="$";
            }else
            {
                if (m_isInsert[v])
                {
                    tuple<int,int> igv = m_insertPosition[v];
                    auto ptr = m_readPosInsert.find(tuple<string,int,int>(read,get<0>(igv),get<1>(igv)));
                    if (ptr==m_readPosInsert.end())
                    {
                        m_vertexReadArray[v][i]="-";
                    }else
                    {
                        m_vertexReadArray[v][i]=m_labels[ptr->second];
                    }
                }else
                {
                    auto ptr = m_readPosVertex.find(tuple<string,int>(read,gv));
                    if (ptr==m_readPosVertex.end())
                    {
                        m_vertexReadArray[v][i]="-";
                    }else
                    {
                        m_vertexReadArray[v][i]=m_labels[ptr->second];
                    }
                }
            }
        }
    }
}

void GenericDagGraph::buildMarkov1Matrix(){
    map<int,map<string,string>> data;
    for (int p=m_begin+1; p<m_end; p++){
        if (m_mismatches.count(p)==0)
            continue;
        for (auto r : m_reads){
            if (m_readStartPos[r]>p) data[p][r] = "^";
            else if (m_readEndPos[r]<p) data[p][r] = "$";
            else data[p][r] = m_readMd[r][p];
        }
    }

    int pu = -1;
    for (int pv=m_begin+1; pv<m_end; pv++){
        if (m_mismatches.count(pv)==0)
            continue;
        if (pu==-1){
            pu = pv;
            continue;
        }
        for (auto r : m_reads){
            string a = data[pu][r];
            string b = data[pv][r];
            if (a=="^" || a=="$" || b=="^" || b=="$")
                continue;
            m_markov1[pu][a][b] = 0;
        }
        for (auto r : m_reads){
            string a = data[pu][r];
            string b = data[pv][r];
            if (a=="^" || a=="$" || b=="^" || b=="$")
                continue;
            m_markov1[pu][a][b] += m_readCount[r];
        }
        pu = pv;
    }

    // normalization
    for (auto p : m_markov1){
        for (auto a : p.second){
            double z = 0;
            for (auto b : a.second){
                z += b.second;
            }
            for (auto b : a.second){
                m_markov1Prob[p.first][a.first][b.first] = b.second/z;
            }
        }
    }
}

void GenericDagGraph::backwardRefineDagGraph()
{
    // topological sorting
    topologicalSorting();

    // convert in-edges to out-edges
    map<Vertex, list<EdgeVertex>> outEdges = m_inEdges;
    // convert out-edges to in-edges
    map<Vertex, list<EdgeVertex>> inEdges  = m_outEdges;

    // traverse the graph in backward direction
    for (auto vptr = m_topoSortVertexs.rbegin(); vptr!=m_topoSortVertexs.rend(); vptr++)
    {
        // current vertex
        Vertex v = *vptr;

        if (m_skip[v])
            continue;

        // check merge or not
        bool merge = false;
        // collect the insert vertices
        map<string,vector<Vertex>> idenVertices;
        for (auto uu : outEdges[v])
        {
            // next vertex
            Vertex u = uu.m_id;

            if (m_isInsert[u])
            {
                auto pptr = idenVertices.find(m_labels[u]);
                if (pptr==idenVertices.end())
                {
                    idenVertices[m_labels[u]] = vector<Vertex>(u);
                }else
                {
                    idenVertices[m_labels[u]].emplace_back(u);
                }
            }
        }
        // collect the vertices to be merged
        vector<vector<Vertex>> mergeVertices;
        for (auto x: idenVertices)
        {
            if (x.second.size()>=2){
                sort(x.second.begin(),x.second.end());
                mergeVertices.emplace_back(x.second);
            }
        }

        // aggregate the links of the vertices to the first vertex
        for (auto x: mergeVertices)
        {
            Vertex x0 = x[0];
            for (int i=1; i<x.size(); i++)
            {
                Vertex xi = x[i];
                // move in-edges
                for (auto e: m_inEdges[xi])
                {
                    insertEdge(e.m_id, x0, e.m_count);
                }
                // move out-edges
                for (auto e: m_outEdges[xi])
                {
                    insertEdge(x0, e.m_id, e.m_count);
                }

                m_skip[xi] = true;
            }
        }
    }
}

/*
 *  Search the top-k path by the metrics of the alignment identity
 *
 *  The general idea is that compute the identities of all reads aligned to a path,
 *  and then sort the identities in the descending order, then pick up the 25 quantile
 *  as the score of the path.
 */

typedef struct _read2path
{
    public:
        _read2path()
            :m_align("")
            ,m_start(0)
            ,m_numMatch(0)
            ,m_numMismatch(0)
            ,m_numInsert(0)
            ,m_numDelete(0)
        {}
        _read2path(string &name)
            :m_name(name)
            ,m_align("")
            ,m_numMatch(0)
            ,m_numMismatch(0)
            ,m_numInsert(0)
            ,m_numDelete(0)
        {}
        bool empty(){return (m_numMatch+m_numMismatch+m_numInsert+m_numDelete)==0;}
    public:
        string       m_name;
        string       m_align;

        Vertex       m_start;
        int          m_numMatch;
        int          m_numMismatch;
        int          m_numInsert;
        int          m_numDelete;

        list<Vertex> m_vertices;
}Read2Path;

void GenericDagGraph::topRankPathsByIden(int k, double quantile, vector<string> &pathLabels, vector<list<Vertex> > &pathVertexs, vector<double> &pathWeights)
{
    // make sure the graph is topologically sorted
    topologicalSorting();

    // this record the information of the paths reaching the vertex
    typedef tuple<string, list<Read2Path>, double, list<Vertex>> PathInfo;
    map<Vertex, vector<PathInfo>> pathToVertex;

    // lambda function for sorting
    auto SortByIden = [](const PathInfo& a, const PathInfo &b)->bool{
        return get<2>(a)>get<2>(b);
    };

    // lambda function for computing alignment identity
    auto AlnIden = [](const int &mat, const int &mis, const int &ins, const int &del)->double{
        if (mat+mis+ins+del==0)
            return 0;
        return mat/(mat+mis+ins+del+0.);
    };

    // lambda function for the quantile
    auto Q = [quantile](vector<tuple<double,int>>& x)->double{
        // sort x in descending order
        sort(x.begin(), x.end(),[](const tuple<double,int> &a, const tuple<double,int> &b)->bool{return (get<0>(a)>get<0>(b));});
        // compute the cumulative sum
        vector<double> y(x.size(),0);
        vector<int>    s(x.size(),0);
        #pragma parallel for shared(y,s)
        for (int i=0; i<x.size(); i++)
        {
            if (i==0)
            {
                y[i] = get<0>(x[i])*get<1>(x[i]);
                s[i] = get<1>(x[i]);
            }
            else
            {
                y[i] = y[i-1] + get<0>(x[i])*get<1>(x[i]);
                s[i] = s[i-1] + get<1>(x[i]);
            }
        }
        // total score
        double z = *y.rbegin();
        // find the quantile
        double q;
        int i=0;
        for (; i<y.size(); i++)
        {
            if (y[i]>=z*quantile)
            {
                break;
            }
        }
        if (i==0)
        {
            q = get<0>(x[i]);
//            q = y[i]/s[i];
        }else
        {
            q = (get<0>(x[i])+get<0>(x[i-1]))/2.0;
//            q = (y[i]+y[i-1])/(s[i]+s[i-1]);
        }
        return q;
    };

    // iterate from begin to end
    for (auto v : m_topoSortVertexs)
    {
        if (m_skip[v])
            continue;

        // debug
        cerr << v << "\t" << m_end << endl;

        // if vertex is begin
        if (v==m_begin)
        {
            list<Read2Path> readAlns;
            for (auto r : m_reads)
            {
                readAlns.emplace_back(Read2Path(r));
            }

            pathToVertex[v] = vector<PathInfo>(1,PathInfo("",readAlns,0.,list<Vertex>(1,m_begin)));
            continue;
        }

        // if vertex is end
        if (v==m_end)
        {
            // collect the path
            vector<PathInfo> buffer;
            // iterate over the in-edges
            for (auto eu : m_inEdges[v])
            {
                Vertex u = eu.m_id;
                for (auto p : pathToVertex[u])
                {
                    get<3>(p).emplace_back(v);
                    buffer.emplace_back(p);
                }
            }
            // compute the indentities
            for (auto &p : buffer)
            {
                vector<tuple<double,int>> readIden;
                for (auto r : get<1>(p))
                {
                    double iden = AlnIden(r.m_numMatch, r.m_numMismatch, r.m_numInsert, r.m_numDelete);
                    readIden.emplace_back(tuple<double,int>(iden,m_readCount[r.m_name]));
                }
                double q = Q(readIden);
                get<2>(p) = q;
            }
            // sort path
            sort(buffer.begin(), buffer.end(), SortByIden);
            // save the top-k path
            vector<PathInfo> localResults;
            for (int i=0; i<buffer.size(); i++)
            {
                if (i>=k)
                    break;

                localResults.emplace_back(buffer[i]);
            }
            pathToVertex[v] = localResults;

            continue;
        }

        // if vertex is the internal node
        if (v!=m_begin && v!=m_end)
        {
            // temporary buffer
            vector<PathInfo> buffer;

            // iterate over the in-edges
            for (auto eu : m_inEdges[v])
            {
                // the id of the precede vertex
                Vertex u = eu.m_id;
                // skip if it is masked
                if (m_skip[u])
                    continue;

                // if current vertex is insert
                if (m_isInsert[v])
                {
                    // the position of the insert on the backbone
                    tuple<int,int> vp = m_insertPosition[v];
                    // iterate over all path reaching u
                    for (auto p : pathToVertex[u])
                    {
                        // iterate over all reads
                        for (auto &r : get<1>(p))
                        {
                            // the last vertex of read before reaching v
                            Vertex ru = *r.m_vertices.rbegin();

                            if (ru != *m_readVertices[r.m_name].rbegin()){
                                // if there is delete
                                auto rv_ptr = m_readPosInsert.find(tuple<string,int,int>(r.m_name, get<0>(vp), get<1>(vp)));
                                if (rv_ptr == m_readPosInsert.end())
                                {
                                    r.m_numDelete += 1;
                                    r.m_align += "D";
                                }else
                                {

                                    // current vertex on the read
                                    Vertex rv = rv_ptr->second;

                                    // there is insert
                                    int iru = m_readVertex2Index[tuple<string,Vertex>(r.m_name,ru)];
                                    int irv = m_readVertex2Index[tuple<string,Vertex>(r.m_name,rv)];
                                    if (iru+1<irv)
                                    {
                                        int i0 = m_readVertex2Index[tuple<string,Vertex>(r.m_name,ru)]+1;
                                        int i1 = m_readVertex2Index[tuple<string,Vertex>(r.m_name,rv)]-1;
                                        for (int ii=i0; ii<=i1; ii++)
                                        {
                                            Vertex rx = m_readIndex2Vertex[tuple<string,int>(r.m_name,ii)];
                                            r.m_numInsert += 1;
                                            r.m_align += "I";
                                            r.m_vertices.emplace_back(rx);
                                        }
                                    }

                                    if (m_labels[rv]==m_labels[v])
                                    {
                                        r.m_numMatch += 1;
                                        r.m_align += "M";
                                    }else
                                    {
                                        r.m_numMismatch += 1;
                                        r.m_align += "X";
                                    }

                                    // insert rv to the vertex list
                                    r.m_vertices.emplace_back(rv);
                                }
                            }
                        }
                        get<0>(p) += m_labels[v];
                        get<3>(p).emplace_back(v);
                        buffer.emplace_back(p);
                    }
                }else
                {
                    // the position of v on the backbone
                    int vp = m_genomePosition[v];
                    // iterate over all path reaching u
                    for (auto p : pathToVertex[u])
                    {
                        // iterate over all reads
                        for (auto &r : get<1>(p))
                        {
                            // read starts before v
                            if (!r.empty()){
                                // the last vertex of read before reaching v
                                Vertex ru = *r.m_vertices.rbegin();
                                if (ru != *m_readVertices[r.m_name].rbegin())
                                {
                                    // it is delete
                                    auto rv_ptr = m_readPosVertex.find(tuple<string,int>(r.m_name,vp));
                                    if (rv_ptr==m_readPosVertex.end())
                                    {
                                        int iru = m_readVertex2Index[tuple<string,Vertex>(r.m_name,ru)];
                                        int irv = iru+1;
                                        Vertex rv = m_readIndex2Vertex[tuple<string,int>(r.m_name,irv)];
                                        while (m_genomePosition[rv]<m_genomePosition[v] && m_readEndPos[r.m_name]!=m_genomePosition[rv])
                                        {
                                            rv = m_readIndex2Vertex[tuple<string,int>(r.m_name,++irv)];
                                        }
                                        irv = m_readVertex2Index[tuple<string,int>(r.m_name,rv)];
                                        // there is insert
                                        if (iru+1<irv || m_readEndPos[r.m_name]==rv)
                                        {
                                            int i0 = iru+1;
                                            int i1 = irv-1;
                                            if (m_readEndPos[r.m_name]==rv)
                                                i1 = irv;
                                            for (int ii=i0; ii<=i1; ii++)
                                            {
                                                Vertex rx = m_readIndex2Vertex[tuple<string,int>(r.m_name,ii)];
                                                r.m_numInsert += 1;
                                                r.m_align += "I";
                                                r.m_vertices.emplace_back(rx);
                                            }
                                        }

                                        r.m_numDelete += 1;
                                        r.m_align += "D";
                                    }else
                                    {
                                        // current vertex on the read
                                        Vertex rv = rv_ptr->second;

                                        // there is insert
                                        int iru = m_readVertex2Index[tuple<string,Vertex>(r.m_name,ru)];
                                        int irv = m_readVertex2Index[tuple<string,Vertex>(r.m_name,rv)];
                                        if (iru+1<irv)
                                        {
                                            int i0 = m_readVertex2Index[tuple<string,Vertex>(r.m_name,ru)]+1;
                                            int i1 = m_readVertex2Index[tuple<string,Vertex>(r.m_name,rv)]-1;
                                            for (int ii=i0; ii<=i1; ii++)
                                            {
                                                Vertex rx = m_readIndex2Vertex[tuple<string,int>(r.m_name,ii)];
                                                r.m_numInsert += 1;
                                                r.m_align += "I";
                                                r.m_vertices.emplace_back(rx);
                                            }
                                        }

                                        if (m_labels[rv]==m_labels[v])
                                        {
                                            r.m_numMatch += 1;
                                            r.m_align += "M";
                                        }else
                                        {
                                            r.m_numMismatch += 1;
                                            r.m_align += "X";
                                        }

                                        // insert rv to the vertex list
                                        r.m_vertices.emplace_back(rv);
                                    }
                                }
                            }else
                            {
                                Vertex rv = m_readPosVertex[tuple<string,int>(r.m_name,vp)];
                                // read start at v
                                if (m_readStartPos[r.m_name]==vp)
                                {
                                    if (m_labels[rv]==m_labels[v])
                                    {
                                        r.m_numMatch += 1;
                                        r.m_align += "M";
                                    }else
                                    {
                                        r.m_numMismatch += 1;
                                        r.m_align += "X";
                                    }
                                    // start vertex
                                    r.m_start = v;
                                    // insert rv to the vertex list
                                    r.m_vertices.emplace_back(rv);
                                }
                                // read start before path
                                else if (m_readStartPos[r.m_name]<vp)
                                {
                                    // before vb
                                    int i1 = m_readVertex2Index[tuple<string,Vertex>(r.m_name,rv)]-1;
                                    for (int ii=0; ii<=i1; ii++)
                                    {
                                        Vertex rx = m_readIndex2Vertex[tuple<string,int>(r.m_name,ii)];
                                        r.m_align += "S";
                                        r.m_vertices.emplace_back(rx);
                                    }
                                    // v
                                    if (m_labels[rv]==m_labels[v])
                                    {
                                        r.m_numMatch += 1;
                                        r.m_align += "M";
                                    }else
                                    {
                                        r.m_numMismatch += 1;
                                        r.m_align += "X";
                                    }
                                    // start vertex
                                    r.m_start = v;
                                    // insert rv to the vertex list
                                    r.m_vertices.emplace_back(rv);
                                }
                            }

                        }
                        get<0>(p) += m_labels[v];
                        get<3>(p).emplace_back(v);
                        buffer.emplace_back(p);
                    }
                }
            }

            // compute the quantile
            for (auto &p : buffer)
            {
                vector<tuple<double,int>> readIden;
                for (auto r : get<1>(p))
                {
                    double iden = AlnIden(r.m_numMatch, r.m_numMismatch, r.m_numInsert, r.m_numDelete);
                    readIden.emplace_back(tuple<double,int>(iden,m_readCount[r.m_name]));
                }
                double q = Q(readIden);
                get<2>(p) = q;
            }
            // sort path
            sort(buffer.begin(), buffer.end(), SortByIden);
            // save the top-k path
            vector<PathInfo> localResults;
            for (int i=0; i<buffer.size(); i++)
            {
                if (i>=k)
                    break;

                localResults.emplace_back(buffer[i]);
            }
            pathToVertex[v] = localResults;

            continue;
        }
    }

    int kk = 0;
    for (auto p : pathToVertex[m_end])
    {
        if (kk>=k)
            break;
        pathLabels.emplace_back(get<0>(p));

        if (&pathVertexs!=nullptr)
        {
            pathVertexs.emplace_back(get<3>(p));
        }

        if (&pathWeights!=nullptr)
        {
            pathWeights.emplace_back(get<2>(p));
        }

        kk ++;
    }

}

void GenericDagGraph::topRankPathsByErrMdl(int k, int k2, string scoring, bool markovian, vector<string> &pathLabels, vector<list<Vertex> > &pathVertexs, vector<double> &pathWeights)
{
    double err = 0.01;
    double le = log(err);
    double lc = log(1-err);

    map<Vertex,int> vertexIdx;
    // make sure the graph is topologically sorted
    topologicalSorting();

    // the index of vertex
    for (int i=0; i<m_topoSortVertexs.size(); i++)
        vertexIdx[m_topoSortVertexs[i]] = i;

    // build the matrix of vertex-read
    buildVertexReadArray();

    // build the first order markovian matrix
    if (markovian)
        buildMarkov1Matrix();

    // data structure for allele count
    // the number of matches and mismatches
    typedef tuple<int,int> AlleleCount;

    // data structure for read info
    // AlleleCount, previous difference and markovian feature
    typedef vector<tuple<AlleleCount,int,LDBL>> ReadInfo;

    // data structure for path
    // (label, read info, score, vertices, read scores)
    typedef tuple<string,ReadInfo,double,list<Vertex>> GraphPath;

    // data structure for vertex
    map<Vertex, vector<GraphPath>> pathToVertex;

    // sorting method for path
    auto ByDescendOrder = [](const GraphPath &a, const GraphPath &b)->bool{
        return get<2>(a)>get<2>(b);
    };

    // path scoring method
    // likelihood scoring
    auto ScorePath = [lc,le,this](GraphPath &a)->double{
        double score = 0;
        int i = 0;
        for (auto x : get<1>(a))
        {
            int n = this->m_readCount[this->m_reads[i++]];
            if (get<0>(get<0>(x))+get<1>(get<0>(x))>0)
                score += n*exp(get<0>(get<0>(x))*lc + get<1>(get<0>(x))*le + get<2>(x));
        }
        return score;
    };

    // maximal assignment scoring
    auto ScorePath2 = [lc, le, this, &ScorePath](vector<GraphPath> &a){
        vector<vector<double>> readPathScore(this->m_numRead);
        for (int i=0; i<this->m_numRead; i++)
        {
            for (auto p : a)
            {
                int n0 = get<0>(get<0>(get<1>(p)[i]));
                int n1 = get<1>(get<0>(get<1>(p)[i]));
                if (n0+n1>0)
                    readPathScore[i].emplace_back(exp(n0*lc+n1*le+get<2>(get<1>(p)[i])));
                else
                    readPathScore[i].emplace_back(0);
            }
        }
        vector<double> pathReadCount(a.size(), 0);
        for (int i=0; i<this->m_numRead; i++)
        {
            LDBL maxs = 0;
            int maxj = -1;
            int n = this->m_readCount[this->m_reads[i]];
            for (int j=0; j<a.size(); j++)
            {
                if (maxs < readPathScore[i][j])
                {
                    maxs = readPathScore[i][j];
                    maxj = j;
                }
            }
            if (maxj != -1) pathReadCount[maxj] += n;
        }
        for (int j=0; j<a.size(); j++)
        {
            get<2>(a[j]) = pathReadCount[j];
        }
    };

    // average assignment scoring
    auto ScorePath3 = [lc, le, this, &ScorePath](vector<GraphPath> &a){
        vector<vector<double>> readPathScore(this->m_numRead);
        for (int i=0; i<this->m_numRead; i++)
        {
            for (auto p : a)
            {
                int n0 = get<0>(get<0>(get<1>(p)[i]));
                int n1 = get<1>(get<0>(get<1>(p)[i]));
                if (n0+n1>0)
                    readPathScore[i].emplace_back(exp(n0*lc+n1*le+get<2>(get<1>(p)[i])));
                else
                    readPathScore[i].emplace_back(0);
            }
        }
        vector<double> pathReadCount(a.size(), 0);
        for (int i=0; i<this->m_numRead; i++)
        {
            LDBL maxt = 0;
            int n = this->m_readCount[this->m_reads[i]];
            for (int j=0; j<a.size(); j++)
            {
                maxt += readPathScore[i][j];
            }
            for (int j=0; j<a.size(); j++)
            {
                pathReadCount[j] += n*readPathScore[i][j]/(maxt+1e-9);
            }
        }
        for (int j=0; j<a.size(); j++)
        {
            get<2>(a[j]) = pathReadCount[j];
        }
    };

    // unique graph path
    auto UniquePath = [&ByDescendOrder](vector<GraphPath> &a)
    {
        vector<GraphPath> b(a.begin(), a.end());

        a.clear();
        a.shrink_to_fit();

        set<string> pathSeqs;
        for (auto p : b)
        {
            bool skip = false;
            for (auto s : pathSeqs)
            {
                if (get<0>(p).find(s)!=string::npos)
                {
                    skip = true;
                    break;
                }
                if (s.find(get<0>(p))!=string::npos)
                {
                    skip = true;
                    break;
                }
            }
            if (!skip)
            {
                a.emplace_back(p);
                pathSeqs.insert(get<0>(p));
            }
        }
    };

    // top-k path
    auto TopkPath = [&UniquePath,scoring](vector<GraphPath> &a, vector<GraphPath> &b, int k)
    {
        for (int i=0; i<b.size(); i++)
        {
            if (i>=k)
                break;
            a.emplace_back(b[i]);
        }
    };

    // iterate from begin to end
    for (auto v : m_topoSortVertexs)
    {
        // if vertex is begin
        if (v==m_begin)
        {
            pathToVertex[v] = vector<GraphPath>(1, GraphPath("", ReadInfo(m_numRead,tuple<AlleleCount,int,LDBL>(AlleleCount(0,0),-1,0)), 0, list<Vertex>(1,m_begin)));
            continue;
        }

//        // if vertex is end
//        if (v==m_end)
//        {
//            vector<GraphPath> pathBuffer;
//            // iterate over in-edges
//            for (auto eu : m_inEdges[v])
//            {
//                // previous vertex
//                auto u = eu.m_id;

//                if (m_skip[u]) continue;

//                // iterate over all paths of u
//                for (auto &p : pathToVertex[u])
//                {
//                    get<3>(p).emplace_back(m_end);
//                    pathBuffer.emplace_back(p);
//                }
//            }

//            if (scoring=="assign")
//                ScorePath2(pathBuffer);

//            // sort the path
//            sort(pathBuffer.begin(), pathBuffer.end(), ByDescendOrder);

//            pathToVertex[v] = vector<GraphPath>();
//            TopkPath(pathToVertex[v], pathBuffer, k);

//            continue;
//        }

        // if vertex is the internal node
        if (v!=m_begin)
        {
            auto vidx = vertexIdx[v];
            int pv = m_genomePosition[v];

            vector<GraphPath> vPathBuffer;
            // iterate over the in-edges
            for (auto eu : m_inEdges[v])
            {
                auto u = eu.m_id;
                auto uidx = vertexIdx[u];
                int pu = m_genomePosition[u];

                // skip if it is masked
                if (m_skip[u]) continue;

                set<int> visited_pos;
                set<tuple<int,int>> visited_ins;

                if (m_mismatches.count(pu))
                    visited_pos.insert(pu);
                if (m_isInsert[u])
                    visited_ins.insert(m_insertPosition[u]);

                vector<GraphPath> uPathBuffer(pathToVertex[u].begin(), pathToVertex[u].end());
                // there is a jump
//                if (pu+1<pv || m_inserts.count(pu)){
                    for (int xidx=uidx+1; xidx<vidx; xidx++)
                    {
                        Vertex x = m_topoSortVertexs[xidx];
                        int px = m_genomePosition[x];

                        if (px>=pv)
                            continue;

                        if (m_mismatches.count(px)){
                            if (visited_pos.count(px))
                                continue;
                            visited_pos.insert(px);
                        }
                        if (m_isInsert[x])
                        {
                            tuple<int,int> ipx = m_insertPosition[x];
                            if (visited_ins.count(ipx))
                                continue;
                            visited_ins.insert(ipx);
                        }

                        for (auto &p : uPathBuffer)
                        {
                            for (int r=0; r<m_numRead; r++)
                            {
                                if (m_vertexReadArray[x][r]=="^")
                                    continue;
                                if (m_vertexReadArray[x][r]=="$")
                                    continue;

                                if (m_vertexReadArray[x][r]!="-")
                                    get<1>(get<0>(get<1>(p)[r])) += 1;
                                else
                                    get<0>(get<0>(get<1>(p)[r])) += 1;
                            }

                            // compute the markovian feature
                            if (!m_isInsert[x] && m_mismatches.count(px) && markovian){
                                for (int r=0; r<m_numRead; r++){
                                    if (m_readStartPos[m_reads[r]]>px) continue;
                                    if (m_readEndPos[m_reads[r]]<px) continue;
                                    int t = get<1>(get<1>(p)[r]);
                                    if (t==-1){
                                        get<1>(get<1>(p)[r]) = px;
                                    }else{
                                        string a = m_readMd[m_reads[r]][t];
                                        string b = m_readMd[m_reads[r]][px];
                                        get<2>(get<1>(p)[r]) += log(m_markov1Prob[t][a][b]);
                                        get<1>(get<1>(p)[r]) = px;
                                    }
                                }
                            }
                        }                        
                    }
//                }

                vPathBuffer.insert(vPathBuffer.end(), uPathBuffer.begin(), uPathBuffer.end());
            }

            // update the alignment information
            if ((m_mismatches.count(pv) || m_deletes.count(pv) || m_isInsert[v]) && v!=m_end)
            {
                for (auto &p : vPathBuffer)
                {
                    for (int r=0; r<m_numRead; r++)
                    {
                        if (m_vertexReadArray[v][r]=="^")
                            continue;
                        if (m_vertexReadArray[v][r]=="$")
                            continue;

                        if (m_vertexReadArray[v][r]==m_labels[v])
                            get<0>(get<0>(get<1>(p)[r])) += 1;
                        else
                            get<1>(get<0>(get<1>(p)[r])) += 1;
                    }

                    if (!m_isInsert[v] && m_mismatches.count(pv) && markovian){
                        for (int r=0; r<m_numRead; r++){
                            if (m_readStartPos[m_reads[r]]>pv) continue;
                            if (m_readEndPos[m_reads[r]]<pv) continue;
                            int t = get<1>(get<1>(p)[r]);
                            if (t==-1){
                                get<1>(get<1>(p)[r]) = pv;
                            }else{
                                string a = m_readMd[m_reads[r]][t];
                                string b = m_readMd[m_reads[r]][pv];
                                get<2>(get<1>(p)[r]) += log(m_markov1Prob[t][a][b]);
                                get<1>(get<1>(p)[r]) = pv;
                            }
                        }
                    }
                }
            }

            // compute path weight and update sequence
            for (auto &p : vPathBuffer)
            {
                if (v!=m_end)
                    get<0>(p) += m_labels[v];

                if ((m_inEdges[v].size()>1 || m_mismatches.count(pv) || m_deletes.count(pv) || m_isInsert[v] || v==m_end)
                        && scoring=="likeli")
                {
                    double s = ScorePath(p);
                    get<2>(p) = s;
                }

                get<3>(p).emplace_back(v);
            }

            if ((m_inEdges[v].size()>1 || m_mismatches.count(pv) || m_deletes.count(pv) || m_isInsert[v] || v==m_end)
                    && scoring=="assign")
                ScorePath3(vPathBuffer);

            // sort the path
            sort(vPathBuffer.begin(), vPathBuffer.end(), ByDescendOrder);
            // top-k path
            pathToVertex[v] = vector<GraphPath>();
            if (v!=m_end)
                TopkPath(pathToVertex[v], vPathBuffer, k2);
            else
                TopkPath(pathToVertex[v], vPathBuffer, k);

            continue;
        }
    }

    for (auto p : pathToVertex[m_end])
    {
        pathLabels.emplace_back(get<0>(p));

        if (&pathVertexs!=nullptr)
        {
            pathVertexs.emplace_back(get<3>(p));
        }

        if (&pathWeights!=nullptr)
        {
            pathWeights.emplace_back(get<2>(p));
        }
    }
}

// the variants in a path
// note: it will use VRM to collect read information in the future
void GenericDagGraph::pathVariantCollect(list<Vertex> &pathVertices, map<int, GenericAllele> &pathVariants)
{
    Vertex u;
    for (auto v : pathVertices){
        if (v==m_begin) u = v;
        if (v==m_end) continue;

        // it is mismatch or insert
        if (!m_isOnGenome[v]){
            // it is mismatch
            if (m_isMismatch[v]){
                auto pos = m_genomePosition[v] - 1;
                string ref = m_labels[pos+1];
                string alt = m_labels[v];
                // allele
                GenericAllele allele;
                allele.m_chrID       = 0;
                allele.m_chrPosition = pos;
                allele.m_allele      = alt;
                allele.m_reference   = ref;
                // save
                pathVariants[pos]    = allele;
            }
            // it is insert
            else if (m_isInsert[v]){
                auto pos = m_genomePosition[v] - 1;
                if (pathVariants.count(pos)){
                    // update
                    pathVariants[pos].m_allele += m_labels[v];
                }else{
                    string ref = m_labels[pos+1];
                    string alt = m_labels[v];
                    // allele
                    GenericAllele allele;
                    allele.m_chrID       = 0;
                    allele.m_chrPosition = pos;
                    allele.m_allele      = alt;
                    allele.m_reference   = ref;
                    // save
                    pathVariants[pos]    = allele;
                }
            }
        }else
        // it is delete
        if (m_genomePosition[u]+1<m_genomePosition[v]){
            for (auto pos=m_genomePosition[u]+1; pos<m_genomePosition[v]; pos++){
                string ref = m_labels[pos];
                string alt = "-";
                // allele
                GenericAllele allele;
                allele.m_chrID       = 0;
                allele.m_chrPosition = pos-1;
                allele.m_allele      = alt;
                allele.m_reference   = ref;
                // save
                pathVariants[pos-1]  = allele;
            }
        }

        // update
        u = v;
    }
}
