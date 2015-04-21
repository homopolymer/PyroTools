#include "GenericProbabilisticAlignment.h"
#include "GenericSequenceGlobal.h"
#include "GenericFastaTools.h"
#include "GenericBamAlignmentTools.h"
using namespace GenericSequenceTools;

#include <limits>
#include <cmath>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <set>
#include <iomanip>
using namespace std;

#include "nlopt.hpp"

//---------------------------------------
struct nlopt_aux_t
{
    vector<string>* readSeqs;
    vector<string>* genomeSeqs;
    int numTrainData;
    double band;
    int iter;
    int verbosity;
    GenericProbabilisticAlignment* ProbAligner;

    map<int, tuple<long double, long double>> zeroProbData;

    vector<string> badReads;
    vector<string> badGenomes;
    set<int> badDataIndex;

    double optim_loglik;
    vector<double> optim_w;
};

double norm2(const vector<double>& w)
{
    double n2 = 0;
    for (int i=0; i<w.size(); ++i)
    {
        n2 += w[i]*w[i];
    }

    return (n2);
}

double loglik(const vector<double>& w, vector<double>& grad, void* aux_data)
{
    nlopt_aux_t* settings = (nlopt_aux_t*) aux_data;

    double lambda = 0.1;

    // print message
    {
        cout << "[Probabilistic Aligner Training #" << settings->iter+1 << " iteration]" << endl;
    }

    // --------------------------------------------------------------------------
    // dimension of w
    int n = 0;
    int shift;

    // set parameters of ProbAligner

    // transition feature
    shift = n;
    for (int i=0; i<settings->ProbAligner->NumberTransitionFeatures; ++i, ++n)
    {
        ProbalnTransitionFeature f = settings->ProbAligner->TransitionFeatureSet[i];
        settings->ProbAligner->TransitionFeatureWeight[f] = w[i+shift];
    }

    // read homopolymer feature    
    shift = n;
    for (int i=0; i<settings->ProbAligner->NumberReadHomopolymerFeatures; ++i, ++n)
    {
        ProbalnReadHomopolymerFeature f = settings->ProbAligner->ReadHomopolymerFeatureSet[i];
        settings->ProbAligner->ReadHomopolymerFeatureWeight[f] = w[i+shift];
    }

    // genome homopolymer feature
    shift = n;
    for (int i=0; i<settings->ProbAligner->NumberGenomeHomopolymerFeatures; ++i, ++n)
    {
        ProbalnGenomeHomopolymerFeature f = settings->ProbAligner->GenomeHomopolymerFeatureSet[i];
        settings->ProbAligner->GenomeHomopolymerFeatureWeight[f] = w[i+shift];
    }

    // neighbor context feature
    if (USE_MARKOV_FEATURE)
    {
        shift = n;
        for (int i=0; i<settings->ProbAligner->NumberNeighborContextFeatures; ++i, ++n)
        {
            ProbalnNeighborContextFeature f = settings->ProbAligner->NeighborContextFeatureSet[i];
            settings->ProbAligner->NeighborContextFeatureWeight[f] = w[i+shift];
        }
    }

    // emission feature
    shift = n;
    for (int i=0; i<settings->ProbAligner->NumberEmissionFeatures; ++i, ++n)
    {
        ProbalnEmissionFeature f = settings->ProbAligner->EmissionFeatureSet[i];
        settings->ProbAligner->EmissionFeatureWeight[f] = w[i+shift];
    }

    // -------------------------------------------------------------------------
    // compute feature counts and expectations
    VectorDouble TransitionFeatureCounts(settings->ProbAligner->NumberTransitionFeatures, 0);
    VectorDouble ReadHomopolymerFeatureCounts(settings->ProbAligner->NumberReadHomopolymerFeatures, 0);
    VectorDouble GenomeHomopolymerFeatureCounts(settings->ProbAligner->NumberGenomeHomopolymerFeatures, 0);
    VectorDouble NeighborContextFeatureCounts(settings->ProbAligner->NumberNeighborContextFeatures, 0);
    VectorDouble EmissionFeatureCounts(settings->ProbAligner->NumberEmissionFeatures, 0);

    VectorDouble TransitionFeatureExpectations(settings->ProbAligner->NumberTransitionFeatures, 0);
    VectorDouble ReadHomopolymerFeatureExpectations(settings->ProbAligner->NumberReadHomopolymerFeatures, 0);
    VectorDouble GenomeHomopolymerFeatureExpectations(settings->ProbAligner->NumberGenomeHomopolymerFeatures, 0);
    VectorDouble NeighborContextFeatureExpectations(settings->ProbAligner->NumberNeighborContextFeatures, 0);
    VectorDouble EmissionFeatureExpectations(settings->ProbAligner->NumberEmissionFeatures, 0);


    int numGapInData = 0;
    int numDataWithGap = 0;

    long double logProbSum = 0;

    // debug
    vector<string> badRead, badGenome, badAlnRead, badAlnGenome;

    // iterate over training data
    for (int i=0; i<settings->numTrainData; ++i)
    {

        // skip if data probability was zero
        if (settings->zeroProbData.find(i)!=settings->zeroProbData.end())
            continue;

        string alignRead, alignGenome;
        long double viterbiScore, forwardScore;
        VectorDouble pfc, pfe;
        VectorDouble rhfc, rhfe;
        VectorDouble ghfc, ghfe;
        VectorDouble lcfc, lcfe;
        VectorDouble rcfc, rcfe;
        VectorDouble efc, efe;

        // compute the alignment
        long double score = settings->ProbAligner->ViterbiComputation(
                        (*settings->readSeqs)[i], (*settings->genomeSeqs)[i], settings->band,
                        alignRead, alignGenome, viterbiScore, forwardScore,
                        pfc, rhfc, ghfc, lcfc, efc,
                        pfe, rhfe, ghfe, lcfe, efe);

        // compute the number of gap in the alignment
        numGapInData += GenericBamAlignmentTools::numGap(alignRead, alignGenome);
        numDataWithGap += 1;

        if (score<=0)
        {
            settings->zeroProbData[i] = tuple<long double,long double>(viterbiScore, forwardScore);
            continue;
        }

        logProbSum += log(score);

        // print message
        {

            cout << "\rData:"
                 << i+1 << "/" << settings->numTrainData << " "
                 << "ViterbiScore:" << setw(10) << setprecision(6) << viterbiScore << "   "
                 << "ForwardScore:" << setw(10) << setprecision(6) << log(forwardScore) << "   "
                 << "Score:" << setw(10) << setprecision(6) << score << " "
                 << "   " << flush;
        }

        // debug
        if (score<0.1)
        {
            badRead.push_back((*settings->readSeqs)[i]);
            badGenome.push_back((*settings->genomeSeqs)[i]);
            badAlnRead.push_back(alignRead);
            badAlnGenome.push_back(alignGenome);

            if (settings->badDataIndex.find(i)==settings->badDataIndex.end())
            {
                settings->badReads.push_back((*settings->readSeqs)[i]);
                settings->badGenomes.push_back((*settings->genomeSeqs)[i]);
                settings->badDataIndex.insert(i);
            }
        }

        // position feature
        for (int j=0; j<settings->ProbAligner->NumberTransitionFeatures; ++j)
        {
            TransitionFeatureCounts[j]       += pfc[j];
            TransitionFeatureExpectations[j] += pfe[j];
        }

        // read homopolymer feature
        for (int j=0; j<settings->ProbAligner->NumberReadHomopolymerFeatures; ++j)
        {
            ReadHomopolymerFeatureCounts[j]       += rhfc[j];
            ReadHomopolymerFeatureExpectations[j] += rhfe[j];
        }

        // genome homopolymer feature
        for (int j=0; j<settings->ProbAligner->NumberGenomeHomopolymerFeatures; ++j)
        {
            GenomeHomopolymerFeatureCounts[j]       += ghfc[j];
            GenomeHomopolymerFeatureExpectations[j] += ghfe[j];
        }

        // neighbor context feature
        if (USE_MARKOV_FEATURE)
        {
            for (int j=0; j<settings->ProbAligner->NumberNeighborContextFeatures; ++j)
            {
                NeighborContextFeatureCounts[j]       += lcfc[j];
                NeighborContextFeatureExpectations[j] += lcfe[j];
            }
        }

        // emission feature
        for (int j=0; j<settings->ProbAligner->NumberEmissionFeatures; ++j)
        {
            EmissionFeatureCounts[j]       += efc[j];
            EmissionFeatureExpectations[j] += efe[j];
        }
    }

    // -------------------------------------
    // grad
    if (!grad.empty())
    {
        n = 0;

        // position feature
        shift  = n;
        for (int i=0; i<settings->ProbAligner->NumberTransitionFeatures; ++i, ++n)
        {
            grad[i+shift] = TransitionFeatureCounts[i] - TransitionFeatureExpectations[i];
        }

        // read homopolymer feature
        shift = n;
        for (int i=0; i<settings->ProbAligner->NumberReadHomopolymerFeatures; ++i, ++n)
        {
            grad[i+shift] = ReadHomopolymerFeatureCounts[i] - ReadHomopolymerFeatureExpectations[i];
        }

        // genome homopolymer feature
        shift = n;
        for (int i=0; i<settings->ProbAligner->NumberGenomeHomopolymerFeatures; ++i, ++n)
        {
            grad[i+shift] = GenomeHomopolymerFeatureCounts[i] - GenomeHomopolymerFeatureExpectations[i];
        }

        // neighbor context feature
        if (USE_MARKOV_FEATURE)
        {
            shift = n;
            for (int i=0; i<settings->ProbAligner->NumberNeighborContextFeatures; ++i, ++n)
            {
                grad[i+shift] = NeighborContextFeatureCounts[i] - NeighborContextFeatureExpectations[i];
            }
        }

        // emission feature
        shift = n;
        for (int i=0; i<settings->ProbAligner->NumberEmissionFeatures; ++i, ++n)
        {
            grad[i+shift] = EmissionFeatureCounts[i] - EmissionFeatureExpectations[i];
        }


        // regularize
        for (int i=0; i<w.size(); ++i)
        {
            grad[i] -= lambda*w[i];
        }
    }


    // regularize
    logProbSum -= 0.5*lambda*norm2(w);
    // penalize on the number of gap
    logProbSum -= numGapInData/double(numDataWithGap);

    // print message
    {
        cout << "log-likelihood-score=" << setprecision(numeric_limits<double>::digits10+1) << logProbSum << endl;

        for (map<int, tuple<long double, long double>>::iterator iter=settings->zeroProbData.begin();
             iter!=settings->zeroProbData.end(); ++iter)
        {
            cout << "Datum:"
                 << iter->first << ", "
                 << "ViterbiScore:"
                 << get<0>(iter->second) << ", "
                 << "ForwardScore:"
                 << get<1>(iter->second) << endl;
        }
    }
    if (!grad.empty() && settings->verbosity>=2)
    {
        int n=0,shift;
        // position feature
        shift = n;
        for (int i=0; i<settings->ProbAligner->NumberTransitionFeatures; ++i, ++n)
        {
            ProbalnTransitionFeature f = settings->ProbAligner->TransitionFeatureSet[i];
            cout << settings->ProbAligner->TransitionFeatureLabel[f] << " "
                 << setprecision(5)
                 << TransitionFeatureCounts[i] << " "
                 << TransitionFeatureExpectations[i] << " "
                 << grad[i+shift] << " "
                 << settings->ProbAligner->TransitionFeatureWeight[f]
                 << endl;
        }

        // read homopolymer feature
        shift = n;
        for (int i=0; i<settings->ProbAligner->NumberReadHomopolymerFeatures; ++i, ++n)
        {
            ProbalnReadHomopolymerFeature f = settings->ProbAligner->ReadHomopolymerFeatureSet[i];
            cout << settings->ProbAligner->ReadHomopolymerFeatureLabel[f] << " "
                 << setprecision(5)
                 << ReadHomopolymerFeatureCounts[i] << " "
                 << ReadHomopolymerFeatureExpectations[i] << " "
                 << grad[i+shift] << " "
                 << settings->ProbAligner->ReadHomopolymerFeatureWeight[f]
                    << endl;
        }


        // genome homopolymer feature
        shift = n;
        for (int i=0; i<settings->ProbAligner->NumberGenomeHomopolymerFeatures; ++i, ++n)
        {
            ProbalnGenomeHomopolymerFeature f = settings->ProbAligner->GenomeHomopolymerFeatureSet[i];
            cout << settings->ProbAligner->GenomeHomopolymerFeatureLabel[f] << " "
                 << setprecision(5)
                 << GenomeHomopolymerFeatureCounts[i] << " "
                 << GenomeHomopolymerFeatureExpectations[i] << " "
                 << grad[i+shift] << " "
                 << settings->ProbAligner->GenomeHomopolymerFeatureWeight[f]
                    << endl;
        }


        // neighbor context feature
        if (USE_MARKOV_FEATURE)
        {
            shift = n;
            for (int i=0; i<settings->ProbAligner->NumberNeighborContextFeatures; ++i, ++n)
            {
                ProbalnNeighborContextFeature f = settings->ProbAligner->NeighborContextFeatureSet[i];
                cout << settings->ProbAligner->NeighborContextFeatureLabel[f] << " "
                     << setprecision(5)
                     << NeighborContextFeatureCounts[i] << " "
                     << NeighborContextFeatureExpectations[i] << " "
                     << grad[i+shift] << " "
                     << settings->ProbAligner->NeighborContextFeatureWeight[f]
                     << endl;
            }
        }

        // emission feature
        shift = n;
        for (int i=0; i<settings->ProbAligner->NumberEmissionFeatures; ++i, ++n)
        {
            ProbalnEmissionFeature f = settings->ProbAligner->EmissionFeatureSet[i];
            cout << settings->ProbAligner->EmissionFeatureLabel[f] << " "
                 << setprecision(5)
                 << EmissionFeatureCounts[i] << " "
                 << EmissionFeatureExpectations[i] << " "
                 << grad[i+shift] << " "
                 << settings->ProbAligner->EmissionFeatureWeight[f]
                 << endl;
        }
    }

    // debug
    if (settings->verbosity>=3)
    {
        for (int i=0; i<settings->badReads.size(); ++i)
        {
            long double t_score;
            string t_alnRead, t_alnGenome;
            t_score = settings->ProbAligner->ViterbiComputation(settings->badReads[i],
                                                                settings->badGenomes[i],
                                                                settings->band,
                                                                t_alnRead, t_alnGenome);

            BamMD t_md;
            GenericBamAlignmentTools::calculateMD(t_alnRead, t_alnGenome, t_md);

            cout << ">" << i+1 << " " << t_score << " " << t_md << endl;
            cout << t_alnRead << endl;
            cout << t_alnGenome << endl;
        }
    }

    settings->iter++;

    // optimal result so far
    if (settings->optim_loglik>logProbSum)
    {
        settings->optim_loglik = logProbSum;
        settings->optim_w.assign(w.begin(), w.end());
    }

    return logProbSum;
}


//-----------------------------------------------------------

GenericProbabilisticAlignment::GenericProbabilisticAlignment()
    : NumberTransitionFeatures(0)
    , NumberReadHomopolymerFeatures(0)
    , NumberGenomeHomopolymerFeatures(0)
    , NumberNeighborContextFeatures(0)
    , NumberEmissionFeatures(0)
{

    // --------------------------------------------------
    // transition features
    initTransitionFeature();

    // --------------------------------------------------
    // read homopolymer features
    initReadHomopolymerFeature();

    // --------------------------------------------------
    // genome homopolymer features
    initGenomeHomopolymerFeature();

    if (USE_MARKOV_FEATURE)
    {
        // --------------------------------------------------
        // neighbor context feature
        initNeighborContextFeature();
    }

    // --------------------------------------------------
    // emission feature
    initEmissionFeature();
}


// -------------------------------------------------------
// featue initialization

// transition feature
void GenericProbabilisticAlignment::initTransitionFeature()
{

    ProbalnState previousState;
    ProbalnState currentState;

    // --------------------------------------------------
    // transition features

    // begin hidden state
    for (currentState=PROBALN_MATCH; currentState<=PROBALN_INSERT; ++currentState)
    {
        setTransitionFeature(PROBALN_BEGIN, currentState);
    }

    // transition from previous to current
    for (previousState=PROBALN_MATCH; previousState<=PROBALN_DELETE; ++previousState)
    {
        for (currentState=PROBALN_MATCH; currentState<=PROBALN_DELETE; ++currentState)
        {
            setTransitionFeature(previousState, currentState);
        }
    }

    // end hidden state
    for (previousState=PROBALN_MATCH; previousState<=PROBALN_DELETE; ++previousState)
    {
        setTransitionFeature(previousState, PROBALN_END);
    }
}

// read homopolymer feature
void GenericProbabilisticAlignment::initReadHomopolymerFeature()
{

    ProbalnState currentState = PROBALN_INSERT;
    int homopolymerLength;

    for (homopolymerLength=midHomopolymerSize; homopolymerLength<=maxHomopolymerSize; ++homopolymerLength)
    {
        for (int delta=1; delta<homopolymerLength; ++delta)
        {
            if (delta>5)
                continue;

            setReadHomopolymerFeature(currentState, homopolymerLength, delta);
        }
    }
}

// genome homopolymer feature
void GenericProbabilisticAlignment::initGenomeHomopolymerFeature()
{

    ProbalnState currentState = PROBALN_DELETE;
    int homopolymerLength;

    for (homopolymerLength=midHomopolymerSize; homopolymerLength<=maxHomopolymerSize; ++homopolymerLength)
    {
        for (int delta=1; delta<homopolymerLength; ++delta)
        {
            if (delta>5)
                continue;

            setGenomeHomopolymerFeature(currentState, homopolymerLength, delta);
        }
    }
}

// neighbor context feature
void GenericProbabilisticAlignment::initNeighborContextFeature()
{

    ProbalnState currentState = PROBALN_INSERT;
    int homopolymerLength;

    for (homopolymerLength=midHomopolymerSize; homopolymerLength<=maxHomopolymerSize; ++homopolymerLength)
    {
        setNeighborContextFeature(currentState, homopolymerLength);
    }
}

// emission feature
void GenericProbabilisticAlignment::initEmissionFeature()
{

    ProbalnState currentState;
    SequenceBase referenceAllele;
    SequenceBase readAllele;

    // match
    currentState = PROBALN_MATCH;
    for (referenceAllele=ADENINE; referenceAllele<=THYMINE; ++referenceAllele)
    {
        setEmissionFeature(currentState, referenceAllele, referenceAllele);
    }

    // mismatch
    currentState = PROBALN_MISMATCH;
    for (referenceAllele=ADENINE; referenceAllele<=THYMINE; ++referenceAllele)
    {
        for (readAllele=ADENINE; readAllele<=THYMINE; ++readAllele)
        {
            if (referenceAllele==readAllele)
                continue;

            setEmissionFeature(currentState, referenceAllele, readAllele);
        }
    }

    // insert
    currentState = PROBALN_INSERT;
    referenceAllele = SPACE;
    for (readAllele=ADENINE; readAllele<=THYMINE; ++readAllele)
    {
        setEmissionFeature(currentState, referenceAllele, readAllele);
    }

    // delete
    currentState = PROBALN_DELETE;
    readAllele = SPACE;
    for (referenceAllele=ADENINE; referenceAllele<=THYMINE; ++referenceAllele)
    {
        setEmissionFeature(currentState, referenceAllele, readAllele);
    }
}

// ------------------------------------------------------
// feature setting

// transition feature
void GenericProbabilisticAlignment::setTransitionFeature(ProbalnState previousState, ProbalnState currentState)
{
    // construct feature
    ProbalnTransitionFeature f(previousState, currentState);

    // insert to the set
    TransitionFeatureSet.push_back(f);

    // set feature index
    TransitionFeatureIndex[f] = NumberTransitionFeatures;

    // set feature name
    stringstream strBuf;

    strBuf << "TransitionFeature "    << NumberTransitionFeatures++        << " ";
    strBuf << "PreviousState:"        << ProbalnStateName[previousState]   << " ";
    strBuf << "CurrentState:"         << ProbalnStateName[currentState];

    TransitionFeatureLabel[f] = strBuf.str();

    // set feature weight
    TransitionFeatureWeight[f] = 0;
}


// read homopolymer feature
void GenericProbabilisticAlignment::setReadHomopolymerFeature(ProbalnState currentState, int homopolymerLength, int delta)
{
    // construct feature
    ProbalnReadHomopolymerFeature f(currentState, homopolymerLength, delta);

    // insert to the set
    ReadHomopolymerFeatureSet.push_back(f);

    // set feature index
    ReadHomopolymerFeatureIndex[f] = NumberReadHomopolymerFeatures;

    // set feature name
    stringstream strBuf;

    strBuf << "ReadHomopolymerFeature "    << NumberReadHomopolymerFeatures++ << " ";
    strBuf << "CurrentState:"              << ProbalnStateName[currentState]  << " ";
    strBuf << "Length:"                    << homopolymerLength               << " ";
    strBuf << "Delta:"                     << delta;

    ReadHomopolymerFeatureLabel[f] = strBuf.str();

    // set feature weight
    ReadHomopolymerFeatureWeight[f] = 0;
}

// genome homopolymer feature
void GenericProbabilisticAlignment::setGenomeHomopolymerFeature(ProbalnState currentState, int homopolymerLength, int delta)
{
    // construct feature
    ProbalnGenomeHomopolymerFeature f(currentState, homopolymerLength, delta);

    // insert to the set
    GenomeHomopolymerFeatureSet.push_back(f);

    // set feature index
    GenomeHomopolymerFeatureIndex[f] = NumberGenomeHomopolymerFeatures;

    // set feature name
    stringstream strBuf;

    strBuf << "GenomeHomopolymerFeature "    << NumberGenomeHomopolymerFeatures++       << " ";
    strBuf << "CurrentState:"                << ProbalnStateName[currentState]          << " ";
    strBuf << "Length:"                      << homopolymerLength                       << " ";
    strBuf << "Delta:"                       << delta;

    GenomeHomopolymerFeatureLabel[f] = strBuf.str();

    // set feature weight
    GenomeHomopolymerFeatureWeight[f] = 0;

}

// neighbor context feature
void GenericProbabilisticAlignment::setNeighborContextFeature(ProbalnState currentState, int homopolymerLength)
{
    // construct feature
    ProbalnNeighborContextFeature f(currentState, homopolymerLength);

    // insert to the set
    NeighborContextFeatureSet.push_back(f);

    // set feature index
    NeighborContextFeatureIndex[f] = NumberNeighborContextFeatures;

    // set feature name
    stringstream strBuf;

    strBuf << "NeighborContextFeature "    << NumberNeighborContextFeatures++    << " ";
    strBuf << "CurrentState:"              << ProbalnStateName[currentState] << " ";
    strBuf << "Length:"                    << homopolymerLength;

    NeighborContextFeatureLabel[f] = strBuf.str();

    // set feature weight
    NeighborContextFeatureWeight[f] = 0;
}


// emission feature
void GenericProbabilisticAlignment::setEmissionFeature(ProbalnState currentState, SequenceBase genomeAllele, SequenceBase readAllele)
{
    // construct feature
    ProbalnEmissionFeature f(currentState, genomeAllele, readAllele);

    // insert to the set
    EmissionFeatureSet.push_back(f);

    // set feature index
    EmissionFeatureIndex[f] = NumberEmissionFeatures;

    // set feature name
    stringstream strBuf;

    strBuf << "EmissionFeature "    << NumberEmissionFeatures++       << " ";
    strBuf << "CurrentState:"       << ProbalnStateName[currentState] << " ";
    strBuf << "GenomeAllele:"       << AlphabetChar[genomeAllele]     << " ";
    strBuf << "ReadAllele:"         << AlphabetChar[readAllele];

    EmissionFeatureLabel[f] = strBuf.str();

    // set feature weight
    EmissionFeatureWeight[f] = 0;
}

// ------------------------------------------------------------------------
// Viterbi computation

void getFeature(string& readSeq, string& genomeSeq,
                int readPointer, int genomePointer,
                ProbalnState previousState, ProbalnState currentState,
                vector<int>& readSuccedeHomopolymerLength,
                vector<int>& genomeSuccedeHomopolymerLength,
                vector<int>& leftNeighborHomopolymerLength,
                vector<int>& rightNeighborHomopolymerLength,
                ProbalnTransitionFeature&        transitionFeature,
                ProbalnReadHomopolymerFeature&   readHomopolymerFeature,
                ProbalnGenomeHomopolymerFeature& genomeHomopolymerFeature,
                ProbalnNeighborContextFeature&   neighborContextFeature,
                ProbalnEmissionFeature&          emissionFeature
                )
{
    // ----------------------------------------
    // transition feature
    transitionFeature = ProbalnTransitionFeature(previousState, currentState);

    // ----------------------------------------
    // read homopolymer feature

    readHomopolymerFeature = ProbalnReadHomopolymerFeature();
    if (currentState==PROBALN_INSERT && genomePointer+1<genomeSeq.length())
    {
        if (readSeq[readPointer]==genomeSeq[genomePointer+1])
        {
            int rhl, thl, delta;
            rhl = readSuccedeHomopolymerLength[readPointer]+1;
            thl = genomeSuccedeHomopolymerLength[genomePointer+1]+1;

            delta = rhl-thl;

            if (delta>5)
                delta = 5;

            if (rhl>=midHomopolymerSize && delta>=1 && delta<=5)
                readHomopolymerFeature = ProbalnReadHomopolymerFeature(PROBALN_INSERT, rhl, delta);
        }
    }

    // ----------------------------------------
    // genome homopolymer feature

    genomeHomopolymerFeature = ProbalnGenomeHomopolymerFeature();
    if (currentState==PROBALN_DELETE && readPointer+1<readSeq.length())
    {
        if (readSeq[readPointer+1]==genomeSeq[genomePointer])
        {
            int rhl, thl, delta;
            rhl = readSuccedeHomopolymerLength[readPointer+1]+1;
            thl = genomeSuccedeHomopolymerLength[genomePointer]+1;
            delta = thl-rhl;;

            if (delta>5)
                delta=5;

            if (thl>=midHomopolymerSize && delta>=1 && delta<=5)
                genomeHomopolymerFeature = ProbalnGenomeHomopolymerFeature(PROBALN_DELETE, thl, delta);
        }
    }

    // ----------------------------------------
    // neighbor context feature

    if (currentState==PROBALN_INSERT && USE_MARKOV_FEATURE &&
            readSeq[readPointer]!=genomeSeq[genomePointer] &&
            (readSeq[readPointer]!=genomeSeq[genomePointer+1] && genomePointer+1<genomeSeq.length()))
    {
        int len = leftNeighborHomopolymerLength[readPointer];
        if (len<rightNeighborHomopolymerLength[readPointer])
            len = rightNeighborHomopolymerLength[readPointer];

        neighborContextFeature = ProbalnNeighborContextFeature(currentState, len);
    }else
    {
        neighborContextFeature = ProbalnNeighborContextFeature();
    }


    // ----------------------------------------
    // emission feature

    if (currentState==PROBALN_MATCH)
    {
        if (Alpha2Base[readSeq[readPointer]]!=AMBIGUOUS && Alpha2Base[genomeSeq[genomePointer]]!=AMBIGUOUS)
        {
            emissionFeature = ProbalnEmissionFeature(currentState,
                                                     Alpha2Base[genomeSeq[genomePointer]],
                                                     Alpha2Base[readSeq[readPointer]]);
        }else if (Alpha2Base[readSeq[readPointer]]==AMBIGUOUS)
        {
            emissionFeature = ProbalnEmissionFeature(currentState,
                                                     Alpha2Base[genomeSeq[genomePointer]],
                                                     Alpha2Base[genomeSeq[genomePointer]]);
        }else if (Alpha2Base[genomeSeq[genomePointer]]==AMBIGUOUS)
        {
            emissionFeature = ProbalnEmissionFeature(currentState,
                                                     Alpha2Base[readSeq[readPointer]],
                                                     Alpha2Base[readSeq[readPointer]]);
        }
    }

    if (currentState==PROBALN_MISMATCH)
    {
        emissionFeature = ProbalnEmissionFeature(currentState,
                                                 Alpha2Base[genomeSeq[genomePointer]],
                                                 Alpha2Base[readSeq[readPointer]]);
    }

    if (currentState==PROBALN_INSERT)
    {
        if (Alpha2Base[readSeq[readPointer]]!=AMBIGUOUS)
        {
            emissionFeature = ProbalnEmissionFeature(currentState,
                                                     SPACE,
                                                     Alpha2Base[readSeq[readPointer]]);
        }else
        {
            emissionFeature = ProbalnEmissionFeature(currentState,
                                                     SPACE,
                                                     ADENINE);
        }
    }

    if (currentState==PROBALN_DELETE)
    {
        if (Alpha2Base[genomeSeq[genomePointer]]!=AMBIGUOUS)
        {
            emissionFeature = ProbalnEmissionFeature(currentState,
                                                     Alpha2Base[genomeSeq[genomePointer]],
                                                     SPACE);
        }else
        {
            emissionFeature = ProbalnEmissionFeature(currentState,
                                                     ADENINE,
                                                     SPACE);
        }
    }
}


int getGenomeStartPosition(string& readSeq, string& genomeSeq)
{

    return 0;

    int genomeStartPosition;
    int genomeSeqLength = genomeSeq.length();

    // find the start position in genome
    vector<int> genomeHitCount(genomeSeqLength, 0);
    for (int i=0; i<10; i++)
    {
        for (int k=1; k<=10; k++)
        {
            if (i+k>=readSeq.length())
                continue;

            string pattern = readSeq.substr(i, k);
            int pos = genomeSeq.find(pattern);
            if (pos!=string::npos)
                genomeHitCount[pos] ++;
        }
    }

    genomeStartPosition = 0;
    int maxHitCount = 0;
    for (int i=0; i<genomeHitCount.size(); ++i)
    {
        if (maxHitCount < genomeHitCount[i])
        {
            maxHitCount = genomeHitCount[i];
            genomeStartPosition = i;
        }
    }


    return genomeStartPosition;
}

void setViterbiBoundaries(string& readSeq, string& genomeSeq, double band, vector<int>& LeftBoundary, vector<int>& RightBoundary)
{
    int genomeStartPosition;
    int genomeSeqLength = genomeSeq.length();
    int readSeqLength   = readSeq.length();
    int deltaLength = abs(genomeSeqLength-readSeqLength);

    // find the start position in genome
    genomeStartPosition = getGenomeStartPosition(readSeq, genomeSeq);
    genomeStartPosition = 0;

    // initialize left boundary and right boundary
    LeftBoundary.assign(genomeSeqLength, 0);
    RightBoundary.assign(genomeSeqLength, readSeqLength);

    // update left boundary and right boundary if band < 1.0
    if (band<1.0)
    {
        int windowSize = band*readSeqLength;

        int halfWindowSize = ceil(0.5*windowSize);
        if (halfWindowSize<=deltaLength)
            halfWindowSize=deltaLength+3;

        int left  = -genomeStartPosition - halfWindowSize;
        int right = halfWindowSize - genomeStartPosition;

        // iterate from first position to last position in genome
        for (int i=0; i<genomeSeqLength; ++i, ++left, ++right)
        {
            if (left>=0)
            {
                if (left<readSeqLength)
                    LeftBoundary[i] = left;
                else
                    LeftBoundary[i] = readSeqLength;
            }
            else
                LeftBoundary[i] = 0;

            if (right>=0)
            {
                if (right<readSeqLength)
                    RightBoundary[i] = right;
                else
                    RightBoundary[i] = readSeqLength;
            }
            else
                RightBoundary[i] = 0;
        }
    }
}

// transition feature
long double transitionFeatureWeight(ProbalnTransitionFeature& f, ProbalnTransitionFeatureValue& W)
{
    ProbalnTransitionFeatureValue::iterator iter = W.find(f);

    if (iter==W.end())
    {
        cerr << "GenericProbabilisticAlignment ERROR: undefined TransitionFeature"
             << "("
             << ProbalnStateName[get<0>(f)] << ","
             << ProbalnStateName[get<1>(f)]
             << ")"
             << "... Aborting" << endl;
        exit(0);
    }

    return iter->second;
}

// read homopolymer feature
long double readHomopolymerFeatureWeight(ProbalnReadHomopolymerFeature& f, ProbalnReadHomopolymerFeatureValue& W)
{
    ProbalnReadHomopolymerFeatureValue::iterator iter = W.find(f);

    if (iter==W.end())
    {
        cerr << "GenericProbabilisticAlignment ERROR: undefined ReadHomopolymerFeature"
             << "("
             << ProbalnStateName[get<0>(f)] << ","
             << get<1>(f) << ","
             << get<2>(f)
             << ")"
             << "... Aborting" << endl;
        exit(0);
    }

    return iter->second;
}

// genome homopolymer feature
long double genomeHomopolymerFeatureWeight(ProbalnGenomeHomopolymerFeature& f, ProbalnGenomeHomopolymerFeatureValue& W)
{
    ProbalnGenomeHomopolymerFeatureValue::iterator iter = W.find(f);

    if (iter==W.end())
    {
        cerr << "GenericProbabilisticAlignment ERROR: undefined GenomeHomopolymerFeature"
             << "("
             << ProbalnStateName[get<0>(f)] << ","
             << get<1>(f) << ","
             << get<2>(f)
             << ")"
             << "... Aborting" << endl;
        exit(0);
    }

    return iter->second;
}

// neighbor context feature
long double neighborContextFeatureWeight(ProbalnNeighborContextFeature& f, ProbalnNeighborContextFeatureValue& W)
{
    ProbalnNeighborContextFeatureValue::iterator iter = W.find(f);

    if (iter==W.end())
    {
        cerr << "GenericProbabilisticAlignment ERROR: undefined NeighborContextFeature"
             << "("
             << ProbalnStateName[get<0>(f)] << ","
             << get<1>(f)
             << ")"
             << "... Aborting" << endl;
        exit(0);
    }

    return iter->second;
}


// emission feature
long double emissionFeatureWeight(ProbalnEmissionFeature& f, ProbalnEmissionFeatureValue& W)
{
    ProbalnEmissionFeatureValue::iterator iter = W.find(f);

    if (iter==W.end())
    {
        cerr << "GenericProbabilisticAlignment ERROR: undefined EmissionFeature" << " "
             << "("
             << ProbalnStateName[get<0>(f)] << "," << AlphabetChar[get<1>(f)] << "," << AlphabetChar[get<2>(f)]
             << ")"
             <<"... Aborting" << endl;
        exit(0);
    }

    return iter->second;
}

long double GenericProbabilisticAlignment::featureWeightSum(ProbalnTransitionFeature &TransitionFeature, ProbalnReadHomopolymerFeature &ReadHomopolymerFeature, ProbalnGenomeHomopolymerFeature &GenomeHomopolymerFeature, ProbalnNeighborContextFeature &NeighborContextFeature, ProbalnEmissionFeature &EmissionFeature)
{
    long double score = 0;

    // transition feature
    if (TransitionFeature!=ProbalnTransitionFeature())
        score += transitionFeatureWeight(TransitionFeature, TransitionFeatureWeight);

    // read homopolymer feature
    if (ReadHomopolymerFeature!=ProbalnReadHomopolymerFeature())
    {
        ProbalnReadHomopolymerFeature f(ReadHomopolymerFeature);

        if (get<1>(f)>=midHomopolymerSize && get<2>(f)>=1)
        {
            if (get<1>(f)>maxHomopolymerSize)
            {
                get<1>(f) = maxHomopolymerSize;
            }
            if (get<2>(f)>5)
            {
                get<2>(f) = 5;
            }
            score += readHomopolymerFeatureWeight(f, ReadHomopolymerFeatureWeight);

//            if (get<1>(ReadHomopolymerFeature)>maxHomopolymerSize)
//                score += get<1>(ReadHomopolymerFeature)-maxHomopolymerSize;
        }
    }

    // genome homopolymer feature
    if (GenomeHomopolymerFeature!=ProbalnGenomeHomopolymerFeature())
    {
        ProbalnGenomeHomopolymerFeature f(GenomeHomopolymerFeature);

        if (get<1>(f)>=midHomopolymerSize && get<2>(f)>=1)
        {
            if (get<1>(f)>maxHomopolymerSize)
            {
                get<1>(f) = maxHomopolymerSize;
            }
            if (get<2>(f)>5)
            {
                get<2>(f) = 5;
            }
            score += genomeHomopolymerFeatureWeight(f, GenomeHomopolymerFeatureWeight);

//            if (get<1>(GenomeHomopolymerFeature)>maxHomopolymerSize)
//                score += get<1>(GenomeHomopolymerFeature)-maxHomopolymerSize;
        }
    }

    // neighbor context feature
    if (NeighborContextFeature!=ProbalnNeighborContextFeature() && USE_MARKOV_FEATURE)
    {
        ProbalnNeighborContextFeature f(NeighborContextFeature);

        if (get<1>(f)>=midHomopolymerSize)
        {
            if (get<1>(f)>maxHomopolymerSize)
            {
                score += get<1>(f)-maxHomopolymerSize;
                get<1>(f) = maxHomopolymerSize;
            }
            score += neighborContextFeatureWeight(f, NeighborContextFeatureWeight);
        }
    }


    // emission feature
    if (EmissionFeature!=ProbalnEmissionFeature())
        score += emissionFeatureWeight(EmissionFeature, EmissionFeatureWeight);

    return score;
}

long double GenericProbabilisticAlignment::ViterbiComputation(string &readSeq, string &genomeSeq, double band,

                                                              string &alignReadSeq, string &alignGenomeSeq,

                                                              long double &ViterbiScore,
                                                              long double &ForwardScore,

                                                              vector<long double> &TransitionFeatureCounts,
                                                              vector<long double> &ReadHomopolymerFeatureCounts,
                                                              vector<long double> &GenomeHomopolymerFeatureCounts,
                                                              vector<long double> &NeighborContextFeatureCounts,
                                                              vector<long double> &EmissionFeatureCounts,

                                                              vector<long double> &TransitionFeatureExpectations,
                                                              vector<long double> &ReadHomopolymerFeatureExpectations,
                                                              vector<long double> &GenomeHomopolymerFeatureExpectations,
                                                              vector<long double> &NeighborContextFeatureExpectations,
                                                              vector<long double> &EmissionFeatureExpectations)
{

    // alignment probability
    long double ProbabilisticScore = 0.0;

    // genome length
    int genomeSeqSize = genomeSeq.length();

    // read length
    int readSeqSize = readSeq.length();

    // left boundary and right boundary
    VectorInteger leftBoundary;
    VectorInteger rightBoundary;

    setViterbiBoundaries(readSeq, genomeSeq, band, leftBoundary, rightBoundary);

    // Matrices and Vectors related to Viterbi algorithm

    Matrix3Double  ViterbiMatrix;
    setMatrixValue(ViterbiMatrix,
                   genomeSeqSize, readSeqSize, PROBALN_STATE_SPACE_SIZE,
                   DOUBLE_NEGATIVE_INFINITY);
    Matrix3Integer ViterbiPreviousStateMatrix;
    setMatrixValue(ViterbiPreviousStateMatrix,
                   genomeSeqSize, readSeqSize, PROBALN_STATE_SPACE_SIZE,
                   PROBALN_UNDEFINE);

    VectorDouble   ViterbiEndArray;
    setVectorValue(ViterbiEndArray,
                   genomeSeqSize,
                   DOUBLE_NEGATIVE_INFINITY);
    VectorInteger  ViterbiEndState;
    setVectorValue(ViterbiEndState,
                   genomeSeqSize,
                   PROBALN_UNDEFINE);

    // Matrices and Vectors related to Forward algorithm
    Matrix3Double  ForwardMatrix;
    setMatrixValue(ForwardMatrix,
                   genomeSeqSize, readSeqSize, PROBALN_STATE_SPACE_SIZE,
                   0);
    VectorDouble   ForwardEndArray;
    setVectorValue(ForwardEndArray,
                   genomeSeqSize,
                   0);

    // Matrices and Vectors related to Feature
    Matrix4Double  TransitionFeatureMatrix;
    setMatrixValue(TransitionFeatureMatrix,
                   genomeSeqSize, readSeqSize, PROBALN_STATE_SPACE_SIZE, NumberTransitionFeatures,
                   0);
    Matrix2Double  TransitionFeatureEnd;
    setMatrixValue(TransitionFeatureEnd,
                   genomeSeqSize, NumberTransitionFeatures,
                   0);

    Matrix4Double  ReadHomopolymerFeatureMatrix;
    setMatrixValue(ReadHomopolymerFeatureMatrix,
                   genomeSeqSize, readSeqSize, PROBALN_STATE_SPACE_SIZE, NumberReadHomopolymerFeatures,
                   0);
    Matrix2Double  ReadHomopolymerFeatureEnd;
    setMatrixValue(ReadHomopolymerFeatureEnd,
                   genomeSeqSize, NumberReadHomopolymerFeatures,
                   0);

    Matrix4Double  GenomeHomopolymerFeatureMatrix;
    setMatrixValue(GenomeHomopolymerFeatureMatrix,
                   genomeSeqSize, readSeqSize, PROBALN_STATE_SPACE_SIZE, NumberGenomeHomopolymerFeatures,
                   0);
    Matrix2Double  GenomeHomopolymerFeatureEnd;
    setMatrixValue(GenomeHomopolymerFeatureEnd,
                   genomeSeqSize, NumberGenomeHomopolymerFeatures,
                   0);

    Matrix4Double  NeighborContextFeatureMatrix;
    setMatrixValue(NeighborContextFeatureMatrix,
                   genomeSeqSize, readSeqSize, PROBALN_STATE_SPACE_SIZE, NumberNeighborContextFeatures,
                   0);
    Matrix2Double  NeighborContextFeatureEnd;
    setMatrixValue(NeighborContextFeatureEnd,
                   genomeSeqSize, NumberNeighborContextFeatures,
                   0);

    Matrix4Double  EmissionFeatureMatrix;
    setMatrixValue(EmissionFeatureMatrix,
                   genomeSeqSize, readSeqSize, PROBALN_STATE_SPACE_SIZE, NumberEmissionFeatures,
                   0);
    Matrix2Double  EmissionFeatureEnd;
    setMatrixValue(EmissionFeatureEnd,
                   genomeSeqSize, NumberEmissionFeatures,
                   0);

    // sequence characteristics

    VectorInteger ReadPrecedeHomopolymerLength;
    VectorInteger ReadSuccedeHomopolymerLength;
    VectorInteger GenomePrecedeHomopolymerLength;
    VectorInteger GenomeSuccedeHomopolymerLength;
    VectorInteger LeftNeighborHomopolymerLength;
    VectorInteger RightNeighborHomopolymerLength;

    GenericFastaTools::markPrecedeHomopolymer(genomeSeq, GenomePrecedeHomopolymerLength);
    GenericFastaTools::markSuccedeHomopolymer(genomeSeq, GenomeSuccedeHomopolymerLength);
    GenericFastaTools::markPrecedeHomopolymer(readSeq, ReadPrecedeHomopolymerLength);
    GenericFastaTools::markSuccedeHomopolymer(readSeq, ReadSuccedeHomopolymerLength);
    GenericFastaTools::markLeftNeighborHomopolymer(readSeq, LeftNeighborHomopolymerLength);
    GenericFastaTools::markRightNeighborHomopolymer(readSeq, RightNeighborHomopolymerLength);

    // Viterbi iteration

    long double prevViterbiScore;
    long double prevForwardScore;

    long double maxViterbiScore;
    int         maxViterbiPrevState;
    int         maxViterbidi;
    int         maxViterbidj;

    long double wG;

    int ftIdx;
    ProbalnState prevState;
    ProbalnState currState;
    ProbalnTransitionFeature posFt;
    ProbalnReadHomopolymerFeature readHomoFt;
    ProbalnGenomeHomopolymerFeature genoHomoFt;
    ProbalnNeighborContextFeature neigContFt;
    ProbalnEmissionFeature emisFt;

    int left, right;

    for (int i=0; i<genomeSeqSize; ++i)
    {

        left  = leftBoundary[i];
        right = rightBoundary[i];

        for (int j=left; j<right; ++j)
        {

            //-------------------------------------------------
            // previous state is begin
            if (j==0)
            {
                prevState = PROBALN_BEGIN;
                prevViterbiScore = 0;
                prevForwardScore = 1;

                for (currState=PROBALN_MATCH; currState<=PROBALN_INSERT; ++currState)
                {

                    if (currState==PROBALN_MATCH)
                    {
                        if (readSeq[j]!=genomeSeq[i] && readSeq[j]!=Amb && genomeSeq[i]!=Amb)
                            continue;
                    }

                    if (currState==PROBALN_MISMATCH)
                    {
                        if (readSeq[j]==genomeSeq[i] || readSeq[j]==Amb || genomeSeq[i]==Amb)
                            continue;
                    }


                    getFeature(readSeq, genomeSeq, j, i, prevState, currState,
                               ReadSuccedeHomopolymerLength, GenomeSuccedeHomopolymerLength,
                               LeftNeighborHomopolymerLength, RightNeighborHomopolymerLength,
                               posFt, readHomoFt, genoHomoFt, neigContFt, emisFt);

                    wG = featureWeightSum(posFt, readHomoFt, genoHomoFt, neigContFt, emisFt);

                    // viterbi matrix
                    maxViterbiScore = wG;
                    maxViterbiPrevState = prevState;
                    ViterbiMatrix[i][j][currState] = maxViterbiScore;
                    ViterbiPreviousStateMatrix[i][j][currState] = maxViterbiPrevState;

                    // forward matrix
                    ForwardMatrix[i][j][currState] = exp(wG);

                    // feature expectation
                    if (&TransitionFeatureExpectations!=0)
                    {
                        // transition feature
                        if (posFt!=ProbalnTransitionFeature())
                        {
                            ftIdx = TransitionFeatureIndex[posFt];
                            TransitionFeatureMatrix[i][j][currState][ftIdx] = exp(wG);
                        }

                        // read homopolymer feature
                        if (readHomoFt!=ProbalnReadHomopolymerFeature())
                        {
                            ftIdx = ReadHomopolymerFeatureIndex[readHomoFt];
                            ReadHomopolymerFeatureMatrix[i][j][currState][ftIdx] = exp(wG);
                        }

                        // genome homopolymer feature
                        if (genoHomoFt!=ProbalnGenomeHomopolymerFeature())
                        {
                            ftIdx = GenomeHomopolymerFeatureIndex[genoHomoFt];
                            GenomeHomopolymerFeatureMatrix[i][j][currState][ftIdx] = exp(wG);
                        }

                        // neighbor context feature
                        if (neigContFt!=ProbalnNeighborContextFeature() && USE_MARKOV_FEATURE)
                        {
                            ftIdx = NeighborContextFeatureIndex[neigContFt];
                            NeighborContextFeatureMatrix[i][j][currState][ftIdx] = exp(wG);
                        }

                        // emission feature
                        if (emisFt!=ProbalnEmissionFeature())
                        {
                            ftIdx = EmissionFeatureIndex[emisFt];
                            EmissionFeatureMatrix[i][j][currState][ftIdx] = exp(wG);
                        }
                    }
                }
            }
            //----------------------------------------------
            // previous state is not begin
            else
            {
                for (currState=PROBALN_MATCH; currState<=PROBALN_DELETE; ++currState)
                {

                    if (currState==PROBALN_MATCH)
                    {
                        if (readSeq[j]!=genomeSeq[i] && readSeq[j]!=Amb && genomeSeq[i]!=Amb)
                            continue;
                    }

                    if (currState==PROBALN_MISMATCH)
                    {
                        if (readSeq[j]==genomeSeq[i] || readSeq[j]==Amb || genomeSeq[i]==Amb)
                            continue;
                    }

                    int di, dj;
                    if (currState==PROBALN_MATCH)
                    {
                        di = 1;
                        dj = 1;
                    }
                    if (currState==PROBALN_MISMATCH)
                    {
                        di = 1;
                        dj = 1;
                    }
                    if (currState==PROBALN_INSERT)
                    {
                        di = 0;
                        dj = 1;
                    }
                    if (currState==PROBALN_DELETE)
                    {
                        di = 1;
                        dj = 0;
                    }

                    if (i-di<0 || j-dj<0)
                        continue;

                    maxViterbiScore = DOUBLE_NEGATIVE_INFINITY;
                    maxViterbiPrevState = PROBALN_UNDEFINE;


                    for (prevState=PROBALN_MATCH; prevState<=PROBALN_DELETE; ++prevState)
                    {
                        getFeature(readSeq, genomeSeq, j, i, prevState, currState,
                                   ReadSuccedeHomopolymerLength, GenomeSuccedeHomopolymerLength,
                                   LeftNeighborHomopolymerLength, RightNeighborHomopolymerLength,
                                   posFt, readHomoFt, genoHomoFt, neigContFt, emisFt);

                        wG = featureWeightSum(posFt, readHomoFt, genoHomoFt, neigContFt, emisFt);

                        // viterbi computation

                        prevViterbiScore = ViterbiMatrix[i-di][j-dj][prevState];

                        if (maxViterbiScore<prevViterbiScore + wG)
                        {
                            maxViterbiScore = prevViterbiScore + wG;
                            maxViterbiPrevState = prevState;
                            maxViterbidi = di;
                            maxViterbidj = dj;
                        }

                        // forward computation

                        ForwardMatrix[i][j][currState] += exp(wG) * ForwardMatrix[i-di][j-dj][prevState];


                        // check inf or nan
                        if (ForwardMatrix[i][j][currState]!=ForwardMatrix[i][j][currState] ||
                                isinf(ForwardMatrix[i][j][currState]))
                        {
                            cout << "GenericProbabilisticAlignment ERROR: inf or nan value comes out... Aborting" << endl;
                            cout << "" << endl
                                 << "More information" << endl
                                 << "size(Read)=" << readSeqSize << endl
                                 << "size(Genome)=" << genomeSeqSize << endl
                                 << "state(Current)=" << ProbalnStateName[currState] << endl
                                 << "ptr(Read)=" << j << endl
                                 << "ptr(Genome)=" << i << endl
                                 << "wG=" << wG << ", " << "exp(wG)=" << exp(wG) << endl;

                            // transition feature
                            if (posFt!=ProbalnTransitionFeature())
                            {
                                cout << TransitionFeatureLabel[posFt] << " " << TransitionFeatureWeight[posFt] << endl;
                            }

                            // read homopolymer feature
                            if (readHomoFt!=ProbalnReadHomopolymerFeature())
                            {
                                cout << ReadHomopolymerFeatureLabel[readHomoFt] << " " << ReadHomopolymerFeatureWeight[readHomoFt] << endl;
                            }

                            // genome homopolymer feature
                            if (genoHomoFt!=ProbalnGenomeHomopolymerFeature())
                            {
                                cout << GenomeHomopolymerFeatureLabel[genoHomoFt] << " " << GenomeHomopolymerFeatureWeight[genoHomoFt] << endl;
                            }

                            // neighbor context feature
                            if (neigContFt!=ProbalnNeighborContextFeature() && USE_MARKOV_FEATURE)
                            {
                                cout << NeighborContextFeatureLabel[neigContFt] << " " << NeighborContextFeatureWeight[neigContFt] << endl;
                            }

                            // emission feature
                            if (emisFt!=ProbalnEmissionFeature())
                            {
                                cout << EmissionFeatureLabel[emisFt] << " " << EmissionFeatureWeight[emisFt] << endl;
                            }

                            cout << "dPtr(Read)=" << j-dj << endl
                                 << "dPtr(Genome)=" << i-di << endl
                                 << "state(Previous)=" << ProbalnStateName[prevState] << endl
                                 << "ForwardScore(Previous)=" << ForwardMatrix[i-di][j-dj][prevState] << endl;

                            cout << endl;

                            exit(EXIT_FAILURE);
                        }


                        // feature expectations

                        if (&TransitionFeatureExpectations!=0)
                        {
                            long double tempCount;

                            // transition feature
                            for (ftIdx=0; ftIdx<NumberTransitionFeatures; ++ftIdx)
                            {
                                tempCount = TransitionFeatureMatrix[i-di][j-dj][prevState][ftIdx];
                                if (posFt!=ProbalnTransitionFeature() && posFt==TransitionFeatureSet[ftIdx])
                                    tempCount += ForwardMatrix[i-di][j-dj][prevState];
                                TransitionFeatureMatrix[i][j][currState][ftIdx] += tempCount * exp(wG);
                            }

                            // read homopolymer feature
                            if (true)
                            {
                                for (ftIdx=0; ftIdx<NumberReadHomopolymerFeatures; ++ftIdx)
                                {
                                    tempCount =  ReadHomopolymerFeatureMatrix[i-di][j-dj][prevState][ftIdx];
                                    if (readHomoFt!=ProbalnReadHomopolymerFeature() && readHomoFt==ReadHomopolymerFeatureSet[ftIdx])
                                        tempCount += ForwardMatrix[i-di][j-dj][prevState];
                                    ReadHomopolymerFeatureMatrix[i][j][currState][ftIdx] += tempCount * exp(wG);
                                }
                            }

                            // genome homopolymer feature
                            if (true)
                            {
                                for (ftIdx=0; ftIdx<NumberGenomeHomopolymerFeatures; ++ftIdx)
                                {
                                    tempCount =  GenomeHomopolymerFeatureMatrix[i-di][j-dj][prevState][ftIdx];
                                    if (genoHomoFt!=ProbalnGenomeHomopolymerFeature() && genoHomoFt==GenomeHomopolymerFeatureSet[ftIdx])
                                        tempCount += ForwardMatrix[i-di][j-dj][prevState];
                                    GenomeHomopolymerFeatureMatrix[i][j][currState][ftIdx] += tempCount * exp(wG);
                                }
                            }

                            // neighbor context feature
                            if (USE_MARKOV_FEATURE)
                            {
                                for (ftIdx=0; ftIdx<NumberNeighborContextFeatures; ++ftIdx)
                                {
                                    tempCount =  NeighborContextFeatureMatrix[i-di][j-dj][prevState][ftIdx];
                                    if (neigContFt!=ProbalnNeighborContextFeature() && neigContFt==NeighborContextFeatureSet[ftIdx])
                                        tempCount += ForwardMatrix[i-di][j-dj][prevState];
                                    NeighborContextFeatureMatrix[i][j][currState][ftIdx] += tempCount * exp(wG);
                                }
                            }

                            // emission feature
                            for (ftIdx=0; ftIdx<NumberEmissionFeatures; ++ftIdx)
                            {
                                tempCount =  EmissionFeatureMatrix[i-di][j-dj][prevState][ftIdx];
                                if (emisFt!=ProbalnEmissionFeature() && emisFt==EmissionFeatureSet[ftIdx])
                                    tempCount += ForwardMatrix[i-di][j-dj][prevState];
                                EmissionFeatureMatrix[i][j][currState][ftIdx] += tempCount * exp(wG);
                            }
                        }
                    }

                    // save Viterbi results
                    ViterbiMatrix[i][j][currState] = maxViterbiScore;
                    ViterbiPreviousStateMatrix[i][j][currState] = maxViterbiPrevState;
                }
            }
        }
    }

    //-----------------------------------------------
    // current state is end

    for (int i=0; i<genomeSeqSize; ++i)
    {
        right = rightBoundary[i];
        if (right < readSeqSize)
            continue;

        // Viterbi computation
        currState = PROBALN_END;

        maxViterbiScore = DOUBLE_NEGATIVE_INFINITY;
        maxViterbiPrevState = PROBALN_UNDEFINE;

        for (prevState=PROBALN_MATCH; prevState<=PROBALN_DELETE; ++prevState)
        {
            posFt = ProbalnTransitionFeature(prevState, currState);
            wG    = transitionFeatureWeight(posFt, TransitionFeatureWeight);

            // viterbi
            prevViterbiScore = ViterbiMatrix[i][readSeqSize-1][prevState];
            if (maxViterbiScore<prevViterbiScore+wG)
            {
                maxViterbiScore = prevViterbiScore+wG;
                maxViterbiPrevState = prevState;
            }

            // forward
            ForwardEndArray[i] += ForwardMatrix[i][readSeqSize-1][prevState] * exp(wG);

            // feature
            if (&TransitionFeatureExpectations!=0)
            {
                long double tempCount;
                // transition feature
                for (ftIdx=0; ftIdx<NumberTransitionFeatures; ++ftIdx)
                {
                    tempCount = TransitionFeatureMatrix[i][readSeqSize-1][prevState][ftIdx];
                    if (posFt!=ProbalnTransitionFeature() && posFt==TransitionFeatureSet[ftIdx])
                        tempCount += ForwardMatrix[i][readSeqSize-1][prevState];
                    TransitionFeatureEnd[i][ftIdx] += tempCount * exp(wG);
                }

                // read homopolymer feature
                if (true)
                {
                    for (ftIdx=0; ftIdx<NumberReadHomopolymerFeatures; ++ftIdx)
                    {
                        tempCount = ReadHomopolymerFeatureMatrix[i][readSeqSize-1][prevState][ftIdx];
                        ReadHomopolymerFeatureEnd[i][ftIdx] += tempCount * exp(wG);
                    }
                }

                // genome homopolymer feature
                if (true)
                {
                    for (ftIdx=0; ftIdx<NumberGenomeHomopolymerFeatures; ++ftIdx)
                    {
                        tempCount = GenomeHomopolymerFeatureMatrix[i][readSeqSize-1][prevState][ftIdx];
                        GenomeHomopolymerFeatureEnd[i][ftIdx] += tempCount * exp(wG);
                    }
                }

                // neighbor context feature
                if (USE_MARKOV_FEATURE)
                {
                    for (ftIdx=0; ftIdx<NumberNeighborContextFeatures; ++ftIdx)
                    {
                        tempCount = NeighborContextFeatureMatrix[i][readSeqSize-1][prevState][ftIdx];
                        NeighborContextFeatureEnd[i][ftIdx] += tempCount * exp(wG);
                    }
                }

                // emission feature
                for (ftIdx=0; ftIdx<NumberEmissionFeatures; ++ftIdx)
                {
                    tempCount = EmissionFeatureMatrix[i][readSeqSize-1][prevState][ftIdx];
                    EmissionFeatureEnd[i][ftIdx] += tempCount * exp(wG);
                }
            }
        }

        // save viterbi result
        ViterbiEndArray[i] = maxViterbiScore;
        ViterbiEndState[i] = maxViterbiPrevState;
    }

    long double t_ViterbiScore, t_ForwardScore;
    // Viterbi Score
    t_ViterbiScore = DOUBLE_NEGATIVE_INFINITY;
    int t_ViterbiReadPtr   = readSeqSize-1;
    int t_ViterbiGenomePtr = genomeSeqSize-1;
    ProbalnState t_ViterbiState = PROBALN_UNDEFINE;

    for (int i=0; i<genomeSeqSize; ++i)
    {
        if (t_ViterbiScore<ViterbiEndArray[i])
        {
            t_ViterbiScore     = ViterbiEndArray[i];
            t_ViterbiGenomePtr = i;
            t_ViterbiState     = ViterbiEndState[i];
        }
    }

    // Forward Score
    t_ForwardScore = 0;
    for (VectorDouble::iterator iter=ForwardEndArray.begin(); iter!=ForwardEndArray.end(); ++iter)
    {
        t_ForwardScore += (*iter);
    }

    // alignment probability
    ProbabilisticScore = exp(t_ViterbiScore-log(t_ForwardScore));

    if (&ViterbiScore!=0)
        ViterbiScore = t_ViterbiScore;

    if (&ForwardScore!=0)
        ForwardScore = t_ForwardScore;

    // backtracking

    if (&alignReadSeq!=0 || &TransitionFeatureCounts!=0)
    {
        alignReadSeq   = string();
        alignGenomeSeq = string();

        while(t_ViterbiState!=PROBALN_BEGIN && t_ViterbiState!=PROBALN_UNDEFINE)
        {
            stringstream alignReadBuf;
            stringstream alignGenomeBuf;

            int di,dj;

            if (t_ViterbiState==PROBALN_MATCH)
            {
                di=1;
                dj=1;

                alignReadBuf   << readSeq.substr(t_ViterbiReadPtr,dj) << alignReadSeq;
                alignGenomeBuf << genomeSeq.substr(t_ViterbiGenomePtr,di) << alignGenomeSeq;
            }

            if (t_ViterbiState==PROBALN_MISMATCH)
            {
                di=1;
                dj=1;

                alignReadBuf   << readSeq.substr(t_ViterbiReadPtr,dj) << alignReadSeq;
                alignGenomeBuf << genomeSeq.substr(t_ViterbiGenomePtr,di) << alignGenomeSeq;
            }

            if (t_ViterbiState==PROBALN_INSERT)
            {
                di = 0;
                dj = 1;

                alignReadBuf   << readSeq.substr(t_ViterbiReadPtr,dj) << alignReadSeq;
                alignGenomeBuf << "-" << alignGenomeSeq;
            }

            if (t_ViterbiState==PROBALN_DELETE)
            {
                di = 1;
                dj = 0;

                alignReadBuf   << "-" << alignReadSeq;
                alignGenomeBuf << genomeSeq.substr(t_ViterbiGenomePtr,di) << alignGenomeSeq;
            }

            // update
            if (t_ViterbiReadPtr>=0)
                t_ViterbiState = ViterbiPreviousStateMatrix[t_ViterbiGenomePtr][t_ViterbiReadPtr][t_ViterbiState];
            else
                t_ViterbiState = PROBALN_BEGIN;

            t_ViterbiReadPtr   -= dj;
            t_ViterbiGenomePtr -= di;

            alignReadSeq     = alignReadBuf.str();
            alignGenomeSeq   = alignGenomeBuf.str();
        }

        // move indel to right
        GenericBamAlignmentTools::moveInDelRight(alignReadSeq, alignGenomeSeq);
    }


    // empirical counting
    if (&TransitionFeatureCounts!=0)
    {

        // debug
        if (alignReadSeq.empty() || alignGenomeSeq.empty())
        {
            cout << endl;
            cout << ">" << readSeqSize << " "
                 << genomeSeqSize << " "
                 << ProbabilisticScore << " "
                 << ViterbiScore << " "
                 << ForwardScore << " "
                 << genomeSeqSize << " "
                 << readSeqSize
                 << endl;
            cout << readSeq << endl;
            cout << genomeSeq << endl;
            cout << alignReadSeq << endl;
            cout << alignGenomeSeq << endl;
            cout << endl;
        }

        getAlignmentFeatureCount(alignReadSeq, alignGenomeSeq,
                                 TransitionFeatureCounts,
                                 ReadHomopolymerFeatureCounts,
                                 GenomeHomopolymerFeatureCounts,
                                 NeighborContextFeatureCounts,
                                 EmissionFeatureCounts);

        // debug transition feature
        if (0)
        {
            printTransitionFeatureCount(TransitionFeatureCounts);
        }
    }

    // expected counting
    if (&TransitionFeatureExpectations!=0)
    {
        TransitionFeatureExpectations.assign(NumberTransitionFeatures, 0);
        ReadHomopolymerFeatureExpectations.assign(NumberReadHomopolymerFeatures, 0);
        GenomeHomopolymerFeatureExpectations.assign(NumberGenomeHomopolymerFeatures, 0);
        NeighborContextFeatureExpectations.assign(NumberNeighborContextFeatures, 0);
        EmissionFeatureExpectations.assign(NumberEmissionFeatures, 0);

        for (int i=0; i<genomeSeqSize; ++i)
        {
            // transition feature
            for (int idx=0; idx<NumberTransitionFeatures; ++idx)
            {
                TransitionFeatureExpectations[idx] += TransitionFeatureEnd[i][idx]/ForwardScore;
            }

            // read homopolymer feature
            if (true)
            {
                for (int idx=0; idx<NumberReadHomopolymerFeatures; ++idx)
                {
                    ReadHomopolymerFeatureExpectations[idx] += ReadHomopolymerFeatureEnd[i][idx]/ForwardScore;
                }
            }

            // genome homopolymer feature
            if (true)
            {
                for (int idx=0; idx<NumberGenomeHomopolymerFeatures; ++idx)
                {
                    GenomeHomopolymerFeatureExpectations[idx] += GenomeHomopolymerFeatureEnd[i][idx]/ForwardScore;
                }
            }

            // neighbor context feature
            if (USE_MARKOV_FEATURE)
            {
                for (int idx=0; idx<NumberNeighborContextFeatures; ++idx)
                {
                    NeighborContextFeatureExpectations[idx] += NeighborContextFeatureEnd[i][idx]/ForwardScore;
                }
            }

            // emission feature
            for (int idx=0; idx<NumberEmissionFeatures; ++idx)
            {
                EmissionFeatureExpectations[idx] += EmissionFeatureEnd[i][idx]/ForwardScore;
            }
        }
    }

    // return alignment probability
    return ProbabilisticScore;
}


void GenericProbabilisticAlignment::getAlignmentFeatureCount(string& alignReadSeq, string& alignGenomeSeq,
                                                             vector<long double> &transitionFeatureCount,
                                                             vector<long double> &readHomopolymerFeatureCount,
                                                             vector<long double> &genomeHomopolymerFeatureCount,
                                                             vector<long double> &neighborContextFeatureCount,
                                                             vector<long double> &emissionFeatureCount)
{
    int readPointer, genomePointer;
    int readLength;
    int i;
    int idx;

    // read sequence and genome sequence
    string readSeq(""), genomeSeq("");
    for (int i=0; i<alignReadSeq.length(); ++i)
    {
        if (alignReadSeq[i]!=Spa)
            readSeq += alignReadSeq[i];
    }
    for (int i=0; i<alignGenomeSeq.length(); ++i)
    {
        if (alignGenomeSeq[i]!=Spa)
            genomeSeq += alignGenomeSeq[i];
    }

    vector<int> readPrecedeHomopolymerLen;
    vector<int> readPrecedeHomopolymerLen2;
    vector<int> genomePrecedeHomopolymerLen;
    vector<int> genomePrecedeHomopolymerLen2;
    vector<int> leftMarkovianHomopolymerLen;
    vector<int> rightMarkovianHOmopolymerLen;

    GenericFastaTools::markPrecedeHomopolymer(alignReadSeq, readPrecedeHomopolymerLen);
    GenericFastaTools::markPrecedeHomopolymer(alignGenomeSeq, genomePrecedeHomopolymerLen);
    GenericFastaTools::markLeftNeighborHomopolymer(alignReadSeq, leftMarkovianHomopolymerLen);
    GenericFastaTools::markRightNeighborHomopolymer(alignReadSeq, rightMarkovianHOmopolymerLen);

    GenericFastaTools::markPrecedeHomopolymer(readSeq, readPrecedeHomopolymerLen2);
    GenericFastaTools::markPrecedeHomopolymer(genomeSeq, genomePrecedeHomopolymerLen2);

    // initialize zero
    transitionFeatureCount.assign(NumberTransitionFeatures, 0);
    readHomopolymerFeatureCount.assign(NumberReadHomopolymerFeatures, 0);
    genomeHomopolymerFeatureCount.assign(NumberGenomeHomopolymerFeatures, 0);
    neighborContextFeatureCount.assign(NumberNeighborContextFeatures, 0);
    emissionFeatureCount.assign(NumberEmissionFeatures, 0);

    // alignment state sequence
    vector<ProbalnState> alignState;
    for (int i=0; i<alignReadSeq.length(); ++i)
    {
        if (alignReadSeq[i]==Spa)
        {
            alignState.push_back(PROBALN_DELETE);
            continue;
        }
        if (alignGenomeSeq[i]==Spa)
        {
            alignState.push_back(PROBALN_INSERT);
            continue;
        }

        if (alignReadSeq[i]!=alignGenomeSeq[i] && alignReadSeq[i]!=Amb && alignGenomeSeq[i]!=Amb)
        {
            alignState.push_back(PROBALN_MISMATCH);
            continue;
        }

        alignState.push_back(PROBALN_MATCH);
    }

    // read length
    readLength = 0;
    for (string::iterator iter=alignReadSeq.begin(); iter!=alignReadSeq.end(); ++iter)
    {
        if ((*iter)!=Spa)
            readLength++;
    }

    // iterate from first position to last position in alignment
    ProbalnTransitionFeature        posFt;
    ProbalnReadHomopolymerFeature   rhFt;
    ProbalnGenomeHomopolymerFeature ghFt;
    ProbalnNeighborContextFeature   lcFt;
    ProbalnEmissionFeature          eFt;

    readPointer   = -1;
    genomePointer = -1;
    for (i=0; i<alignReadSeq.length(); ++i)
    {

        // update read pointer
        if (alignReadSeq[i]!=Spa)
            readPointer++;

        // update genome pointer
        if (alignGenomeSeq[i]!=Spa)
            genomePointer++;


        // transition feature
        if (i==0)
        {
            posFt = ProbalnTransitionFeature(PROBALN_BEGIN, alignState[i]);
            idx = TransitionFeatureIndex[posFt];
            transitionFeatureCount[idx] ++;
        }else
        {
            posFt = ProbalnTransitionFeature(alignState[i-1], alignState[i]);
            idx = TransitionFeatureIndex[posFt];
            transitionFeatureCount[idx] ++;
        }

        if (USE_MARKOV_FEATURE)
        {
            // neighbor context feature
            if (alignState[i]==PROBALN_INSERT &&
                    readSeq[readPointer]!=genomeSeq[genomePointer] &&
                    (readSeq[readPointer]!=genomeSeq[genomePointer+1] && genomePointer+1<genomeSeq.length()))
            {
                int len = leftMarkovianHomopolymerLen[i];
                if (len<rightMarkovianHOmopolymerLen[i])
                    len = rightMarkovianHOmopolymerLen[i];

                if (len>=midHomopolymerSize && len<=maxHomopolymerSize)
                {
                    lcFt = ProbalnNeighborContextFeature(alignState[i], len);
                    idx = NeighborContextFeatureIndex[lcFt];
                    neighborContextFeatureCount[idx] ++;
                }
            }
        }

        // emission feature
        if (alignState[i]==PROBALN_MATCH)
        {
            eFt = ProbalnEmissionFeature(alignState[i], Alpha2Base[alignGenomeSeq[i]], Alpha2Base[alignReadSeq[i]]);
            if (alignGenomeSeq[i]==Amb)
                eFt = ProbalnEmissionFeature(alignState[i], Alpha2Base[alignReadSeq[i]], Alpha2Base[alignReadSeq[i]]);
            if (alignReadSeq[i]==Amb)
                eFt = ProbalnEmissionFeature(alignState[i], Alpha2Base[alignGenomeSeq[i]], Alpha2Base[alignGenomeSeq[i]]);

            idx = EmissionFeatureIndex[eFt];
            emissionFeatureCount[idx] ++;
        }
        if (alignState[i]==PROBALN_MISMATCH)
        {
            eFt = ProbalnEmissionFeature(alignState[i], Alpha2Base[alignGenomeSeq[i]], Alpha2Base[alignReadSeq[i]]);
            idx = EmissionFeatureIndex[eFt];
            emissionFeatureCount[idx] ++;
        }
        if (alignState[i]==PROBALN_INSERT)
        {
            if (alignReadSeq[i]!=Amb)
                eFt = ProbalnEmissionFeature(alignState[i], SPACE, Alpha2Base[alignReadSeq[i]]);
            else
                eFt = ProbalnEmissionFeature(alignState[i], SPACE, ADENINE);

            idx = EmissionFeatureIndex[eFt];
            emissionFeatureCount[idx] ++;
        }
        if (alignState[i]==PROBALN_DELETE)
        {
            if (alignGenomeSeq[i]!=Amb)
                eFt = ProbalnEmissionFeature(alignState[i], Alpha2Base[alignGenomeSeq[i]], SPACE);
            else
                eFt = ProbalnEmissionFeature(alignState[i], ADENINE, SPACE);

            idx = EmissionFeatureIndex[eFt];
            emissionFeatureCount[idx] ++;
        }

    }

    posFt = ProbalnTransitionFeature(alignState[alignState.size()-1], PROBALN_END);
    idx = TransitionFeatureIndex[posFt];
    transitionFeatureCount[idx] ++;

    // divide alignment into homopolymer-wise block
    if (true)
    {
        vector<BamAlignmentBlock> blocks;
        GenericBamAlignmentTools::divideAlignmentToBlocks(alignReadSeq, alignGenomeSeq, blocks);
        for (int i=0; i<blocks.size(); ++i)
        {
            int rhl, thl, delta;
            // insert
            if (blocks[i].m_read.length()>blocks[i].m_genome.length() && blocks[i].m_genome.length()>0)
            {
                thl   = blocks[i].m_genome.length();

                for (int len=thl+1; len<=blocks[i].m_read.length(); ++len)
                {

                    rhl = len;
                    delta = rhl-thl;

                    if (rhl>maxHomopolymerSize)
                        rhl = maxHomopolymerSize;
                    if (delta>5)
                        delta = 5;

                    if (rhl>=midHomopolymerSize && rhl<=maxHomopolymerSize && delta>=1 && delta<=5)
                    {
                        rhFt = ProbalnReadHomopolymerFeature(PROBALN_INSERT, rhl, delta);
                        idx = ReadHomopolymerFeatureIndex[rhFt];
                        readHomopolymerFeatureCount[idx] ++;
                    }
                }
            }

            // delete
            if (blocks[i].m_genome.length()>blocks[i].m_read.length() && blocks[i].m_read.length()>0)
            {
                rhl   = blocks[i].m_read.length();

                for (int len=rhl+1; len<=blocks[i].m_genome.length(); ++len)
                {
                    thl   = len;
                    delta = thl-rhl;

                    if (thl>maxHomopolymerSize)
                        thl=maxHomopolymerSize;
                    if (delta>5)
                        delta=5;

                    if (thl>=midHomopolymerSize && thl<=maxHomopolymerSize && delta>=1 && delta<=5)
                    {
                        ghFt = ProbalnGenomeHomopolymerFeature(PROBALN_DELETE, thl, delta);
                        idx = GenomeHomopolymerFeatureIndex[ghFt];
                        genomeHomopolymerFeatureCount[idx] ++;
                    }
                }
            }
        }
    }
}


void GenericProbabilisticAlignment::initModelSetting(vector<string> &alignReadSeqs, vector<string> &alignGenomeSeqs)
{
    VectorDouble transitionFeatureCounts(NumberTransitionFeatures, 0);
    VectorDouble readHomopolymerFeatureCounts(NumberReadHomopolymerFeatures, 0);
    VectorDouble genomeHomopolymerFeatureCounts(NumberGenomeHomopolymerFeatures, 0);
    VectorDouble neighborContextFeatureCounts(NumberNeighborContextFeatures, 0);
    VectorDouble emissionFeatureCounts(NumberEmissionFeatures, 0);

    VectorDouble transitionFeatureCount;
    VectorDouble readHomopolymerFeatureCount;
    VectorDouble genomeHomopolymerFeatureCount;
    VectorDouble neighborContextFeatureCount;
    VectorDouble emissionFeatureCount;

    VectorDouble ghlens(maxHomopolymerSize+1, 0);
    VectorDouble rhlens(maxHomopolymerSize+1, 0);
    VectorDouble lchlens(maxHomopolymerSize+1, 0);
    VectorDouble rchlens(maxHomopolymerSize+1, 0);

    for (int i=0; i<alignReadSeqs.size(); i++)
    {
        string alnReadSeq   = alignReadSeqs[i];
        string alnGenomeSeq = alignGenomeSeqs[i];

        // move indel to right
        GenericBamAlignmentTools::moveInDelRight(alnReadSeq, alnGenomeSeq);

        // get feature count
        getAlignmentFeatureCount(alnReadSeq, alnGenomeSeq,
                                 transitionFeatureCount,
                                 readHomopolymerFeatureCount,
                                 genomeHomopolymerFeatureCount,
                                 neighborContextFeatureCount,
                                 emissionFeatureCount);

        // update transition feature
        for (int j=0; j<NumberTransitionFeatures; ++j)
            transitionFeatureCounts[j] += transitionFeatureCount[j];

        // update read homopolymer feature
        if (true)
        {
            for (int j=0; j<NumberReadHomopolymerFeatures; ++j)
                readHomopolymerFeatureCounts[j] += readHomopolymerFeatureCount[j];
        }

        // update genome homopolymer feature
        if (true)
        {
            for (int j=0; j<NumberGenomeHomopolymerFeatures; ++j)
                genomeHomopolymerFeatureCounts[j] += genomeHomopolymerFeatureCount[j];
        }

        if (USE_MARKOV_FEATURE)
        {
            // update neighbor context feature
            for (int j=0; j<NumberNeighborContextFeatures; ++j)
                neighborContextFeatureCounts[j] += neighborContextFeatureCount[j];
        }

        // update emission feature
        for (int j=0; j<NumberEmissionFeatures; j++)
            emissionFeatureCounts[j] += emissionFeatureCount[j];

        // homopolymer length
        string t_read("");
        string t_genome("");
        for (string::iterator iter=alnReadSeq.begin(); iter!=alnReadSeq.end(); ++iter)
        {
            if ((*iter)!=Spa)
                t_read += (*iter);
        }
        for (string::iterator iter=alnGenomeSeq.begin(); iter!=alnGenomeSeq.end(); ++iter)
        {
            if ((*iter)!=Spa)
                t_genome += (*iter);
        }

        // homopolymer feature
        vector<int> t_ghlen;
        GenericFastaTools::markPrecedeHomopolymer(t_genome, t_ghlen);

        vector<int> t_rhlen;
        GenericFastaTools::markPrecedeHomopolymer(t_read, t_rhlen);

        vector<int> t_lchlen;
        GenericFastaTools::markLeftNeighborHomopolymer(t_read, t_lchlen);

        vector<int> t_rchlen;
        GenericFastaTools::markRightNeighborHomopolymer(t_read, t_rchlen);

        for (vector<int>::iterator iter=t_ghlen.begin(); iter!=t_ghlen.end(); ++iter)
        {
            if ((*iter)>=midHomopolymerSize && (*iter)<=maxHomopolymerSize)
                ghlens[*iter] ++;
        }

        for (vector<int>::iterator iter=t_rhlen.begin(); iter!=t_rhlen.end(); ++iter)
        {
            if ((*iter)>=midHomopolymerSize && (*iter)<=maxHomopolymerSize)
                rhlens[*iter] ++;
        }

        if (USE_MARKOV_FEATURE)
        {

            // neighbor context feature
            for (int k=0; k<t_lchlen.size(); ++k)
            {
                int len = t_lchlen[k];
                if (len<t_rchlen[k])
                    len = t_rchlen[k];

                if (len>=midHomopolymerSize && len<=maxHomopolymerSize)
                    lchlens[len] ++;
            }
        }
    }

    // model setting initialization

    // transition feature
    initTransitionFeatureWeight(transitionFeatureCounts);

    // read homopolymer feature
    if (true)
        initReadHomopolymerFeatureWeight(readHomopolymerFeatureCounts, rhlens);

    // genome homopolymer feature
    if (true)
        initGenomeHomopolymerFeatureWeight(genomeHomopolymerFeatureCounts, ghlens);

    if (USE_MARKOV_FEATURE)
    {
        // neighbor context feature
        initNeighborContextFeatureWeight(neighborContextFeatureCounts, lchlens);
    }

    // emission feature
    initEmissionFeatureWeight(emissionFeatureCounts);
}

// transition feature
void GenericProbabilisticAlignment::initTransitionFeatureWeight(VectorDouble &counts)
{
    ProbalnTransitionFeature f;
    ProbalnState prev;
    ProbalnState curr;
    double z;
    int idx;

    // add pseudo counts
    for (VectorDouble::iterator iter=counts.begin(); iter!=counts.end(); ++iter)
        (*iter) += 1;

    // begin
    z = 0;
    prev = PROBALN_BEGIN;
    for (curr=PROBALN_MATCH; curr<=PROBALN_INSERT; ++curr)
    {
        f = ProbalnTransitionFeature(prev, curr);
        idx = TransitionFeatureIndex[f];
        z += counts[idx];
    }
    for (curr=PROBALN_MATCH; curr<=PROBALN_INSERT; ++curr)
    {
        f = ProbalnTransitionFeature(prev, curr);
        idx = TransitionFeatureIndex[f];
        if (curr==PROBALN_MATCH)
        {
            TransitionFeatureWeight[f] = 0;
        }else
        {
            TransitionFeatureWeight[f] = log(counts[idx]) - log(z);
        }
    }

    // end
    z = 0;
    curr = PROBALN_END;
    for (prev=PROBALN_MATCH; prev<=PROBALN_DELETE; ++prev)
    {
        f = ProbalnTransitionFeature(prev, curr);
        idx = TransitionFeatureIndex[f];
        z += counts[idx];
    }
    for (prev=PROBALN_MATCH; prev<=PROBALN_DELETE; ++prev)
    {
        f = ProbalnTransitionFeature(prev, curr);
        idx = TransitionFeatureIndex[f];
        TransitionFeatureWeight[f] = log(counts[idx]) - log(z);
    }

    // middle
    for (prev=PROBALN_MATCH; prev<=PROBALN_DELETE; ++prev)
    {
        z = 0;
        for (curr=PROBALN_MATCH; curr<=PROBALN_DELETE; ++curr)
        {
            f = ProbalnTransitionFeature(prev, curr);
            idx = TransitionFeatureIndex[f];
            z += counts[idx];
        }

        for (curr=PROBALN_MATCH; curr<=PROBALN_DELETE; ++curr)
        {

            f = ProbalnTransitionFeature(prev, curr);
            idx = TransitionFeatureIndex[f];
            if (curr==PROBALN_MATCH)
            {
                TransitionFeatureWeight[f] = 0;
            }else
            {
                TransitionFeatureWeight[f] = log(counts[idx]) - log(z);

                if (prev==PROBALN_DELETE && curr==PROBALN_DELETE)
                    TransitionFeatureWeight[f] += 0.5;
                if (prev==PROBALN_INSERT && curr==PROBALN_INSERT)
                    TransitionFeatureWeight[f] += 0.5;
                if (prev==PROBALN_DELETE && curr==PROBALN_INSERT)
                    TransitionFeatureWeight[f] -= 0.5;
                if (prev==PROBALN_INSERT && curr==PROBALN_DELETE)
                    TransitionFeatureWeight[f] -= 0.5;
            }
        }
    }
}

// read homopolymer feature
void GenericProbabilisticAlignment::initReadHomopolymerFeatureWeight(VectorDouble &counts, VectorDouble &hlen)
{
    ProbalnState curr = PROBALN_INSERT;
    ProbalnReadHomopolymerFeature f;
    int i;

    for (int l=midHomopolymerSize; l<=maxHomopolymerSize; l++)
    {
        for (int delta=1; delta<=l; delta++)
        {
            if (delta>5)
                continue;

            f = ProbalnReadHomopolymerFeature(curr, l, delta);
            i = ReadHomopolymerFeatureIndex[f];
            ReadHomopolymerFeatureWeight[f] = 0.0;
        }
    }
}

// genome homopolymer feature
void GenericProbabilisticAlignment::initGenomeHomopolymerFeatureWeight(VectorDouble &counts, VectorDouble &hlen)
{
    ProbalnState curr = PROBALN_DELETE;
    ProbalnGenomeHomopolymerFeature f;
    int i;

    for (int l=midHomopolymerSize; l<=maxHomopolymerSize; l++)
    {
        for (int delta=1; delta<=l; delta++)
        {
            if (delta>5)
                continue;

            f = ProbalnGenomeHomopolymerFeature(curr, l, delta);
            i = GenomeHomopolymerFeatureIndex[f];
            GenomeHomopolymerFeatureWeight[f] = 0.0;
        }
    }
}

void GenericProbabilisticAlignment::initNeighborContextFeatureWeight(VectorDouble &counts, VectorDouble &hlen)
{
    ProbalnState curr = PROBALN_INSERT;
    ProbalnNeighborContextFeature f;
    int i;

    for (int l=midHomopolymerSize; l<=maxHomopolymerSize; l++)
    {
        f = ProbalnNeighborContextFeature(curr, l);
        i = NeighborContextFeatureIndex[f];
        NeighborContextFeatureWeight[f] = 0.005*l;
    }
}


// emission feature
void GenericProbabilisticAlignment::initEmissionFeatureWeight(VectorDouble &counts)
{
    ProbalnState curr;
    ProbalnEmissionFeature f;
    SequenceBase a, b;
    long double z;
    int idx;

    // add pseudo counts
    for (VectorDouble::iterator iter=counts.begin(); iter!=counts.end(); ++iter)
        (*iter) += 1;

    // total count
    z = 0;
    for (curr=PROBALN_MATCH; curr<=PROBALN_DELETE; ++curr)
    {
        for (a=ADENINE; a<=THYMINE; ++a)
        {
            for (b=ADENINE; b<=THYMINE; ++b)
            {
                f = ProbalnEmissionFeature(curr, a, b);
                if (EmissionFeatureWeight.find(f)!=EmissionFeatureWeight.end())
                {
                    idx = EmissionFeatureIndex[f];
                    z   += counts[idx];
                }

            }
        }
    }

    // match
    curr = PROBALN_MATCH;
    for (a=ADENINE; a<=THYMINE; ++a)
    {
        f = ProbalnEmissionFeature(curr, a, a);
        EmissionFeatureWeight[f] = 1;
    }

    // mismatch
    curr = PROBALN_MISMATCH;
    for (a=ADENINE; a<=THYMINE; ++a)
    {
        for (b=ADENINE; b<=THYMINE; ++b)
        {
            if (a==b) continue;

            f = ProbalnEmissionFeature(curr, a, b);
            idx = EmissionFeatureIndex[f];

            EmissionFeatureWeight[f] = log(counts[idx]) -log(z);
        }
    }

    // insert
    curr = PROBALN_INSERT;
    a = SPACE;
    for (b=ADENINE; b<=THYMINE; ++b)
    {
        f = ProbalnEmissionFeature(curr, a, b);
        idx = EmissionFeatureIndex[f];

        EmissionFeatureWeight[f] = log(counts[idx]) -log(z) + .5;
    }

    // delete
    curr = PROBALN_DELETE;
    b = SPACE;
    for (a=ADENINE; a<=THYMINE; ++a)
    {
        f = ProbalnEmissionFeature(curr, a, b);
        idx = EmissionFeatureIndex[f];

        EmissionFeatureWeight[f] = log(counts[idx]) -log(z) + .5;
    }
}

void GenericProbabilisticAlignment::print(ostream &output)
{
    // transition feature
    for (int i=0; i<NumberTransitionFeatures; ++i)
    {
        ProbalnTransitionFeature f = TransitionFeatureSet[i];
        output << TransitionFeatureLabel[f] << "\t"
               << TransitionFeatureWeight[f] << endl;
    }

    // read homopolymer feature
    for (int i=0; i<NumberReadHomopolymerFeatures; ++i)
    {
        ProbalnReadHomopolymerFeature f = ReadHomopolymerFeatureSet[i];
        output << ReadHomopolymerFeatureLabel[f] << "\t"
               << ReadHomopolymerFeatureWeight[f] << endl;
    }

    // genome homopolymer feature
    for (int i=0; i<NumberGenomeHomopolymerFeatures; ++i)
    {
        ProbalnGenomeHomopolymerFeature f = GenomeHomopolymerFeatureSet[i];
        output << GenomeHomopolymerFeatureLabel[f] << "\t"
               << GenomeHomopolymerFeatureWeight[f] << endl;
    }

    // neighbor context feature
    for (int i=0; i<NumberNeighborContextFeatures; ++i)
    {
        ProbalnNeighborContextFeature f = NeighborContextFeatureSet[i];
        output << NeighborContextFeatureLabel[f] << "\t"
               << NeighborContextFeatureWeight[f] << endl;
    }

    // emission feature
    for (int i=0; i<NumberEmissionFeatures; ++i)
    {
        ProbalnEmissionFeature f = EmissionFeatureSet[i];
        output << EmissionFeatureLabel[f] << "\t"
               << EmissionFeatureWeight[f] << endl;
    }
}


void GenericProbabilisticAlignment::train(vector<string> &alignReadSeqs, vector<string> &alignGenomeSeqs, TrainSetting& settings)
{
    // parameter setting initialization
    initModelSetting(alignReadSeqs, alignGenomeSeqs);

    // erase SPACE from sequences
    vector<string> readSeqs;

    for (vector<string>::iterator iter=alignReadSeqs.begin(); iter!=alignReadSeqs.end(); ++iter)
    {
        string seq("");
        for (string::iterator it=iter->begin(); it!=iter->end(); ++it)
        {
            if ((*it)!=Spa)
                seq += (*it);
        }
        readSeqs.push_back(seq);
    }

    vector<string> genomeSeqs;

    for (vector<string>::iterator iter=alignGenomeSeqs.begin(); iter!=alignGenomeSeqs.end(); ++iter)
    {
        string seq("");
        for (string::iterator it=iter->begin(); it!=iter->end(); ++it)
        {
            if ((*it)!=Spa)
                seq += (*it);
        }
        genomeSeqs.push_back(seq);
    }

    // parameter estimation
    nlopt_aux_t aux_data;
    aux_data.readSeqs     = &readSeqs;
    aux_data.genomeSeqs   = &genomeSeqs;
    aux_data.band         = settings.m_band;
    aux_data.numTrainData = readSeqs.size();
    aux_data.iter         = 0;
    aux_data.verbosity    = settings.m_verbosity;
    aux_data.ProbAligner  = this;
    aux_data.optim_loglik = DOUBLE_NEGATIVE_INFINITY;


    int dim = NumberTransitionFeatures + NumberEmissionFeatures;
    if (true)
        dim += NumberReadHomopolymerFeatures+NumberGenomeHomopolymerFeatures;
    if (USE_MARKOV_FEATURE)
        dim += NumberNeighborContextFeatures;


    nlopt::opt iterative_trainer(((USE_NLOPT_LBFGS) ? nlopt::LD_LBFGS : nlopt::LD_SLSQP),
                                 dim);


    iterative_trainer.set_max_objective(loglik, &aux_data);
    iterative_trainer.set_xtol_abs(settings.m_stopVal);
    iterative_trainer.set_ftol_abs(settings.m_stopVal);
    iterative_trainer.set_ftol_rel(settings.m_stopVal);
    iterative_trainer.set_xtol_rel(settings.m_stopVal);
    iterative_trainer.set_maxeval(settings.m_maxIter);

    // initial weight
    vector<double> w;

    // transition feature
    for (int i=0; i<NumberTransitionFeatures; ++i)
    {
        ProbalnTransitionFeature f = TransitionFeatureSet[i];
        w.push_back(TransitionFeatureWeight[f]);
    }

    // read homopolymer feature
    if (true)
    {
        for (int i=0; i<NumberReadHomopolymerFeatures; ++i)
        {
            ProbalnReadHomopolymerFeature f = ReadHomopolymerFeatureSet[i];
            w.push_back(ReadHomopolymerFeatureWeight[f]);
        }
    }

    // genome homopolymer feature
    if (true)
    {
        for (int i=0; i<NumberGenomeHomopolymerFeatures; ++i)
        {
            ProbalnGenomeHomopolymerFeature f = GenomeHomopolymerFeatureSet[i];
            w.push_back(GenomeHomopolymerFeatureWeight[f]);
        }
    }

    // neighbor context feature
    if (USE_MARKOV_FEATURE)
    {
        for (int i=0; i<NumberNeighborContextFeatures; ++i)
        {
            ProbalnNeighborContextFeature f = NeighborContextFeatureSet[i];
            w.push_back(NeighborContextFeatureWeight[f]);
        }
    }

    // emission feature
    for (int i=0; i<NumberEmissionFeatures; ++i)
    {
        ProbalnEmissionFeature f = EmissionFeatureSet[i];
        w.push_back(EmissionFeatureWeight[f]);
    }

    // train
    double score;
    try
    {
        iterative_trainer.optimize(w, score);
    }
    catch(nlopt::roundoff_limited&)
    {
        w.assign(aux_data.optim_w.begin(), aux_data.optim_w.end());
        iterative_trainer.force_stop();
    }
    catch(nlopt::forced_stop&)
    {
        w.assign(aux_data.optim_w.begin(), aux_data.optim_w.end());
        iterative_trainer.force_stop();
    }
    catch(std::runtime_error&)
    {
        w.assign(aux_data.optim_w.begin(), aux_data.optim_w.end());
        iterative_trainer.force_stop();
    }

    // save trained parameters
    int shift,n=0;

    // position feature
    shift = n;
    for (int i=0; i<NumberTransitionFeatures; ++i, ++n)
    {
        ProbalnTransitionFeature f = TransitionFeatureSet[i];
        TransitionFeatureWeight[f] = w[i+shift];
    }

    // read homopolymer feature
    if (true)
    {
        shift = n;
        for (int i=0; i<NumberReadHomopolymerFeatures; ++i, ++n)
        {
            ProbalnReadHomopolymerFeature f = ReadHomopolymerFeatureSet[i];
            ReadHomopolymerFeatureWeight[f] = w[i+shift];
        }
    }

    // genome homopolymer feature
    if (true)
    {
        shift = n;
        for (int i=0; i<NumberGenomeHomopolymerFeatures; ++i, ++n)
        {
            ProbalnGenomeHomopolymerFeature f = GenomeHomopolymerFeatureSet[i];
            GenomeHomopolymerFeatureWeight[f] = w[i+shift];
        }
    }

    // neighbor context feature
    if (USE_MARKOV_FEATURE)
    {
        shift = n;
        for (int i=0; i<NumberNeighborContextFeatures; ++i, ++n)
        {
            ProbalnNeighborContextFeature f = NeighborContextFeatureSet[i];
            NeighborContextFeatureWeight[f] = w[i+shift];
        }
    }

    // emission feature
    shift = n;
    for (int i=0; i<NumberEmissionFeatures; ++i, ++n)
    {
        ProbalnEmissionFeature f = EmissionFeatureSet[i];
        EmissionFeatureWeight[f] = w[i+shift];
    }

}


void GenericProbabilisticAlignment::printTransitionFeatureCount(VectorDouble &count)
{
    for (int i=0; i<NumberTransitionFeatures; ++i)
    {
        ProbalnTransitionFeature f = TransitionFeatureSet[i];
        cout << TransitionFeatureLabel[f] << " "
             << count[i] << endl;
    }
}


void GenericProbabilisticAlignment::read(istream &input)
{
    string TRANSITION        = "TransitionFeature";
    string READHOMOPOLYMER   = "ReadHomopolymerFeature";
    string GENOMEHOMOPOLYMER = "GenomeHomopolymerFeature";
    string EMISSION          = "EmissionFeature";

    string line;
    while(!input.eof() && input.good())
    {
        line.clear();
        getline(input, line);

        if (line.empty())
            continue;

        stringstream buffer;
        buffer << line;

        string       featType;
        int          featID;
        string       featCurrState;
        string       featPrevState;
        string       featHomoLen;
        string       featHomoDelta;
        string       featEmisAlpha;
        string       featEmisBeta;
        double       featWeight;

        buffer >> featType;

        // transition feature
        if (featType==TRANSITION)
        {
            buffer >> featID;
            buffer >> featPrevState;
            buffer >> featCurrState;
            buffer >> featWeight;

            ProbalnTransitionFeature f = TransitionFeatureSet[featID];
            TransitionFeatureWeight[f] = featWeight;
        }

        // read homopolymer feature
        if (featType==READHOMOPOLYMER)
        {
            buffer >> featID;
            buffer >> featCurrState;
            buffer >> featHomoLen;
            buffer >> featHomoDelta;
            buffer >> featWeight;

            ProbalnReadHomopolymerFeature f = ReadHomopolymerFeatureSet[featID];
            ReadHomopolymerFeatureWeight[f] = featWeight;
        }

        // genome homopolymer feature
        if (featType==GENOMEHOMOPOLYMER)
        {
            buffer >> featID;
            buffer >> featCurrState;
            buffer >> featHomoLen;
            buffer >> featHomoDelta;
            buffer >> featWeight;

            ProbalnGenomeHomopolymerFeature f = GenomeHomopolymerFeatureSet[featID];
            GenomeHomopolymerFeatureWeight[f] = featWeight;
        }

        // emission feature
        if (featType==EMISSION)
        {
            buffer >> featID;
            buffer >> featCurrState;
            buffer >> featEmisAlpha;
            buffer >> featEmisBeta;
            buffer >> featWeight;

            ProbalnEmissionFeature f = EmissionFeatureSet[featID];
            EmissionFeatureWeight[f] = featWeight;
        }
    }
}
