#include "GenericIndividualSnpCall.h"
#include "GenericBamAlignmentTools.h"
#include "GenericFastaTools.h"
#include "GenericGraphTools.h"
using namespace GenericSequenceTools;

#include "api/BamAux.h"
#include "api/BamAlignment.h"
#include "utils/bamtools_pileup_engine.h"
using namespace BamTools;

#include <set>
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <random>
#include <chrono>
#include <ctime>
#include <functional>
using namespace std;

GenericIndividualSnpCall::GenericIndividualSnpCall(int downSample, int minReadLength, int minMapQuality, double maxMismatchFrac, int minSnpRead, int minSnpFrac, int verbosity)
    : m_minReadLength(minReadLength)
    , m_minMapQuality(minMapQuality)
    , m_maxMismatchFrac(maxMismatchFrac)
    , m_downSample(downSample)
    , m_minSnpRead(minSnpRead)
    , m_minSnpFrac(minSnpFrac)
    , m_verbosity(verbosity)
{
}

int setupGenomeBlock(RefVector& chromosomes, BamRegion& roi, vector<int>& BlockChrID, vector<int>& BlockLeftPos, vector<int>& BlockRightPos)
{
    int BlockSize = 1e6;
    // ROI is a chromosome
    if (roi.LeftRefID==roi.RightRefID)
    {
        int regionSize = roi.RightPosition-roi.LeftPosition;

        if (regionSize>BlockSize)
        {
            for (int pos=roi.LeftPosition; pos<roi.RightPosition; pos+=BlockSize)
            {
                BlockChrID.push_back(roi.LeftRefID);
                BlockLeftPos.push_back(pos);
                if (pos+BlockSize<roi.RightPosition)
                {
                    BlockRightPos.push_back(pos+BlockSize);
                }else
                {
                    BlockRightPos.push_back(roi.RightPosition);
                }
            }
        }else
        {
            BlockChrID.push_back(roi.LeftRefID);
            BlockLeftPos.push_back(roi.LeftPosition);
            BlockRightPos.push_back(roi.RightPosition);
        }
    }

    // ROI is a set of chromosomes
    if (roi.LeftRefID<roi.RightRefID)
    {
        for (int chrID=roi.LeftRefID; chrID<=roi.RightRefID; ++chrID)
        {
            int chrLeftPos, chrRightPos;
            if (chrID==roi.LeftRefID)
            {
                chrLeftPos = roi.LeftPosition;
                chrRightPos = chromosomes[roi.LeftRefID].RefLength;
            }else if (chrID==roi.RightRefID)
            {
                chrLeftPos = 0;
                chrRightPos = roi.RightPosition;
            }else
            {
                chrLeftPos = 0;
                chrRightPos = chromosomes[chrID].RefLength;
            }

            int regionSize = chrRightPos - chrLeftPos;
            if (regionSize>BlockSize)
            {
                for (int pos=chrLeftPos; pos<chrRightPos; pos+=BlockSize)
                {
                    BlockChrID.push_back(chrID);
                    BlockLeftPos.push_back(pos);
                    if (pos+BlockSize<chrRightPos)
                    {
                        BlockRightPos.push_back(pos+BlockSize);
                    }else
                    {
                        BlockRightPos.push_back(chrRightPos);
                    }
                }
            }else
            {
                BlockChrID.push_back(chrID);
                BlockLeftPos.push_back(chrLeftPos);
                BlockRightPos.push_back(chrRightPos);
            }
        }
    }

    return BlockChrID.size();
}

struct SortAlleleByPosition
{
    bool operator()(const Allele& a, const Allele& b)
    {
        return (a.m_chrPosition < b.m_chrPosition);
    }
};

void mergeSnpSitesToBlocks(vector<Allele>& SnpSites, vector<tuple<int,int, list<Allele>>>& SnpBlocks)
{
    // sort by position
    sort(SnpSites.begin(), SnpSites.end(), SortAlleleByPosition());

    SnpBlocks.clear();
    for (int i=0; i<SnpSites.size(); i++)
    {
        if (i==0)
        {
            int pos = SnpSites[i].m_chrPosition;
            SnpBlocks.push_back(tuple<int,int,list<Allele>>(pos, pos+1, list<Allele>(1, SnpSites[i])));
            continue;
        }

        int j = SnpBlocks.size()-1;
        int SnpBlockLeftPos  = get<0>(SnpBlocks[j]);
        int SnpBlockRightPos = get<1>(SnpBlocks[j]);

        int pos = SnpSites[i].m_chrPosition;

        if ((pos-SnpBlockRightPos)<0.1*(pos-SnpBlockLeftPos))
        {
            get<1>(SnpBlocks[j]) = pos+1;
            get<2>(SnpBlocks[j]).emplace_back(SnpSites[i]);
            continue;
        }

        SnpBlocks.push_back(tuple<int,int,list<Allele>>(pos, pos+1, list<Allele>(1, SnpSites[i])));
    }
}

int GenericIndividualSnpCall::call(Fasta &fastaObj, BamReader &bamObj, BamRegion &roi, GenericProbabilisticAlignment &probAligner, VariantCallSetting& snpCallSettings, vector<GenericVariant> &variantSet)
{
    RefVector chromosomes = bamObj.GetReferenceData();
    // set up genome blocks
    vector<int> BlockChrID, BlockLeftPos, BlockRightPos;
    int BlockNumber=setupGenomeBlock(chromosomes, roi, BlockChrID, BlockLeftPos, BlockRightPos);

    int numSNP = 0;

    // iterate throught blocks
    for (int i=0; i<BlockNumber; ++i)
    {
        if (m_verbosity>=1)
        {
            cout << "processing " << chromosomes[BlockChrID[i]].RefName << ":" << BlockLeftPos[i]+1 << "-" << BlockRightPos[i] << endl;
        }

        clock_t startTime = clock();

        // genome
        string BlockGenome;
        fastaObj.GetSequence(BlockChrID[i], BlockLeftPos[i], BlockRightPos[i], BlockGenome);

        map<int,list<tuple<char,int,int,double>>> BlockBamData;
        AlleleSet BlockSnpAlleleCandidates;
        // profile SNP sites by the simple method
        simpleSnpCall(BlockGenome, bamObj, BlockChrID[i], BlockLeftPos[i], BlockRightPos[i], BlockSnpAlleleCandidates, BlockBamData);

        // merge SNP sites to SNP blocks
        vector<tuple<int,int,list<Allele>>> BlockSnpLoci;
        mergeSnpSitesToBlocks(BlockSnpAlleleCandidates, BlockSnpLoci);

        // iterate through Snp locus
        for (int j=0; j<BlockSnpLoci.size(); j++)
        {
            int BlockSnpLeftPos  = get<0>(BlockSnpLoci[j]);
            int BlockSnpRightPos = get<1>(BlockSnpLoci[j]);

            // it is a SNP site
            if (BlockSnpRightPos==BlockSnpLeftPos+1)
            {
                simpleBayesianSnpCall(fastaObj, bamObj, BlockChrID[i], BlockSnpLeftPos, BlockSnpRightPos, get<2>(BlockSnpLoci[j]), BlockBamData[BlockSnpLeftPos], snpCallSettings, variantSet);
            }else if (BlockSnpRightPos==BlockSnpLeftPos+2)
            {
                for (int pos=BlockSnpLeftPos; pos<BlockSnpRightPos; pos++)
                {
                    list<Allele> fAlleles = get<2>(BlockSnpLoci[j]);
                    list<Allele> tAlleles;
                    for (list<Allele>::iterator faIter=fAlleles.begin(); faIter!=fAlleles.end(); faIter++)
                    {
                        if (faIter->m_chrPosition==pos)
                            tAlleles.emplace_back(*faIter);
                    }

                    if (!tAlleles.empty())
                        simpleBayesianSnpCall(fastaObj, bamObj, BlockChrID[i], pos, pos+1, tAlleles, BlockBamData[pos], snpCallSettings, variantSet);

                }
            }
            else   // it is a MNP site
            {
                PyroHMMsnp(fastaObj, bamObj, BlockChrID[i], BlockSnpLeftPos, BlockSnpRightPos, probAligner, get<2>(BlockSnpLoci[j]), snpCallSettings, variantSet);
            }
        }

        clock_t endTime = clock();
        if (m_verbosity>=1)
        {
            cout << "time elapsed " << ((endTime-startTime)/(double)CLOCKS_PER_SEC/60.) << " minutes";
            cout << ", ";
            cout << "call " << variantSet.size()-numSNP << " SNPs" << endl;
        }

        numSNP = variantSet.size();
    }

    return variantSet.size();
}


class SimpleSnpCallPileupVisitor : public PileupVisitor
{
    public:
        SimpleSnpCallPileupVisitor(string* fastaObj, int chrID, int chrLeftPos, int chrRightPos, int downsample, set<int>* SnpPositions, vector<Allele>* SnpAlleles, map<int,list<tuple<char,int,int,double>>>* BamData)
            : PileupVisitor()
            , m_fasta(fastaObj)
            , m_snpPositions(SnpPositions)
            , m_snpAlleles(SnpAlleles)
            , m_bamData(BamData)
            , m_chrID(chrID)
            , m_chrLeftPos(chrLeftPos)
            , m_chrRightPos(chrRightPos)
            , m_downSample(downsample)
        {
        }
    public:
        void Visit(const PileupPosition& pileupData);

    public:
        string*         m_fasta;
        set<int>*       m_snpPositions;
        vector<Allele>* m_snpAlleles;
        map<int,list<tuple<char,int,int,double>>>* m_bamData;
        int             m_chrID;
        int             m_chrLeftPos;
        int             m_chrRightPos;
        int             m_downSample;
};

void SimpleSnpCallPileupVisitor::Visit(const PileupPosition &pileupData)
{

    if (pileupData.Position<m_chrLeftPos || pileupData.Position>=m_chrRightPos)
        return;

    if (m_snpPositions->find(pileupData.Position)==m_snpPositions->end())
        return;

    // random dicer
    unsigned UniformSeed = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 UniformGenerator(UniformSeed);
    uniform_int_distribution<int> UniformDist;

    UniformDist = uniform_int_distribution<int>(0, pileupData.PileupAlignments.size());
    auto UniformDicer = bind(UniformDist, UniformGenerator);

    // reference base
    char refBase = m_fasta->at(pileupData.Position-m_chrLeftPos);

    if (refBase=='N')
        return;

    // information buffer
    int globalDepth         = 0;
    double globalMapAvgQual = 0;
    int globalStrandPos     = 0;
    int globalStrandNeg     = 0;

    list<tuple<char,int,int,double>> bamData;
    set<char> alleleBases;
    map<char, int> alleleDepth;
    map<char, double> alleleMapAvgQual;
    map<char, int> alleleStrandPos;
    map<char, int> alleleStrandNeg;
    // iterate through alignments
    for (int i=0; i<pileupData.PileupAlignments.size(); i++)
    {
        char readBase;
        BamAlignment al = pileupData.PileupAlignments[i].Alignment;

        // update global info
        globalDepth         += 1;
        globalMapAvgQual    += al.MapQuality*al.MapQuality;
        if (al.IsReverseStrand())
            globalStrandNeg += 1;
        else
            globalStrandPos += 1;

        // downsampling
        if (UniformDicer()>m_downSample)
            continue;

        // data info
        int nx = GenericBamAlignmentTools::numBamAlignmentMismatches(al);
        int nm = GenericBamAlignmentTools::numBamAlignmentMatches(al);
        double nr = nm/(nm+nx+0.);

        if (pileupData.PileupAlignments[i].IsCurrentDeletion)
        {
            bamData.push_back(tuple<char,int,int,double>('-',0,al.MapQuality,1.));
        }else
        {
            readBase = al.QueryBases.at(pileupData.PileupAlignments[i].PositionInAlignment);
            int baseQuality = int(al.Qualities.at(pileupData.PileupAlignments[i].PositionInAlignment))-33;

            // debug
            if (baseQuality<0)
            {
                cout << int(al.Qualities.at(pileupData.PileupAlignments[i].PositionInAlignment)) << endl;
            }

            if (readBase==refBase || readBase=='N')
            {
                bamData.push_back(tuple<char,int,int,double>(readBase, baseQuality, al.MapQuality,1.));
            }
            else
            {
                bamData.push_back(tuple<char,int,int,double>(readBase, baseQuality, al.MapQuality,nr));
            }
        }

        // here is deletion
        if (pileupData.PileupAlignments[i].IsCurrentDeletion)
            continue;

        // update allele info
        if (readBase!='N' && readBase!=refBase)
        {
            if (pileupData.PileupAlignments[i].IsSegmentBegin || pileupData.PileupAlignments[i].IsSegmentEnd)
                continue;

            if (alleleBases.find(readBase)==alleleBases.end())
            {
                alleleDepth[readBase]         = 1;
                alleleMapAvgQual[readBase]    = al.MapQuality*al.MapQuality;
                if (al.IsReverseStrand())
                {
                    alleleStrandPos[readBase] = 0;
                    alleleStrandNeg[readBase] = 1;
                }
                else
                {
                    alleleStrandPos[readBase] = 1;
                    alleleStrandNeg[readBase] = 0;
                }
                alleleBases.insert(readBase);
            }else
            {
                alleleDepth[readBase]         += 1;
                alleleMapAvgQual[readBase]    += al.MapQuality*al.MapQuality;
                if (al.IsReverseStrand())
                    alleleStrandNeg[readBase] += 1;
                else
                    alleleStrandPos[readBase] += 1;
            }
        }
    }

    // save allele info
    for (set<char>::iterator iter=alleleBases.begin(); iter!=alleleBases.end(); iter++)
    {
        Allele allele;
        allele.m_chrID            = pileupData.RefId;
        allele.m_chrPosition      = pileupData.Position;
        allele.m_allele           = *iter;
        allele.m_reference        = refBase;
        allele.m_alleleDepth      = alleleDepth[*iter];
        allele.m_alleleMapAvgQual = sqrt(alleleMapAvgQual[*iter])/allele.m_alleleDepth;
        allele.m_alleleStrandPos  = alleleStrandPos[*iter];
        allele.m_alleleStrandNeg  = alleleStrandNeg[*iter];

        // global info
        allele.m_globalDepth      = globalDepth;
        allele.m_globalMapAvgQual = sqrt(globalMapAvgQual)/globalDepth;
        allele.m_globalStrandPos  = globalStrandPos;
        allele.m_globalStrandNeg  = globalStrandNeg;

        m_snpAlleles->push_back(allele);
    }

    (*m_bamData)[pileupData.Position] = bamData;
}

void GenericIndividualSnpCall::simpleSnpCall(string &fastaObj, BamReader &bamObj, int chrID, int leftPosition, int rightPosition, vector<Allele> &variantCandidates, map<int,list<tuple<char,int,int,double>>> &bamData)
{
    set<int> BlockSnpPositions;
    vector<Allele> BlockSnpAlleles;

    // rewind
    bamObj.Rewind();
    // set region
    bamObj.SetRegion(chrID, leftPosition, chrID, rightPosition);

    BamAlignment al;
    // search SNP positions in the region
    while (bamObj.GetNextAlignment(al))
    {
        if (!GenericBamAlignmentTools::goodAlignment(al))
            continue;

        if (!al.HasTag("MD"))
            continue;

        vector<long> SnpInAlignment;
        GenericBamAlignmentTools::getBamAlignmentMismatches(al, SnpInAlignment);

        for (int i=0; i<SnpInAlignment.size(); i++)
        {
            BlockSnpPositions.insert(SnpInAlignment[i]);
        }
    }

    // pileup visitor
    SimpleSnpCallPileupVisitor visitor(&fastaObj, chrID, leftPosition, rightPosition, m_downSample, &BlockSnpPositions, &BlockSnpAlleles, &bamData);

    PileupEngine SimpleSnpCallPileupEngine;
    SimpleSnpCallPileupEngine.AddVisitor(&visitor);

    // rewind
    bamObj.Rewind();
    // set region
    bamObj.SetRegion(chrID, leftPosition, chrID, rightPosition);
    // load data
    while(bamObj.GetNextAlignment(al))
    {

        if (!GenericBamAlignmentTools::goodAlignment(al))
            continue;

        if (!GenericBamAlignmentTools::validMapQuality(al, m_minMapQuality))
            continue;

        if (!GenericBamAlignmentTools::validReadIdentity(al, m_maxMismatchFrac))
            continue;

        if (!GenericBamAlignmentTools::validReadLength(al, m_minReadLength))
            continue;

        if (!al.HasTag("MD"))
            continue;

        SimpleSnpCallPileupEngine.AddAlignment(al);
    }
    SimpleSnpCallPileupEngine.Flush();

    // Filter SNP candidiate
    for (int i=0; i<BlockSnpAlleles.size(); i++)
    {
        Allele allele = BlockSnpAlleles[i];

        if (allele.m_alleleDepth < m_minSnpRead)
            continue;
        if (allele.m_alleleDepth < m_minSnpFrac*allele.m_globalDepth)
            continue;

        variantCandidates.push_back(allele);
    }

}

void simpleBayesianSnpCallPrior(VariantCallSetting& snpCallSettings, int numAlleles, vector<long double>& allelePriors)
{
    if (snpCallSettings.m_priorType==SNPCALL_PRIOR_FLAT)
    {
        for (int i=0; i<numAlleles; i++)
            allelePriors[i] = 1./numAlleles;
    }else
    {
        allelePriors[0] = 1-snpCallSettings.m_prior;
        for (int i=1; i<numAlleles; i++)
            allelePriors[i] = snpCallSettings.m_prior/(numAlleles-1);
    }
}

void simpleBayesianSnpCallAlleleLikelihood(vector<string>& alleleChars, int numAlleles, list<tuple<char,int,int,double>>& data, vector<vector<long double>>& alleleDataLikelihoods)
{
    list<tuple<char,int,int,double>>::iterator dataIter;
    for (int i=0; i<numAlleles; i++)
    {
        char alleleChar = alleleChars[i][0];
        for (dataIter=data.begin(); dataIter!=data.end(); dataIter++)
        {
            double nr = get<3>(*dataIter);

            double dataQual, dataProb, dataLike;
            if (get<0>(*dataIter)=='-')
            {
                dataQual = get<2>(*dataIter);
                dataQual *= nr;

                dataProb = 1-pow(10,-0.1*dataQual);
                if (dataProb<1e-6)
                    dataProb = 1e-6;

                dataLike = pow(10,-0.1*dataQual);
            }else
            {
                dataQual = get<1>(*dataIter);
                if (dataQual>get<2>(*dataIter))
                    dataQual = get<2>(*dataIter);
                dataQual *= nr;

                dataProb = 1-pow(10,-0.1*dataQual);
                if (dataProb<1e-6)
                    dataProb = 1e-6;

                if (get<0>(*dataIter)=='N')
                {
                    dataLike = dataProb;
                }else if (get<0>(*dataIter)==alleleChar)
                {
                    dataLike = dataProb;
                }else
                {
                    dataLike = pow(10,-0.1*dataQual);
                }
            }

            alleleDataLikelihoods[i].push_back(dataLike);
        }
    }
}

void simpleBayesianSnpCallGenotypeSet(int ploidy, int allele, int numAlleles, vector<int>& precedeAlleles, vector<vector<int>>& genotypes, set<set<int>>& genotypeDiscovered)
{
    set<int> key;
    for (int i=0; i<precedeAlleles.size(); i++)
    {
        key.insert(precedeAlleles[i]);
    }
    key.insert(allele);

    if (genotypeDiscovered.find(key)!=genotypeDiscovered.end())
        return;

    if (ploidy==1)
    {
        vector<int> G(precedeAlleles);
        G.push_back(allele);
        genotypes.push_back(G);

    }else
    {
        for (int i=0; i<numAlleles; i++)
        {
            vector<int> localPrecedeAlleles(precedeAlleles);
            localPrecedeAlleles.push_back(allele);
            simpleBayesianSnpCallGenotypeSet(ploidy-1, i, numAlleles, localPrecedeAlleles, genotypes, genotypeDiscovered);
        }
    }

    genotypeDiscovered.insert(key);
}

void simpleBayesianSnpCallGenotypePrior(VariantCallSetting& snpCallSettings, vector<vector<int>>& genotypes, vector<long double>& genotypePriors)
{
    if (snpCallSettings.m_priorType==SNPCALL_PRIOR_FLAT)
    {
        int numGenotypes = genotypes.size();
        for (int i=0; i<numGenotypes; i++)
            genotypePriors[i]=(1./numGenotypes);
    }else
    {
        int numVariantGenotypes = genotypes.size()-1;
        genotypePriors[0] = (1-snpCallSettings.m_prior);
        for (int i=1; i<numVariantGenotypes+1; i++)
        {
            genotypePriors[i]=(snpCallSettings.m_prior/numVariantGenotypes);
        }
    }
}

void simpleBayesianSnpCallGenotypeLikelihood(int ploidy, int datasize, vector<vector<int>>& genotypes, vector<vector<long double>>& alleleDataLikelihoods, vector<long double>& genotypeLikelihoods)
{
    int numGenotype = genotypes.size();
    for (int i=0; i<numGenotype; i++)
    {
        vector<int> G = genotypes[i];

        genotypeLikelihoods[i] = 0;

        for (int j=0; j<datasize; j++)
        {
            long double genotypeDataLike = 0;
            for (int k=0; k<ploidy; k++)
            {
                genotypeDataLike += alleleDataLikelihoods[G[k]][j]/ploidy;
            }

            genotypeLikelihoods[i] += log(genotypeDataLike);
        }

        genotypeLikelihoods[i] = exp(genotypeLikelihoods[i]);
    }
}

void simpleBayesianSnpCallGenotypePosterior(int numGenotypes, vector<long double>& genotypePriors, vector<long double>& genotypeLikelihoods, vector<long double>& genotypePosteriors)
{
    long double Z = 0;
    for (int i=0; i<numGenotypes; i++)
    {
        genotypePosteriors[i] = genotypePriors[i]*genotypeLikelihoods[i];
        Z += genotypePosteriors[i];
    }
    for (int i=0; i<numGenotypes; i++)
    {
        genotypePosteriors[i] /= Z;
    }
}

void GenericIndividualSnpCall::simpleBayesianSnpCall(Fasta &fastaObj, BamReader &bamObj, int chrID, int leftPosition, int rightPosition, list<Allele>& alleles, list<tuple<char,int,int,double>>& data, VariantCallSetting& snpCallSettings, vector<GenericVariant> &variantResult)
{
    // number of alleles, include reference
    int numAlleles = 1+alleles.size();

    // array of alleles
    vector<string> alleleChars;
    list<Allele>::iterator alleleListIter = alleles.begin();
    // reference allele
    alleleChars.push_back(alleleListIter->m_reference);
    // variant allele
    for ( ; alleleListIter!=alleles.end(); alleleListIter++)
    {
        alleleChars.push_back(alleleListIter->m_allele);
    }

    // allele information
    vector<Allele> allelePool;
    Allele refAllele;
    refAllele.m_chrID       = chrID;
    refAllele.m_chrPosition = leftPosition;
    refAllele.m_allele      = alleleChars[0];
    allelePool.push_back(refAllele);
    for (list<Allele>::iterator alleleIter=alleles.begin(); alleleIter!=alleles.end(); alleleIter++)
    {
        allelePool.push_back(*alleleIter);
    }

    // likelihood of alleles
    vector<vector<long double>> alleleDataLikelihoods(numAlleles);
    simpleBayesianSnpCallAlleleLikelihood(alleleChars, numAlleles, data, alleleDataLikelihoods);

    // genotypes
    vector<vector<int>> genotypes;
    set<set<int>> genotypeDiscovered;
    for (int i=0; i<numAlleles; i++)
    {
        vector<int> precedeAlleles;
        simpleBayesianSnpCallGenotypeSet(snpCallSettings.m_ploidy, i, numAlleles, precedeAlleles, genotypes, genotypeDiscovered);
    }

    int numGenotype = genotypes.size();

    // genotype prior
    vector<long double> genotypePriors(numGenotype);
    simpleBayesianSnpCallGenotypePrior(snpCallSettings, genotypes, genotypePriors);

    // genotype likelihood
    vector<long double> genotypeLikelihoods(numGenotype);
    simpleBayesianSnpCallGenotypeLikelihood(snpCallSettings.m_ploidy, data.size(), genotypes, alleleDataLikelihoods, genotypeLikelihoods);

    // genotype posterior
    vector<long double> genotypePosteriors(numGenotype);
    simpleBayesianSnpCallGenotypePosterior(numGenotype, genotypePriors, genotypeLikelihoods, genotypePosteriors);

    // search maximal posterior genotype
    long double maxGenotypePosterior = 0;
    int inferGenotype;
    for (int i=1; i<numGenotype; i++)
    {
        if (maxGenotypePosterior<genotypePosteriors[i])
        {
            maxGenotypePosterior = genotypePosteriors[i];
            inferGenotype = i;
        }
    }

    // save result
    GenericVariant result;

    vector<int> G = genotypes[inferGenotype];
    for (int i=0; i<G.size(); i++)
    {
        result.m_alleles.push_back(allelePool[G[i]]);
    }

    result.m_chrID        = chrID;
    result.m_chrPosition  = leftPosition;
    result.m_variantType  = VARIANT_SNP;
    result.m_probScoreRef = genotypePosteriors[0];
    result.m_probScoreVar = maxGenotypePosterior;
    if (fabs(result.m_probScoreRef-1)<1e-300)
        result.m_quality  = 0;
    else if (result.m_probScoreRef<1e-300)
        result.m_quality  = 3000;
    else
        result.m_quality  = -10*log10(result.m_probScoreRef);
    result.m_reference    = alleleChars[0];

    for (int i=0; i<result.m_alleles.size(); i++)
    {
        if (result.m_alleles[i].m_allele==result.m_reference)
            result.m_haploidType.push_back(0);
        else
            result.m_haploidType.push_back(1);
    }


    // filter
    if (result.m_quality>=snpCallSettings.m_variantQualityFilter)
        variantResult.push_back(result);
}



struct PyroHMMsnp_Sequence_t
{
    string          t_ID;
    string          t_sequence;
    Cigar           t_cigar;
    BamMD           t_md;
    int             t_numMismatch;
    int             t_numInDel;
    int             t_length;
    int             t_type;
    string          t_alnLocalRead;
    string          t_alnLocalGenome;
    int             t_mapQualScore;
    int             t_startPositionShift;
    int             t_endPositionShift;
};

void PyroHMMsnpHaplotypeDataLikelihood(GenericProbabilisticAlignment& ProbAligner, int band, int numHaplotypes, vector<string>& haplotypeSeqs, vector<PyroHMMsnp_Sequence_t>& data, vector<vector<long double>>& haplotypeDataLikelihoods)
{
    for (int i=0; i<numHaplotypes; i++)
    {
        string haplotype = haplotypeSeqs[i];

        for (int j=0; j<data.size(); j++)
        {
            string read = data[j].t_sequence;
            long double viterbiScore;
            ProbAligner.ViterbiComputation(read, haplotype, band, *static_cast<string*>(0), *static_cast<string*>(0), viterbiScore);
            haplotypeDataLikelihoods[i].push_back(exp(viterbiScore/read.length()));
        }
    }
}

void PyroHMMsnpGenotypeSet(int ploidy, int haplotype, int numHaplotypes, vector<int>& precedeHaplotypes, vector<vector<int>>& genotypes, set<set<int>>& genotypeDiscovered)
{
    set<int> key;
    for (int i=0; i<precedeHaplotypes.size(); i++)
    {
        key.insert(precedeHaplotypes[i]);
    }
    key.insert(haplotype);

    if (genotypeDiscovered.find(key)!=genotypeDiscovered.end())
        return;

    if (ploidy==1)
    {
        vector<int> G(precedeHaplotypes);
        G.push_back(haplotype);
        genotypes.push_back(G);

    }else
    {
        for (int i=0; i<numHaplotypes; i++)
        {
            vector<int> localPrecedeHaplotypes(precedeHaplotypes);
            localPrecedeHaplotypes.push_back(haplotype);
            PyroHMMsnpGenotypeSet(ploidy-1, i, numHaplotypes, localPrecedeHaplotypes, genotypes, genotypeDiscovered);
        }
    }

    genotypeDiscovered.insert(key);
}


void PyroHMMsnpGenotypePrior(int numGenotypes, vector<vector<int>>& genotypes, VariantCallSetting& snpCallSettings, vector<long double>& genotypePriors)
{
    if (snpCallSettings.m_priorType==SNPCALL_PRIOR_FLAT)
    {
        for (int i=0; i<numGenotypes; i++)
        {
            genotypePriors[i] = 1./numGenotypes;
        }
    }else
    {
        genotypePriors[0] = (1-snpCallSettings.m_prior);
        for (int i=1; i<numGenotypes; i++)
        {
            genotypePriors[i] = snpCallSettings.m_prior / (numGenotypes-1.0);
        }
    }
}

void PyroHMMsnpGenotypeLikelihood(int numGenotypes, vector<vector<int>>& genotypes, int datasize, vector<vector<long double>>& haplotypeDataLikehoods, VariantCallSetting& snpCallingSettings, vector<long double>& genotypeLikehoods)
{
    for (int i=0; i<numGenotypes; i++)
    {
        vector<int> haplotypes = genotypes[i];

        long double genotypeLikehood = 0;
        for (int k=0; k<datasize; k++)
        {
            long double dataLikelihood = 0;
            for (int j=0; j<snpCallingSettings.m_ploidy; j++)
            {
                dataLikelihood += haplotypeDataLikehoods[haplotypes[j]][k]/snpCallingSettings.m_ploidy;
            }

            genotypeLikehood += log(dataLikelihood);
        }

        genotypeLikehoods[i] = exp(genotypeLikehood);
    }
}

void PyroHMMsnpGenotypePosterior(int numGenotypes, vector<long double>& genotypePriors, vector<long double>& genotypeLikelihoods, vector<long double>& genotypePosteriors)
{
    long double Z = 0;
    for (int i=0; i<numGenotypes; i++)
    {
        genotypePosteriors[i] = genotypePriors[i]*genotypeLikelihoods[i];
        Z += genotypePosteriors[i];
    }
    for (int i=0; i<numGenotypes; i++)
    {
        genotypePosteriors[i] /= Z;
    }
}

void GenericIndividualSnpCall::PyroHMMsnp(Fasta &fastaObj, BamReader &bamObj, int chrID, int leftPosition, int rightPosition, GenericProbabilisticAlignment &probAligner, list<Allele>& allelesInBlock, VariantCallSetting& snpCallSettings, vector<GenericVariant> &variantResults)
{
    VariantCallSetting settingForPyroHMMsnp = snpCallSettings;

    // allele pool
    vector<Allele> allelePool;
    for (list<Allele>::iterator allelesInBlockIter=allelesInBlock.begin(); allelesInBlockIter!=allelesInBlock.end(); allelesInBlockIter++)
    {
        allelePool.push_back(*allelesInBlockIter);
    }

    // add 10bp flanking segment at each side
    int windowLeftPosition  = leftPosition  - snpCallSettings.m_flankingSize;
    int windowRightPosition = rightPosition + snpCallSettings.m_flankingSize;

    // genome
    string genome;
    fastaObj.GetSequence(chrID, windowLeftPosition, windowRightPosition, genome);

    int    globalDepth;
    double globalMapQual;
    int    globalStrandPos;
    int    globalStrandNeg;

    vector<PyroHMMsnp_Sequence_t> readsInWindow;

    // rewind BAM reader
    bamObj.Rewind();
    // set BAM region
    bamObj.SetRegion(chrID, windowLeftPosition, chrID, windowRightPosition);
    // read alignment
    BamAlignment al;
    while (bamObj.GetNextAlignment(al))
    {
        // skip if it is not a good alignment
        if (!GenericBamAlignmentTools::goodAlignment(al))
        {
            continue;
        }

        // skip if it is not valid at length
        if (!GenericBamAlignmentTools::validReadLength(al, m_minReadLength))
        {
            continue;
        }

        // skip if it is not valid at map quality
        if (!GenericBamAlignmentTools::validMapQuality(al, m_minMapQuality))
        {
            continue;
        }

        // skip if it is not valid at alignment identity
        if (!GenericBamAlignmentTools::validReadIdentity(al, m_maxMismatchFrac))
        {
            continue;
        }

        // global info
        globalDepth   += 1;
        globalMapQual += al.MapQuality*al.MapQuality;
        if (al.IsReverseStrand())
            globalStrandNeg += 1;
        else
            globalStrandPos += 1;

        // get local alignment
        string t_localRead, t_localGenome;
        Cigar  t_cigar;
        BamMD  t_md;
        int    t_numMismatch, t_numInDel;
        GenericBamAlignmentTools::getLocalAlignment(al, windowLeftPosition, windowRightPosition-windowLeftPosition,
                                                    t_localRead, t_localGenome, t_cigar, t_md,
                                                    t_numMismatch, t_numInDel);

        if (t_localRead.empty() || t_localGenome.empty())
            continue;


        // save into set
        PyroHMMsnp_Sequence_t t_seq;
        t_seq.t_ID           = GenericBamAlignmentTools::getBamAlignmentID(al);
        t_seq.t_sequence     = t_localRead;
        t_seq.t_cigar        = t_cigar;
        t_seq.t_md           = t_md;
        t_seq.t_numMismatch  = t_numMismatch;
        t_seq.t_numInDel     = t_numInDel;
        t_seq.t_mapQualScore = al.MapQuality;


        if (al.Position>windowLeftPosition)
            t_seq.t_startPositionShift = al.Position-windowLeftPosition;
        else
            t_seq.t_startPositionShift = 0;

        if (al.GetEndPosition()<windowRightPosition)
            t_seq.t_endPositionShift = windowRightPosition-al.GetEndPosition();
        else
            t_seq.t_endPositionShift = 0;

        readsInWindow.push_back(t_seq);
    }

    int numData = readsInWindow.size();

    // construct the consensus sequence graph
    GenericDagGraph consensusGraph;
    vector<string>  consensusGraphReads;
    vector<Cigar>   consensusGraphReadCigars;
    vector<int>     consensusGraphReadStarts;

    // set of aligned reads to construct the graph
    for (int i=0; i<numData; ++i)
    {
        consensusGraphReads.push_back(readsInWindow[i].t_sequence);
        consensusGraphReadCigars.push_back(readsInWindow[i].t_cigar);
        consensusGraphReadStarts.push_back(readsInWindow[i].t_startPositionShift);
    }

    // build up the graph
    consensusGraph.buildDagGraph(genome, consensusGraphReads, consensusGraphReadCigars, consensusGraphReadStarts);
    consensusGraph.edgePruning(snpCallSettings.m_graphPruneLevel);

    // search topK paths, excluding reference
    vector<string>       topRankConsensusGraphPaths;
    vector<list<Vertex>> topRankConsensusGraphPathVertexs;
    vector<double>       topRankConsensusGraphPathWeights;
    consensusGraph.topRankPathsExcludeGenome(30, topRankConsensusGraphPaths, topRankConsensusGraphPathVertexs, topRankConsensusGraphPathWeights);

    // change vertex list to vertex set
    vector<set<Vertex>>  topRankConsensusGraphPathVertexSet;
    for (int i=0; i<topRankConsensusGraphPathVertexs.size(); i++)
    {
        list<Vertex>::iterator vertexIter = topRankConsensusGraphPathVertexs[i].begin();
        set<Vertex> vertexSet;
        for (; vertexIter!=topRankConsensusGraphPathVertexs[i].end(); vertexIter++)
        {
            vertexSet.insert(*vertexIter);
        }
        topRankConsensusGraphPathVertexSet.push_back(vertexSet);
    }

    // get variant vertices
    vector<int>    allelePositions;
    vector<string> alleleChars;
    for (list<Allele>::iterator alleleIter=allelesInBlock.begin(); alleleIter!=allelesInBlock.end(); alleleIter++)
    {
        Allele allele = *alleleIter;
        allelePositions.push_back(allele.m_chrPosition-windowLeftPosition);
        alleleChars.push_back(allele.m_allele);
    }
    // map allele to graph vertex
    set<Vertex> variantVertexs;
    map<int,Vertex> mapAlleleToVertex;
    map<Vertex,int> mapVertexToAllele;
    for (int v=0; v<consensusGraph.m_numVertexs; v++)
    {
        if (consensusGraph.m_skip[v])
            continue;

        if (!consensusGraph.m_isMismatch[v])
            continue;

        int gp = consensusGraph.m_genomePosition[v] - 1;


        for (int j=0; j<allelePool.size(); j++)
        {
            int ap = allelePositions[j];
            if (ap==gp)
            {
                if (alleleChars[j]==consensusGraph.m_labels[v])
                {
                    variantVertexs.insert(v);
                    mapAlleleToVertex[j] = v;
                    mapVertexToAllele[v] = j;
                }
            }
        }
    }


    // set up the haplotypes
    vector<string> haplotypes;
    vector<int>    haplotypeToPathIndex;
    vector<set<Vertex>> haplotypeVariantVertexs;

    haplotypes.push_back(genome);
    haplotypeToPathIndex.push_back(-1);
    haplotypeVariantVertexs.push_back(set<Vertex>());

    int kk = 0;
    for (int i=0; i<topRankConsensusGraphPaths.size(); i++)
    {
        if (kk>=snpCallSettings.m_topK)
            continue;

        bool hasVariantVertex = false;
        int  deltaLength = (topRankConsensusGraphPaths[i].length()-genome.length());
        deltaLength = abs(deltaLength);

        if (deltaLength>5)
            continue;

        set<Vertex> pathVertexs = topRankConsensusGraphPathVertexSet[i];
        set<Vertex> pathVariantVertexs;
        for (set<Vertex>::iterator variantIter=variantVertexs.begin(); variantIter!=variantVertexs.end(); variantIter++)
        {
            if (pathVertexs.find(*variantIter)!=pathVertexs.end())
            {
                hasVariantVertex = true;
                pathVariantVertexs.insert(*variantIter);
            }
        }

        int totalNumberVariantVertexInPath = 0;
        for (set<Vertex>::iterator vertexIter=pathVertexs.begin(); vertexIter!=pathVertexs.end(); vertexIter++)
        {
            int v = *vertexIter;
            if (consensusGraph.m_isMismatch[v])
            {
                totalNumberVariantVertexInPath += 1;
            }
        }

        if (hasVariantVertex && totalNumberVariantVertexInPath<=pathVariantVertexs.size())
        {
            haplotypes.push_back(topRankConsensusGraphPaths[i]);
            haplotypeToPathIndex.push_back(i);
            haplotypeVariantVertexs.push_back(pathVariantVertexs);

            kk++;
        }
    }

    int numHaplotypes = haplotypes.size();

    // skip if there is no variant haplotype
    if (numHaplotypes==1)
    {
        return;
    }

    // compute haplotype data likelihood
    vector<vector<long double>> haplotypeDataLikelihoods(numHaplotypes);
    PyroHMMsnpHaplotypeDataLikelihood(probAligner, snpCallSettings.m_band, numHaplotypes, haplotypes, readsInWindow, haplotypeDataLikelihoods);


    // genotype
    vector<vector<int>> genotypes;
    set<set<int>> genotypeDiscovered;
    for (int i=0; i<numHaplotypes; i++)
    {
        vector<int> precedeHaplotypes;
        PyroHMMsnpGenotypeSet(snpCallSettings.m_ploidy, i, numHaplotypes, precedeHaplotypes, genotypes, genotypeDiscovered);
    }

    int numGenotypes = genotypes.size();

    // genotype variant vertex
    vector<set<Vertex>> genotypeVariantVertexs;
    for (int i=0; i<numGenotypes; i++)
    {
        set<Vertex> variantVertexInGenotype;
        for (int j=0; j<settingForPyroHMMsnp.m_ploidy; j++)
        {
            int haplotype = genotypes[i][j];
            set<Vertex> variantVertexInHaplotype = haplotypeVariantVertexs[haplotype];
            variantVertexInGenotype.insert(variantVertexInHaplotype.begin(), variantVertexInHaplotype.end());
        }
        genotypeVariantVertexs.push_back(variantVertexInGenotype);
    }

    // genotype priors
    vector<long double> genotypePriors(numGenotypes);
    PyroHMMsnpGenotypePrior(numGenotypes, genotypes, settingForPyroHMMsnp, genotypePriors);

    // genotype likelihoods
    vector<long double> genotypeLikelihoods(numGenotypes);
    PyroHMMsnpGenotypeLikelihood(numGenotypes, genotypes, readsInWindow.size(), haplotypeDataLikelihoods, snpCallSettings, genotypeLikelihoods);

    // genotype posteriors
    vector<long double> genotypePosteriors(numGenotypes);
    PyroHMMsnpGenotypePosterior(numGenotypes, genotypePriors, genotypeLikelihoods, genotypePosteriors);

    // search maximal genotype posterior
    long double maxGenotypePosterior = 0;
    int inferGenotype;
    for (int i=1; i<numGenotypes; i++)
    {
        if (maxGenotypePosterior<genotypePosteriors[i])
        {
            maxGenotypePosterior = genotypePosteriors[i];
            inferGenotype = i;
        }
    }

    // all variant vertexs in the inferred genotype
    set<Vertex> inferGenotypeVariantVertexs = genotypeVariantVertexs[inferGenotype];

    // count haploid type of variant
    map<Vertex,vector<int>> inferGenotypeVariantHaploidType;
    set<Vertex>::iterator inferVariantIter = inferGenotypeVariantVertexs.begin();
    for (; inferVariantIter!=inferGenotypeVariantVertexs.end(); inferVariantIter++)
    {
        int v = *inferVariantIter;
        vector<int> variantHaploidType;
        for (int j=0; j<settingForPyroHMMsnp.m_ploidy; j++)
        {
            int haplotype = genotypes[inferGenotype][j];
            set<Vertex> variantVertexInHaplotype = haplotypeVariantVertexs[haplotype];
            if (variantVertexInHaplotype.find(v)==variantVertexInHaplotype.end())
            {
                variantHaploidType.push_back(0);
            }else
            {
                variantHaploidType.push_back(1);
            }
        }
        inferGenotypeVariantHaploidType[v] = variantHaploidType;
    }
    // variant score
    map<Vertex,long double> inferGenotypeVariantScore;
    inferVariantIter = inferGenotypeVariantVertexs.begin();
    for (; inferVariantIter!=inferGenotypeVariantVertexs.end(); inferVariantIter++)
    {
        int v = *inferVariantIter;
        long double variantScore = 0;
        for (int i=0; i<numGenotypes; i++)
        {
            set<Vertex> variantVertexInGenotype = genotypeVariantVertexs[i];
            if (variantVertexInGenotype.find(v)!=variantVertexInGenotype.end())
                variantScore += genotypePosteriors[i];
        }

        inferGenotypeVariantScore[v] = variantScore;
    }

    // save variant result
    inferVariantIter = inferGenotypeVariantVertexs.begin();
    for (; inferVariantIter!=inferGenotypeVariantVertexs.end(); inferVariantIter++)
    {
        GenericVariant result;

        int v = *inferVariantIter;
        int a = mapVertexToAllele[v];

        int variantChrID;
        int variantChrPos;

        vector<int> haploidType = inferGenotypeVariantHaploidType[v];
        for (int j=0; j<settingForPyroHMMsnp.m_ploidy; j++)
        {
            if (haploidType[j]==0)
            {
                int g = consensusGraph.m_genomePosition[v];

                Allele allele;
                allele.m_allele = consensusGraph.m_labels[g];
                result.m_alleles.push_back(allele);
            }else
            {
                Allele allele = allelePool[a];
                result.m_alleles.push_back(allele);

                variantChrID  = allele.m_chrID;
                variantChrPos = allele.m_chrPosition;
            }
        }

        result.m_chrID           = variantChrID;
        result.m_chrPosition     = variantChrPos;
        result.m_probScoreRef    = genotypePosteriors[0];
        result.m_probScoreVar    = genotypePosteriors[inferGenotype];
        result.m_variantType     = VARIANT_SNP;
        long double variantScore = inferGenotypeVariantScore[v];
        if (fabs(1-variantScore)<1e-300)
            result.m_quality     = 3000;
        else if (variantScore<1e-300)
            result.m_quality     = 0;
        else
            result.m_quality     = -10*log10(1-variantScore);

        char refBase;
        fastaObj.GetBase(result.m_chrID, result.m_chrPosition, refBase);
        result.m_reference       = refBase;

        for (int i=0; i<result.m_alleles.size(); i++)
        {
            if (result.m_alleles[i].m_allele==result.m_reference)
                result.m_haploidType.push_back(0);
            else
                result.m_haploidType.push_back(1);
        }


        // filter
        if (result.m_quality>=snpCallSettings.m_variantQualityFilter)
            variantResults.push_back(result);

    }

}


