#include "GenericVcfTools.h"
using namespace GenericSequenceTools;

#include <set>
#include <fstream>
#include <iostream>
using namespace std;

GenericVcfTools::GenericVcfTools()
{
}


void writeVcfHeader(ostream& out, VariantCallSetting& variantCallingSettings)
{
    out << "##fileformat=VCFv4.1" << endl;
    out << "##source=PyroTools" << endl;
    out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">" << endl;
    out << "##INFO=<ID=GDF,Number=1,Type=Integer,Description=\"Number of reads on forward strand\">" << endl;
    out << "##INFO=<ID=GDR,Number=1,Type=Integer,Description=\"Number of reads on reverse strand\">" << endl;
    out << "##INFO=<ID=ADF,Number=1,Type=Integer,Description=\"Number of variant-supporting reads on forward strand\">" << endl;
    out << "##INFO=<ID=ADR,Number=1,Type=Integer,Description=\"Number of variant-supporting reads on reverse strand\">" << endl;
    out << "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root-mean-square mapping quality of covering reads\">" << endl;
    out << "##INFO=<ID=MQ1,Number=1,Type=Integer,Description=\"Root-mean-square mapping quality of reads having ALT allele\">" << endl;
    out << "##FORMAT=<ID=GT,Number=1,Type=String, Description=\"Genotype\">" << endl;
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    if (!variantCallingSettings.m_sampleList.empty())
    {
        set<string>::iterator SampleListBegin = variantCallingSettings.m_sampleList.begin();
        set<string>::iterator SampleListEnd   = variantCallingSettings.m_sampleList.end();
        set<string>::iterator sampleIter;
        for (sampleIter=SampleListBegin; sampleIter!=SampleListEnd; sampleIter++)
            out << "\t" << *sampleIter;
    }else
    {
        out << "\tSample1";
    }
    out << endl;
}

void writeVcfRecordInfoDp(ostream& out, GenericVariant& v)
{
    int dp;
    for (int i=0; i<v.m_alleles.size(); i++)
    {
        Allele a = v.m_alleles[i];
        if (a.m_allele==v.m_reference)
            continue;
        dp = a.m_globalDepth;
    }
    out << "DP=" << dp;
}

void writeVcfRecordInfoStrand(ostream& out, GenericVariant& v)
{
    if (v.m_variantType==VARIANT_SNP)
    {
        for (int i=0; i<v.m_alleles.size(); i++)
        {
            Allele a = v.m_alleles[i];
            if (a.m_allele==v.m_reference)
                continue;
            out << "GDF=" << a.m_globalStrandPos << ";";
            out << "GDR=" << a.m_globalStrandNeg;
            break;
        }
    }
}

void writeVcfRecordInfoVariantStrand(ostream& out, GenericVariant& v)
{
    if (v.m_variantType==VARIANT_SNP)
    {
        for (int i=0; i<v.m_alleles.size(); i++)
        {
            Allele a = v.m_alleles[i];
            if (a.m_allele==v.m_reference)
                continue;
            out << "ADF=" << a.m_alleleStrandPos << ";";
            out << "ADR=" << a.m_alleleStrandNeg;
            break;
        }
    }
}

void writeVcfRecordInfoMq(ostream& out, GenericVariant& v)
{
    for (int i=0; i<v.m_alleles.size(); i++)
    {
        Allele a = v.m_alleles[i];
        if (a.m_allele==v.m_reference)
            continue;
        out << "MQ=" << a.m_globalMapAvgQual;
        break;
    }
}

void writeVcfRecordInfoMq1(ostream& out, GenericVariant& v)
{
    if (v.m_variantType==VARIANT_SNP)
    {
        for (int i=0; i<v.m_alleles.size(); i++)
        {
            Allele a = v.m_alleles[i];
            if (a.m_allele==v.m_reference)
                continue;
            out << "MQ1=" << a.m_alleleMapAvgQual;
        }
    }
}

void writeVcfRecordInfo(ostream& out, GenericVariant& v)
{
    // DP
    writeVcfRecordInfoDp(out, v);
    out << ";";

    // Global Read Strand
    writeVcfRecordInfoStrand(out, v);
    out << ";";

    // MQ
    writeVcfRecordInfoMq(out, v);
    out << ";";

    // MQ1
    writeVcfRecordInfoMq1(out, v);

    out << "\t";
}

void writeVcfRecordFormat(ostream& out, GenericVariant& v)
{
    out << "GT" << "\t";
    for (int i=0; i<v.m_haploidType.size(); i++)
    {
        out << v.m_haploidType[i];
        if (i!=v.m_haploidType.size()-1)
            out << "/";
    }
}

void writeVcfRecord(ostream& out, RefVector& chromosomes, GenericVariant& v)
{
    out << chromosomes[v.m_chrID].RefName << "\t";
    out << v.m_chrPosition+1 << "\t";
    out << "." << "\t";
    out << v.m_reference << "\t";

    set<string> hasDiscovered;
    vector<string> uniqAlleles;
    for (int i=0; i<v.m_alleles.size(); i++)
    {
        Allele a = v.m_alleles[i];

        if (a.m_allele==v.m_reference)
            continue;

        if (hasDiscovered.find(a.m_allele)!=hasDiscovered.end())
            continue;

        uniqAlleles.push_back(a.m_allele);

        hasDiscovered.insert(a.m_allele);
    }
    for (int i=0; i<uniqAlleles.size(); i++)
    {
        string a = uniqAlleles[i];

        out << a;

        if (i!=uniqAlleles.size()-1)
            out << ",";
        else
            out << "\t";
    }

    out << v.m_quality << "\t";
    out << "." << "\t";

    // INFO
    writeVcfRecordInfo(out, v);

    // Format
    writeVcfRecordFormat(out, v);

    out << endl;
}

void GenericVcfTools::write(RefVector& chromosomes, string &vcfFileName, VariantCallSetting& variantCallSettings, vector<GenericVariant> &variantsToReport)
{
    ofstream out;
    out.open(vcfFileName);

    writeVcfHeader(out, variantCallSettings);

    for (int i=0; i<variantsToReport.size(); i++)
    {
        writeVcfRecord(out, chromosomes, variantsToReport[i]);
    }

    out.close();
}
