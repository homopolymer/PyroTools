#include <string>
#include <list>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <tuple>
using namespace std;

#include <limits.h>
#include <zlib.h>
#include <omp.h>
#include <getopt.h>
#include <libgen.h>
#include <time.h>

#include "kseq.h"
#include "ssw_cpp.h"

// verbose function
int Verbose(string msg){
    time_t now = time(0);
    char buf[128];
    strftime(buf, sizeof(buf), "%c",localtime(&now));
    stringstream ss;
    ss << "[" << buf << "] " << msg << "\n";
    cerr << ss.str();
    return 0;
}

// help message
int help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    ssw_aligner [OPTIONS] <fastafile1> <fastafile2>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    Compute the pairwise alignments of two sets of sequences" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "    -c,--cores       specify the number of threads (default:1) [INT]" << endl;
    cerr << "    -m,--match       specify the match score (default:2) [INT]" << endl;
    cerr << "    -x,--mismatch    specify the mismatch penalty (default:2) [INT]" << endl;
    cerr << "    -o,--gap-open    specify the gap open penalty (default:3) [INT]" << endl;
    cerr << "    -e,--gap-ext     specify the gap extend penalty (default:1) [INT]" << endl;
    cerr << "    -f,--outfmt      specify the output format, valid arguments: bed and sam (default:bed)" << endl;
    cerr << "    -h,--help        print help message" << endl;
}

int parseCommandLine(int argc, char *argv[], int &cores, int &match, int &mismatch, int &gapopen, int &gapext, string &outfmt, vector<string> &files)
{
    cores=1;
    match=2;
    mismatch=2;
    gapopen=3;
    gapext=1;
    outfmt="bed";

    int c;
    while(1){
        static struct option long_options[] =
        {
            {"match",   required_argument,0,'m'},
            {"mismatch",required_argument,0,'x'},
            {"gap-open",required_argument,0,'o'},
            {"gap-ext", required_argument,0,'e'},
            {"outfmt",  required_argument,0,'f'},
            {"cores",   required_argument,0,'c'},
            {"help",    no_argument,      0,'h'},
            {0,0,0,0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "c:m:x:o:e:f:h", long_options, &option_index);

        // detect the end of the options
        if (c==-1) break;

        switch(c){
            case 'm':
                match = stoi(optarg);
                break;
            case 'x':
                mismatch = stoi(optarg);
                break;
            case 'o':
                gapopen = stoi(optarg);
                break;
            case 'e':
                gapopen = stoi(optarg);
                break;
            case 'h':
                help();
                exit(0);
                break;
            case 'c':
                cores = stoi(optarg);
                break;
            case 'f':
                outfmt=optarg;
                break;
            default:
                Verbose("invalid option...abort");
                abort();
        }
    }


   for (; optind<argc; optind++){
       char result[PATH_MAX];
       realpath(argv[optind],result);
       files.emplace_back(result);
   }
}

void progressbar(unsigned long long int x, unsigned long long int n, float t, int w=30){
    stringstream ss;

//    if ( (x != n) && (x % (n/100+1) != 0) ) return;
    double deta = (n-x)*t;
    double rest = 0;
    string eta;
    if (1){
        stringstream ss;
        unsigned long long int yrs,days,hrs,minus,secs;
        yrs = (unsigned long long int)(deta/365./24./60./60.);
        rest = deta-yrs*365*24*60*60;
        days = (unsigned long long int)(rest/24./60./60.);
        rest = rest-days*24*60*60;
        hrs = (unsigned long long int)(rest/60./60.);
        rest = rest-hrs*60*60;
        minus = (unsigned long long int)(rest/60.);
        rest = rest-minus*60;
        secs = rest;
        ss << yrs << "y" << days << "d" << hrs << "h" << minus << "m" << fixed << setprecision(2) << secs << "s";
        eta = ss.str();
    }

//    float ratio = x/(float)n;
//    unsigned long long int c = ratio * w;

//    ss << setw(3) << (int)(ratio*100) << "% [";
//    for (unsigned long long int y=0; y<c; y++) ss << "=";
//    if (c<w) ss << ">";
//    for (unsigned long long int y=c+1; y<w; y++) ss << " ";
//    ss << "]    ";
    ss << "\r" << flush;
    ss << "alns/sec: " << setw(5) << fixed << setprecision(2) << 1./t
       << "    processed: " << scientific << setw(7) << x
       << "    eta: " << eta << "  \r" << flush;

    cerr << ss.str();
}


int align(string &seqa, string &seqb, int match, int mismatch, int gapopen, int gapext, string &alnseqa, string &alnseqb, int32_t &a0, int32_t &a1, int32_t &b0, int32_t &b1, uint16_t &score1, uint16_t &score2, int32_t &nx, int32_t &ni, int32_t &nd, string &cigar){

    StripedSmithWaterman::Aligner aligner(match, mismatch, gapopen, gapext);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

    aligner.Align(seqb.c_str(), seqa.c_str(), seqa.size(), filter, &alignment);

    // retrieve the aligned sequences
    alnseqa = "";
    alnseqb = "";
    int i = alignment.ref_begin;
    int j = alignment.query_begin;
    for (auto c : alignment.cigar){
        int op = c&0xf;
        int ol = c>>4;
        if (op==0 || op==7){ // match
            for (int l=0; l<ol; l++){
                alnseqa += seqa[i++];
                alnseqb += seqb[j++];
            }
        }else if (op==8){ // mismatch
            for (int l=0; l<ol; l++){
                alnseqa += seqa[i++];
                alnseqb += seqb[j++];
            }
        }else if (op==2){ // delete
            for (int l=0; l<ol; l++){
                alnseqa += seqa[i++];
                alnseqb += "-";
            }
        }else if (op==1){ // insert
            for (int l=0; l<ol; l++){
                alnseqa += "-";
                alnseqb += seqb[j++];
            }
        }
    }

    // compute align length and difference
    nx = 0;
    ni = 0;
    nd = 0;
    for (auto c : alignment.cigar){
        int op=c&0xf;
        int ol=c>>4;
        if (op==8){
            nx += ol;
        }
        else if (op==2){
            nd += ol;
        }
        else if (op==1){
            ni += ol;
        }
    }

    // begin and end position
    a0 = alignment.ref_begin+1;
    a1 = alignment.ref_end+1;
    b0 = alignment.query_begin+1;
    b1 = alignment.query_end+1;

    // cigar
    cigar = alignment.cigar_string;
    if (alignment.query_begin>0)
        cigar = to_string(alignment.query_begin)+"S"+cigar;
    if (alignment.query_end+1<seqb.length())
        cigar = cigar + to_string(seqb.length()-alignment.query_end-1) + "S";

    // alignment score
    score1 = alignment.sw_score;
    score2 = alignment.sw_score_next_best;

    return (nx+ni+nd);

}

KSEQ_INIT(gzFile,gzread)

int main(int argc, char *argv[]){
    // print help message if no arguments provided
    if (argc==1){
        help();
        exit(0);
    }
 
    // parse the command line
    int cores;
    int match,mismatch,gapopen,gapext;
    string outfmt;
    vector<string> fastafiles;
    parseCommandLine(argc,argv,cores,match,mismatch,gapopen,gapext,outfmt,fastafiles);


    // read all sequences in memory
    Verbose("load all sequences in memory");
    vector<int> fastaseqnums;
    list<tuple<string,string>> fastaseqs;
    for (auto fastafile : fastafiles){
        int n = 0;
        int l;
        gzFile fp = gzopen(fastafile.c_str(),"r");
        kseq_t *seq = kseq_init(fp);
        while ((l=kseq_read(seq))>=0){
            fastaseqs.emplace_back(make_tuple(seq->name.s,seq->seq.s));
            n++;
        }
        fastaseqnums.emplace_back(n);
        kseq_destroy(seq);
        gzclose(fp);
    }

    // output sam header
    if (outfmt == "sam"){
        cout << "@HD\tVN:1.4\tSO:queryname" << endl;
        gzFile fp = gzopen(fastafiles[0].c_str(),"r");
        kseq_t *seq = kseq_init(fp);
        int l;
        while ((l=kseq_read(seq))>=0){
            cout << "@SQ" << "\t" << "SN:" << seq->name.s << "\t" << "LN:" << (int32_t)seq->seq.l << endl;
        }
        kseq_destroy(seq);
        gzclose(fp);
    }
    
    Verbose("compute the pairwise alignments");

    int num1 = fastaseqnums[0];
    int num2 = fastaseqnums[1];
    int numpairs = 0;
    unsigned long long int x = 1;
    unsigned long long int n = num1*num2;
    float tt = 0;
    // compute the pairwise alignments
    omp_set_dynamic(0);
    omp_set_num_threads(cores);
    #pragma omp parallel for shared(numpairs,x,n,tt)
    for (int i=0; i<num1; i++){
        for (int j=0; j<num2; j++){
            clock_t t0 = clock();
            auto pi = fastaseqs.begin();
            auto pj = fastaseqs.begin();
            advance(pi,i);
            advance(pj,j+num1);
            string seqa_name = get<0>(*pi);
            string seqa_seq = get<1>(*pi);
            string seqb_name = get<0>(*pj);
            string seqb_seq = get<1>(*pj);

            int32_t a0,a1,b0,b1;
            int32_t nx,ni,nd;
            uint16_t score1,score2;
            string alna,alnb,cigar;
            int nm = align(seqa_seq, seqb_seq, match, mismatch, gapopen, gapext, alna, alnb, a0, a1, b0, b1, score1, score2, nx, ni, nd, cigar);

            uint32_t mapq = -4.343 * log(1 - (double)abs(score1 - score2)/(double)score1);
            mapq = (uint32_t) (mapq + 4.99);
            mapq = mapq < 254 ? mapq : 254;

            stringstream ss;
            if (outfmt=="bed") ss << seqa_name << "\t" << a0 << "\t" << a1 << "\t" << seqb_name << "\t" << score1 << "\t" << nm << "\t" << nx << "\t" << ni << "\t" << nd << "\t" << alna << "\t" << alnb << endl;
            if (outfmt=="sam") ss << seqb_name << "\t" << 0 << "\t" << seqa_name << "\t" << a0 << "\t" << mapq << "\t" << cigar << "\t" << "*" << "\t" << 0 << "\t" << 0 << "\t" << seqb_seq << "\t" << "*" << "\t" << "NM:i:" << nx << endl;
            cout << ss.str();

            clock_t t1 = clock();

            tt += float(t1-t0)/CLOCKS_PER_SEC;
            float at = tt/x;
            progressbar(x++,n,at);

            numpairs++;
        }
    }
    cerr << endl;

    Verbose("totally process "+to_string(numpairs)+" alignments");

    return 0;
}
