#include <string>
#include <forward_list>
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
    cerr << "    seq_ssw_sim [OPTIONS] <fastafile1> [<fastafile2>]" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    Compute the similarity by using the striped smith-waterman alignment" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "    -c,--cores    specify the number of threads (default:1) [INT]" << endl;
    cerr << "    --match       specify the match score (default:2) [INT]" << endl;
    cerr << "    --mismatch    specify the mismatch penalty (default:2) [INT]" << endl;
    cerr << "    --gap-open    specify the gap open penalty (default:3) [INT]" << endl;
    cerr << "    --gap-ext     specify the gap extend penalty (default:1) [INT]" << endl;
    cerr << "    -p,--pos      report the target begin and end positions" << endl;
    cerr << "    -b,--both     compute the similarity by considering both positive and negative direction" << endl;
    cerr << "    -h,--help     print help message" << endl;
}

int parseCommandLine(int argc, char *argv[], int &cores, int &match, int &mismatch, int &gapopen, int &gapext, bool &repo_pos, bool &both_strand, vector<string> &files)
{
    cores=1;
    match=2;
    mismatch=2;
    gapopen=3;
    gapext=1;
    repo_pos=false;
    both_strand=false;

    int c;
    while(1){
        static struct option long_options[] =
        {
            {"match",required_argument,0,0},
            {"mismatch",required_argument,0,0},
            {"gap-open",required_argument,0,0},
            {"gap-ext",required_argument,0,0},
            {"pos",no_argument,0,'p'},
            {"both",no_argument,0,'b'},
            {"cores",required_argument,0,'c'},
            {"help",no_argument,0,'h'},
            {0,0,0,0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "c:pbh", long_options, &option_index);

        // detect the end of the options
        if (c==-1) break;

        switch(c){
            case 0:
                switch(option_index){
                    case 0:
                        match = stoi(optarg);
                        break;
                    case 1:
                        mismatch = stoi(optarg);
                        break;
                    case 2:
                        gapopen = stoi(optarg);
                        break;
                    case 3:
                        gapopen = stoi(optarg);
                        break;
                }
            case 'h':
                help();
                exit(0);
                break;
            case 'c':
                cores = stoi(optarg);
                break;
            case 'p':
                repo_pos=true;
                break;
            case 'b':
                both_strand=true;
                break;
            default:
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

float seqsim(string &seqa, string &seqb, int match, int mismatch, int gapopen, int gapext, int &a0, int &a1, int &b0, int &b1){

    StripedSmithWaterman::Aligner aligner(match, mismatch, gapopen, gapext);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

    aligner.Align(seqb.c_str(), seqa.c_str(), seqa.size(), filter, &alignment);

    // compute align length and differenct
    float alnlen = 0;
    float alndiff = 0;
    for (auto c : alignment.cigar){
        int op=c&0xf;
        int ol=c>>4;
        if (op==0 || op==7) {
            alnlen += ol;
        }
        else if (op==8){
            alnlen += ol;
            alndiff += ol;
        }
        else if (op==2){
            alnlen += ol;
            alndiff += ol;
        }
        else if (op==1){
            alnlen += ol;
            alndiff += ol;
        }
    }

    a0 = alignment.ref_begin+1;
    a1 = alignment.ref_end+1;
    b0 = alignment.query_begin+1;
    b1 = alignment.query_end+1;

    return (1-(alndiff/alnlen));

}

void DnaReverseComplement(string &dna,string &dna_rev){
    dna_rev="";
    for (auto ptr=dna.rbegin(); ptr<dna.rend(); ptr++){
        if (*ptr=='A' || *ptr=='a')  dna_rev+='T';
        else if (*ptr=='C' || *ptr=='c')  dna_rev+='G';
        else if (*ptr=='G' || *ptr=='g')  dna_rev+='C';
        else if (*ptr=='T' || *ptr=='t')  dna_rev+='A';
        else dna_rev+=*ptr;
    }
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
    bool repo_pos, both_strand;
    vector<string> fastafiles;
    parseCommandLine(argc,argv,cores,match,mismatch,gapopen,gapext,repo_pos,both_strand,fastafiles);


    // read all sequences in memory
    Verbose("load all sequences in memory");
    vector<int> fastaseqnums;
    forward_list<tuple<string,string>> fastaseqs;
    for (auto fastafile : fastafiles){
        int n = 0;
        int l;
        gzFile fp = gzopen(fastafile.c_str(),"r");
        kseq_t *seq = kseq_init(fp);
        while ((l=kseq_read(seq))>=0){
            fastaseqs.emplace_front(make_tuple(seq->name.s,seq->seq.s));
            n++;
        }
        fastaseqnums.emplace_back(n);
        kseq_destroy(seq);
        gzclose(fp);
    }
    
    Verbose("compute the pairwise alignments");

    int num1 = fastaseqnums[0];
    int numpairs = 0;
    unsigned long long int x = 1;
    unsigned long long int n = num1*(num1-1)/2;
    if (fastaseqnums.size()==2)
        n = num1*fastaseqnums[1];
    float tt = 0;
    // compute the similarity
    omp_set_dynamic(0);
    omp_set_num_threads(cores);
    #pragma omp parallel for shared(numpairs,x,n,tt)
    for (int i=0; i<num1; i++){
        // mode 1
        if (fastaseqnums.size()==1){
            for (int j=i+1; j<num1; j++){
                clock_t t0 = clock();

                auto pi = fastaseqs.begin();
                auto pj = fastaseqs.begin();
                advance(pi,i);
                advance(pj,j);
                string seqa_name = get<0>(*pi);
                string seqa_seq = get<1>(*pi);
                string seqb_name = get<0>(*pj);
                string seqb_seq = get<1>(*pj);


                int a0,a1,b0,b1;
                // forward direction
                float sim = seqsim(seqa_seq, seqb_seq, match, mismatch, gapopen, gapext, a0, a1, b0, b1);
                // reverse direction
                if (both_strand){
                    string seqb_seq_rev;
                    DnaReverseComplement(seqb_seq, seqb_seq_rev);
                    float sim2 = seqsim(seqa_seq, seqb_seq_rev, match, mismatch, gapopen, gapext, a0, a1, b0, b1);
                    if (sim<sim2) sim = sim2;
                }
                stringstream ss;
                ss << seqa_name << "\t" << seqb_name << "\t" << fixed << setprecision(3) << setw(4) << sim;
                if (repo_pos) ss << "\t" << a0 << "\t" << a1 << "\t" << b0 << "\t" << b1;
                ss << "\n";
                cout << ss.str();

                clock_t t1 = clock();

                tt += float(t1-t0)/CLOCKS_PER_SEC;
                float at = tt/x;
                progressbar(x++,n,at);

                numpairs++;
            }
        }

        // mode 2
        if (fastaseqnums.size()==2){
            int num2 = fastaseqnums[1];
            for (int j=0; j<num2; j++){
                
                clock_t t0 = clock();
                auto pi = fastaseqs.begin();
                auto pj = fastaseqs.begin();
                advance(pi,i+num2);
                advance(pj,j);
                string seqa_name = get<0>(*pi);
                string seqa_seq = get<1>(*pi);
                string seqb_name = get<0>(*pj);
                string seqb_seq = get<1>(*pj);


                int a0,a1,b0,b1;
                // forward direction
                float sim = seqsim(seqa_seq, seqb_seq, match, mismatch, gapopen, gapext, a0, a1, b0, b1);
                // reverse direction
                if (both_strand){
                    string seqb_seq_rev;
                    DnaReverseComplement(seqb_seq, seqb_seq_rev);
                    float sim2 = seqsim(seqa_seq, seqb_seq_rev, match, mismatch, gapopen, gapext, a0, a1, b0, b1);
                    if (sim<sim2) sim = sim2;
                }
                stringstream ss;
                ss << seqa_name << "\t" << seqb_name << "\t" << fixed << setprecision(3) << setw(4) << sim;
                if (repo_pos) ss << "\t" << a0 << "\t" << a1 << "\t" << b0 << "\t" << b1;
                ss << "\n";
                cout << ss.str();

                clock_t t1 = clock();

                tt += float(t1-t0)/CLOCKS_PER_SEC;
                float at = tt/x;
                progressbar(x++,n,at);

                numpairs++;
            }
        }
    }
    cerr << endl;

    Verbose("totally process "+to_string(numpairs)+" alignments");

    return 0;
}
