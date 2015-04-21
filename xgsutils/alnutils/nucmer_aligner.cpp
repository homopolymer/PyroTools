#include <regex>
#include <ctime>
#include <iostream>
#include <string>
#include <tuple>
using namespace std;

#include "../c++/kseq.h"
#include "../c++/utils.h"
using namespace XgsUtils;

#include <omp.h>
#include <zlib.h>
#include <getopt.h>

// help message
int help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    nucmer_aligner [OPTIONS] <fas1> <fas2>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    Invoke nucmer to compute the alignments of all pairs of sequences in fas1 and fas2" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "    -a,--all        compute all alignments" << endl;
    cerr << "    -f,--format     specify the output format, valid argument: bed and sam (default:bed) [STR]" << endl;
    cerr << "    -c,--cores      specify the number of threads (default:1) [INT]" << endl;
    cerr << "    -v,--verbose    specify the verbose level" << endl;
    cerr << "    -h,--help       print help message" << endl;
    return 0;
}

// parse the command line
int parseCommandLine(int argc, char *argv[], string &fas1, string &fas2, string &fmt, bool& all, int &cores, int &verbose){
    fmt = "bed";
    cores = 1;
    verbose = 0;
    all = false;

    int c;
    while(1){
        static struct option long_options[] =
        {
            {"all",     no_argument,       0, 'a'},
            {"format",  required_argument, 0, 'f'},
            {"cores",   required_argument, 0, 'c'},
            {"verbose", required_argument, 0, 'v'},
            {"help",    no_argument,       0, 'h'},
            {0,0,0,0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "af:v:h", long_options, &option_index);

        // detect the end of the options
        if (c==-1) break;

        switch(c){
            case 'h':
                help();
                exit(0);
                break;
            case 'a':
                all = true;
                break;
            case 'f':
                fmt = optarg;
                break;
            case 'c':
                cores = stoi(optarg);
                break;
            case 'v':
                verbose = stoi(optarg);
                break;
            default:
                abort();
        }
    }

    fas1 = "";
    fas2 = "";
    // file for fas1
    if (optind < argc)
        fas1 = DirectoryUtils::AbsolutePath(argv[optind++]);
    // file for fas2
    if (optind < argc)
        fas2 = DirectoryUtils::AbsolutePath(argv[optind++]);

    if (fas1.empty() || fas2.empty()){
        VerboseUtils::Verbose("nucmer_aligner", "required input(s) not provided");
        VerboseUtils::Verbose("nucmer_aligner", "abort");
        exit(-1);
    }
    return 0;
}

// parse the nucmer output
void parseNucmerResult(string &aln, string &a, string &b, string &cigar, uint32_t &p0, uint32_t &p1, uint32_t &q0, uint32_t &q1, int &nm, int &nx, int &ni, int &nd, double &iden){
    a = "";
    b = "";
       
    // regex pattern
    regex alninfo("^-- BEGIN alignment \\[ [+-]1 (\\d+) - (\\d+) \\| [+-]1 (\\d+) - (\\d+) \\]");
    regex alnaln("^\\d+\\s+([AaCcGgTtNn\\.]+)");

    // stringstream 
    stringstream ss(aln);
    string line;
    smatch matches;
    ssub_match ssm;
    int i = 0;
    while(getline(ss,line)){
        if (regex_search(line, matches, alninfo)){
            ssm = matches[1];
            p0 = stoi(ssm.str())-1;
            ssm = matches[2];
            p1 = stoi(ssm.str())-1;
            ssm = matches[3];
            q0 = stoi(ssm.str())-1;
            ssm = matches[4];
            q1 = stoi(ssm.str())-1;           
        }
        if (regex_search(line, matches, alnaln)){
            if (i%2==0){
                ssm = matches[1];
                a += ssm.str();
            }
            if (i%2==1){
                ssm = matches[1];
                b += ssm.str();
            }
            i++;
        }
    }

    // compute cigar and nm and iden
    nm = 0;
    nx = 0;
    ni = 0;
    nd = 0;
    string status = "";
    for (int j=0; j<a.length(); j++){
        if (a[j]==b[j]) { status += "M"; }
        else if (a[j]=='.') { a[j] = '-'; status += "I"; nm += 1; ni += 1; }
        else if (b[j]=='.') { b[j] = '-'; status += "D"; nm += 1; nd += 1; }
        else { status += "X"; nm += 1; nx += 1; }
    }
    iden = 1-nm/(a.length()+0.0);

    vector<tuple<char,int>> cigar_array;
    for (auto s : status){
        if (s=='X') s='M';
        if (cigar_array.empty()){
            cigar_array.emplace_back(tuple<char,int>(s,1));
        }else{
            if (get<0>(*cigar_array.rbegin())==s){
                get<1>(*cigar_array.rbegin()) += 1;
            }else{
                cigar_array.emplace_back(tuple<char,int>(s,1));
            }
        }
    }
    cigar = "";
    for (auto c : cigar_array){
        cigar += to_string(get<1>(c))+get<0>(c);
    }
    
}

// convert nucmer output to bed format
void writeBedOutput(string &nucmer_output, string &seq1name, string &seq2name, string &seq1seq, string &seq2seq, uint32_t l1, uint32_t l2){
    string a,b;
    string cigar;
    uint32_t p0,p1,q0,q1;
    int nm,nx,ni,nd;
    double iden;
 
    parseNucmerResult(nucmer_output, a, b, cigar, p0, p1, q0, q1, nm, nx, ni, nd, iden);

    stringstream ss;
    if (iden == iden)
        ss << seq1name << "\t" << p0 << "\t" << p1 << "\t" << seq2name << "\t" << iden << "\t" << nm << "\t"  << nx << "\t" << ni << "\t" << nd << "\t" << a << "\t" << b << endl;
    else
        ss << seq1name << "\t" << "*" << "\t" << "*" << "\t" << seq2name << "\t" << 0 << "\t" << nm << "\t"  << nx << "\t" << ni << "\t" << nd << "\t" << a << "\t" << b << endl;
    cout << ss.str();
}

// convert nucmer output to sam format
void writeSamOutput(string &nucmer_output, string &seq1name, string &seq2name, string &seq1seq, string &seq2seq, uint32_t l1, uint32_t l2){
    string a,b;
    string cigar;
    uint32_t p0,p1,q0,q1;
    int nm,nx,ni,nd;
    double iden;

    parseNucmerResult(nucmer_output, a, b, cigar, p0, p1, q0, q1, nm, nx, ni, nd, iden);

    // add soft clip site
    if (q0>0) cigar = to_string(q0) + "S" + cigar;
    if (q1+1<l2) cigar = cigar + to_string(l2-q1-1) + "S";

    if (iden == iden){
        stringstream ss;
        ss << seq2name << "\t" << 0 << "\t" << seq1name << "\t" << p0+1 << "\t" << 254 << "\t" << cigar << "\t" << "*\t0\t0" << "\t" << seq2seq << "\t" << "*" << "\t" << "NM:i:" << nm << endl;
        cout << ss.str();
    }
}

KSEQ_INIT(gzFile,gzread)

int main(int argc, char *argv[])
{
    // print help message if no argument provided
    if (argc==1){
        help();
        exit(0);
    }
    
    string fas1, fas2;
    string fmt;
    bool all;
    int vbos;
    int cores;
    // parse the command line
    parseCommandLine(argc, argv, fas1, fas2, fmt, all, cores, vbos);

    if (fmt == "sam"){
        cout << "@HD\tVN:1.4\tSO:queryname" << endl;
    }

    // sequences
    int n1=0, n2=0;
    typedef tuple<string,string,uint32_t> SEQ;
    vector<SEQ> fas1seqs, fas2seqs;
    gzFile fp;
    kseq_t *seq;
    int l;
    // sequences in fas1
    fp = gzopen(fas1.c_str(),"r");
    seq = kseq_init(fp);
    while((l=kseq_read(seq))>=0){
        if (fmt == "sam") cout << "@SQ" << "\t" << "SN:" << seq->name.s << "\t" << "LN:" << (int32_t)seq->seq.l << endl;
        fas1seqs.emplace_back(SEQ(seq->name.s,seq->seq.s,strlen(seq->seq.s)));
        n1++;
    }
    gzclose(fp);
    kseq_destroy(seq);
    // sequences in fas2
    fp = gzopen(fas2.c_str(),"r");
    seq = kseq_init(fp);
    while((l=kseq_read(seq))>=0){
        fas2seqs.emplace_back(SEQ(seq->name.s,seq->seq.s,strlen(seq->seq.s)));
        n2++;
    }
    gzclose(fp);
    kseq_destroy(seq);

    // get the current working directory
    string cwd = DirectoryUtils::CurrentWorkDir();
    // make a temporary directory
    string twd = cwd + "/temporary" + to_string(clock());
    DirectoryUtils::MkDir(twd);
    // change to the temporary directory
    DirectoryUtils::ChDir(twd);

    stringstream ss;
    string cmd;
    string maxmatch="";
    if (all) maxmatch="--maxmatch";
    // command to invoke nucmer
    ss << "nucmer " << maxmatch << " " << fas1 << " " << fas2 << " >/dev/null 2>&1";
    cmd = ss.str();
    // debug
    cerr << cmd << endl;
    // run the command
    SystemUtils::Run(cmd);
    // clean the stringstream
    ss.clear();
    ss.str("");
    
    // iteratively parse the alignments
    omp_set_dynamic(0);
    omp_set_num_threads(cores);
    #pragma omp parallel for
    for (int i=0; i<n1; i++){
        string s1 = get<0>(fas1seqs[i]);
        string s1s = get<1>(fas1seqs[i]);
        uint32_t l1 = get<2>(fas1seqs[i]);
        for (int j=0; j<n2; j++){
            string s2 = get<0>(fas2seqs[j]);
            string s2s = get<1>(fas2seqs[j]);
            uint32_t l2 = get<2>(fas2seqs[j]);
            ss << "show-aligns out.delta " << s1 << " " << s2 << " 2>/dev/null";
            cmd = ss.str();
            string result;
            SystemUtils::Run(cmd,result);
            if (fmt == "bed") writeBedOutput(result,s1,s2,s1s,s2s,l1,l2);
            if (fmt == "sam") writeSamOutput(result,s1,s2,s1s,s2s,l1,l2);
            ss.clear();
            ss.str("");
        }
    }
    
    DirectoryUtils::ChDir(cwd);
    DirectoryUtils::RmDir(twd);

    return 0;
}
