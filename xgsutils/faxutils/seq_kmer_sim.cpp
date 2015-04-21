#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <tuple>
#include <vector>
#include <iterator>
using namespace std;

#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <sys/stat.h>
#include <libgen.h>
#include <omp.h>
#include <getopt.h>

#ifdef __APPLE__
    #include <mach-o/dyld.h>
#endif

// help message
int help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    seq_kmer_sim <fastafile1> [<fastafile2>]" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    compute the sequence similarity based on the kmer spectrums" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "    -k,--ksize    specify the size of the kmer" << endl;
    cerr << "    -o,--osize    specify the overlap size of two nearby kmers" << endl;
    cerr << "    -c,--cores    specify the number of threads" << endl;
    cerr << "    -h,--help     print help message" << endl;
    return 0;
}

int parseCommandLine(int argc, char *argv[], int &ksize, int &osize, int &cores, vector<string> &files)
{
    ksize=31;
    osize=30;
    cores=1;

    int c;
    while(1){
        static struct option long_options[] =
        {
            {"ksize",required_argument,0,'k'},
            {"osize",required_argument,0,'o'},
            {"cores",required_argument,0,'c'},
            {"help",no_argument,0,'h'},
            {0,0,0,0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "k:o:c:h", long_options, &option_index);

        // detect the end of the options
        if (c==-1) break;

        switch(c){
            case 'h':
                help();
                exit(0);
                break;
            case 'k':
                ksize = stoi(optarg);
                break;
            case 'o':
                osize = stoi(optarg);
                break;
            case 'c':
                cores = stoi(optarg);
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

string getWorkDirectory(){
    char result[PATH_MAX];
    getcwd(result,PATH_MAX);
    return string(result);
}

// get the path direct to xgsutils
// reference: http://stackoverflow.com/questions/143174/how-do-i-get-the-directory-that-a-program-is-running-from
string getexepath()
{
    string path;
    char result[PATH_MAX];
    uint32_t size=sizeof(result);
    #ifdef __linux
        ssize_t count = readlink( "/proc/self/exe", result, size );
        path = string( result, (count > 0) ? count : 0 );
    #elif __APPLE__
        _NSGetExecutablePath( result, &size );
        path = string(result);
    #endif
    return string(dirname((char*)path.c_str()));
}

// run the command and get the output
// http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c
string exec(const char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
    	if(fgets(buffer, 128, pipe) != NULL)
    		result += buffer;
    }
    pclose(pipe);
    return result;
}

int run_mode1(string file,int ksize,int osize,int cores){
    string cmd;

    // get the current working directory
    string cwd = getWorkDirectory();

    // make a temporary directory
    string twd = cwd + "/temporary_" + to_string(clock());
    mkdir(twd.c_str(),S_IRWXU);
    // change the working directory
    chdir(twd.c_str());
    
    // the exec file directory
    string exd = getexepath();

    // create a symbolic link to the original file
    #ifdef __linux
        cmd="cp -s " + file + " ./tf.fasta";
    #elif __APPLE__
        cmd="gcp -s " + file + " ./tf.fasta";
    #endif
    system(cmd.c_str());

    // index the file
    cmd="samtools faidx tf.fasta";
    system(cmd.c_str());

    // read the sequence list
    vector<string> seqs;
    ifstream infile("tf.fasta.fai");
    string line;
    while(getline(infile,line)){
        if (line.empty()) continue;
        stringstream ss(line);
        string seqname;
        ss >> seqname;
        seqs.emplace_back(seqname);
    }
    
    // extract the sequences and compute kmers
    for (auto seq : seqs){
        string subseq = "tf."+seq+".fasta";
        cmd="samtools faidx tf.fasta " + seq + " > " + subseq;
        system(cmd.c_str()); 
        cmd=exd+"/count_kmer -k " + to_string(ksize) + " -o " + to_string(osize) + " " + subseq + " > tf." + seq + ".kmers 2>/dev/null";
        system(cmd.c_str());
        
    }


    // compute the cartesian product of all sequences
    vector<tuple<string,string>> seqPairs;
    for (size_t i=0; i<seqs.size(); i++){
        for (size_t j=i+1; j<seqs.size(); j++){
            seqPairs.emplace_back(make_tuple(seqs[i],seqs[j]));
        }
    }


    omp_set_dynamic(0);
    omp_set_num_threads(cores);
    #pragma omp parallel for
    for (size_t i=0; i<seqPairs.size(); i++){
        string kmer1 = "tf." + get<0>(seqPairs[i]) + ".kmers";
        string kmer2 = "tf." + get<1>(seqPairs[i]) + ".kmers";
        cmd=exd+"/pair_kmer_sim " + kmer1 + " " + kmer2;
        string val = exec(cmd.c_str());
        stringstream ss;
        ss << val;
        ss >> val;

        ss.clear();
        ss.str(string(""));

        ss << get<0>(seqPairs[i]) << "\t" << get<1>(seqPairs[i]) << "\t" << val << "\n";
        cout << ss.str();
    }


    // change the previous working directory
    chdir(cwd.c_str());
   
    cmd = "rm -rf " + twd;
    system(cmd.c_str());

    return 0;
}

int run_mode2(string file1,string file2,int ksize,int osize,int cores){
    string cmd;

    // get the current working directory
    string cwd = getWorkDirectory();

    // make a temporary directory
    string twd = cwd + "/temporary_" + to_string(clock());
    mkdir(twd.c_str(),S_IRWXU);
    // change the working directory
    chdir(twd.c_str());

    // the exec file directory
    string exd = getexepath();

    // create a symbolic link to the original file1
    #ifdef __linux
        cmd="cp -s " + file1 + " ./tf1.fasta";
    #elif __APPLE__
        cmd="gcp -s " + file1 + " ./tf1.fasta";
    #endif
    system(cmd.c_str());

    // create a symbolic link to the original file2
    #ifdef __linux
        cmd="cp -s " + file2 + " ./tf2.fasta";
    #elif __APPLE__
        cmd="gcp -s " + file2 + " ./tf2.fasta";
    #endif
    system(cmd.c_str());

    // index the file
    cmd="samtools faidx tf1.fasta";
    system(cmd.c_str());
    cmd="samtools faidx tf2.fasta";
    system(cmd.c_str());

    // read the sequence list
    vector<string> seqs1,seqs2;
    ifstream infile1("tf1.fasta.fai");
    string line;
    while(getline(infile1,line)){
        if (line.empty()) continue;
        stringstream ss(line);
        string seqname;
        ss >> seqname;
        seqs1.emplace_back(seqname);
    }
    ifstream infile2("tf2.fasta.fai");
    while(getline(infile2,line)){
        if (line.empty()) continue;
        stringstream ss(line);
        string seqname;
        ss >> seqname;
        seqs2.emplace_back(seqname);
    }

    // extract the sequences and compute kmers
    for (auto seq : seqs1){
        string subseq = "tf1."+seq+".fasta";
        cmd="samtools faidx tf1.fasta " + seq + " > " + subseq;
        system(cmd.c_str());
        cmd=exd+"/count_kmer -k " + to_string(ksize) + " -o " + to_string(osize) + " " + subseq + " > tf1." + seq + ".kmers 2>/dev/null";
        system(cmd.c_str());
    }
    for (auto seq : seqs2){
        string subseq = "tf2."+seq+".fasta";
        cmd="samtools faidx tf2.fasta " + seq + " > " + subseq;
        system(cmd.c_str());
        cmd=exd+"/count_kmer -k " + to_string(ksize) + " -o " + to_string(osize) + " " + subseq + " > tf2." + seq + ".kmers 2>/dev/null";
        system(cmd.c_str());
    }


    // compute the cartesian product of all sequences
    vector<tuple<string,string>> seqPairs;
    for (size_t i=0; i<seqs1.size(); i++){
        for (size_t j=0; j<seqs2.size(); j++){
            seqPairs.emplace_back(make_tuple(seqs1[i],seqs2[j]));
        }
    }


    omp_set_dynamic(0);
    omp_set_num_threads(cores);
    #pragma omp parallel for
    for (size_t i=0; i<seqPairs.size(); i++){
        string kmer1 = "tf1." + get<0>(seqPairs[i]) + ".kmers";
        string kmer2 = "tf2." + get<1>(seqPairs[i]) + ".kmers";
        cmd=exd+"/pair_kmer_sim " + kmer1 + " " + kmer2;
        string val = exec(cmd.c_str());
        stringstream ss;
        ss << val;
        ss >> val;

        ss.clear();
        ss.str(string(""));

        ss << get<0>(seqPairs[i]) << "\t" << get<1>(seqPairs[i]) << "\t" << val << "\n";
        cout << ss.str();
    }


    // change the previous working directory
    chdir(cwd.c_str());

    cmd = "rm -rf " + twd;
    system(cmd.c_str());

    return 0;
}


int main(int argc, char *argv[]){
    // print help message if no arguments provided
    if (argc==1){
        help();
        exit(0);
    }

    // parse the command line
    int ksize, osize, cores;
    vector<string> files;
    parseCommandLine(argc,argv,ksize,osize,cores,files);

    // mode1
    if (files.size()==1){
        run_mode1(files[0],ksize,osize,cores);
    }
    
    // mode 2
    if (files.size()==2){
        run_mode2(files[0],files[1],ksize,osize,cores);
    }

    return 0;
}
