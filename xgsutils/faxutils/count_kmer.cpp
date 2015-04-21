#include <zlib.h>
#include <getopt.h>
#include <string>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_map>
using namespace std;

#include "kseq.h"

// help message
int help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    count_kmer [OPTIONS] <FASTA>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    count the kmers in the fasta sequence" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "    -k,--ksize    specify the length of the k-mer" << endl;
    cerr << "    -o,--osize    specify the overlap length of two contiguous k-mers" << endl;
    cerr << "    -h,--help     print help message" << endl;
    return 0;
}

// parse the command line
int parseCommandLine(int argc, char *argv[], int &ksize, int &osize, string &fastafile)
{
    ksize=16;
    osize=15;

    int c;
    while(1){
        static struct option long_options[] =
        {
            {"ksize",required_argument,0,'k'},
            {"osize",required_argument,0,'o'},
            {"help",no_argument,0,'h'},
            {0,0,0,0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "k:o:h", long_options, &option_index);

        // detect the end of the options
        if (c==-1) break;

        switch(c){
            case 'k':
                ksize = stoi(optarg);
                break;
            case 'o':
                osize = stoi(optarg);
                break;
            case 'h':
                help();
                exit(0);
            default:
                abort();
        }
    }
   
    if (ksize-osize<=0) {
        cerr << "invalid value for the option --osize" << endl;
        exit(1);
    }

    fastafile = argv[optind];
}


KSEQ_INIT(gzFile,gzread)

int main(int argc, char *argv[])
{
    int ksize,osize,shift;
    string fastafile;

    // print help message if no argument provides
    if (argc==1) { help(); exit(0); }
    
    // parse the command line
    parseCommandLine(argc,argv,ksize,osize,fastafile);
    shift = ksize-osize;


    gzFile fp;
    kseq_t *fasta;
    int l;
    
    fp = gzopen(fastafile.c_str(),"r");
    fasta = kseq_init(fp);


    unordered_map<string,long> kmers;
    while ((l=kseq_read(fasta)) >= 0){
        for (int i=0; i<strlen(fasta->seq.s); i+=shift){
            char kmer[ksize+1];
            memcpy(kmer,fasta->seq.s+i,ksize);
            kmer[ksize]=0;
            if (strlen(kmer)<ksize) continue;
            auto ptr = kmers.find(string(kmer));
            if (ptr==kmers.end()){
                kmers[string(kmer)] = 1;
            }else{
                ptr->second += 1;
            }
        }
    }
 
    // sort the kmers
    map<string,long> sorted_kmers(kmers.begin(),kmers.end());
    for (auto ptr=sorted_kmers.begin(); ptr!=sorted_kmers.end(); ptr++){
        cout << setw(ksize) << ptr->first << "\t" << setw(12) << ptr->second << endl;
    }

    cerr << "total " << sorted_kmers.size() << " kmers" << endl;
}
