#include <map>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

#include <getopt.h>
#include <zlib.h>
#include "kseq.h"

int help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    filter_contig_by_count <fastafile> <readcount>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    Filter out the contigs with few read hits" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "    -t,--thres    specify the threshold of read hit" << endl;
    cerr << "    -h,--help     print help message" << endl;
    return 0;
}

int parseCommandLine(int argc, char *argv[], string &fasta, string &readhit, int &t){
    t = 0;
    
    int c;
    while(1){
        static struct option long_options[] =
        {
            {"thres",required_argument,0,'t'},
            {"help",no_argument,0,'h'},
            {0,0,0,0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "t:h", long_options, &option_index);

        // detect the end of the options
        if (c==-1) break;
        
        switch(c){
            case 't':
                t=stoi(optarg);
                break;
            case 'h':
                help();
                exit(0);
            default:
                abort();
        }
    }

    fasta = argv[optind++];
    readhit = argv[optind];

    return 0;
}


KSEQ_INIT(gzFile,gzread)

int main(int argc, char *argv[]){
    // print help if no arguments provided
    if (argc==1){
        help();
        exit(0);
    }

    string fasta,readhit;
    int T;
    // parse the command line
    parseCommandLine(argc,argv,fasta,readhit,T);

    gzFile fp;
    kseq_t *seq;
    int l;

    // load hit information
    map<string,int> contig_counts;
    ifstream infile(readhit);
    string contig;
    int count;
    while(infile >> contig >> count){
        contig_counts[contig] = count;
    }
    infile.close();

    // read fasta file
    fp=gzopen(fasta.c_str(),"r");
    seq=kseq_init(fp);
    while((l=kseq_read(seq))>=0){
        string contig(seq->name.s);
        if (contig_counts[contig]>T){
            cout << ">" << seq->name.s << endl;
            cout << seq->seq.s << endl;
        }
    }
    kseq_destroy(seq);
    gzclose(fp);

    return 0;
}
