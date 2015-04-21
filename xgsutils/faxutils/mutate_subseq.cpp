// reference
// http://lh3lh3.users.sourceforge.net/parsefastq.shtml

#include <zlib.h>
#include <getopt.h>
#include <cstdlib>
#include <regex>
#include <string>
#include <iostream>
using namespace std;

#include "kseq.h"

// help message
int help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    mutate_subseq <gf> <mf> <roi> <name>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    replace the genome with the mutant in the region of interest" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "" << endl;
    return 0;
}


KSEQ_INIT(gzFile,gzread)

int main(int argc, char *argv[])
{
    gzFile fp,fp2;
    kseq_t *genome, *mutant;
    int l;

    if (argc==1)
    {
        help();
        exit(0);
    }
   
    // parse the region of interest
    regex pattern("(\\w+):(\\d+)-(\\d+)");
    smatch match;
    regex_match(string(argv[3]),match,pattern);
    // chromosome name
    ssub_match sub_match = match[1];
    string chr = sub_match.str();
    // left position
    sub_match = match[2];
    int left = stoi(sub_match.str())-1;
    // right position
    sub_match = match[3];
    int right = stoi(sub_match.str())-1;

    // load mutant 
    fp2 = gzopen(argv[2],"r");
    mutant = kseq_init(fp2);
    kseq_read(mutant);
    
    // load genome
    fp = gzopen(argv[1],"r");
    genome = kseq_init(fp);
 
    while ((l=kseq_read(genome)) >= 0){
        if ( chr == string(genome->name.s) ){
            string head,tail;
            for (int i=0; i<left; i++){
                head += genome->seq.s[i];
            }
            for (int i=right+1; i<strlen(genome->seq.s); i++){
                tail += genome->seq.s[i];
            }
            string mutseq;
            mutseq = head + string(mutant->seq.s) + tail;
//            for (int i=left; i<=right; i++){
//                genome->seq.s[i] = mutant->seq.s[i-left];
//            }
            cout << ">" << argv[4] << endl;
            cout << mutseq << endl;
            break;
        }
    }

    // clean everything
    kseq_destroy(genome);
    gzclose(fp);

    kseq_destroy(mutant);
    gzclose(fp2);
 
    return 0;
}
