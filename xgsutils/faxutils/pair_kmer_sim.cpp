#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
using namespace std;


// help message
int help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    pair_kmer_sim <kmer1> <kmer2> " << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    compute the cosine similarity of two kmer spectrums" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "" << endl;
    return 0;
}


int main(int argc, char *argv[])
{
    // print help message if no arguments provided
    if (argc==1) { help(); exit(0); }

    // kmer spectrum one
    string ksf1(argv[1]);
    // kmer spectrum two
    string ksf2(argv[2]);

    // load kmer spectrum
    unordered_map<string,float> ks1,ks2;
    // kmer spectrum one
    string kmer;
    int    count;
    ifstream infile1(ksf1);
    while (infile1 >> kmer >> count){
        ks1[kmer] = count;
    }

    ifstream infile2(ksf2);
    while (infile2 >> kmer >> count){
        ks2[kmer] = count;
    }

    // normalized
    float z = 0.0;
    for (auto x : ks1){ z += x.second; }
    for (auto &x: ks1){ x.second /= z; }

    z = 0;
    for (auto x : ks2){ z += x.second; }
    for (auto &x: ks2){ x.second /= z; }

    // kmer dictionary
    unordered_set<string> kmer_dict;
    for (auto x : ks1) { kmer_dict.emplace(x.first); }
    for (auto x : ks2) { kmer_dict.emplace(x.first); }

    // compute the sum
    float s1=0,s2=0,ss=0;
    auto ks1end = ks1.end();
    auto ks2end = ks2.end();
    for (auto k : kmer_dict){
        int occ=0;
        if (ks1.find(k)!=ks1end) {
            s1 += ks1[k]*ks1[k];
            occ++;
        }
        if (ks2.find(k)!=ks2end) {
            s2 += ks2[k]*ks2[k];
            occ++;
        }
        if (occ==2){
            ss += ks1[k]*ks2[k];
        }
    }

    // print result
    cout << (ss/sqrt(s1*s2)) << endl;

    return 0;
}
