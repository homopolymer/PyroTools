#include "../c++/kseq.h"
#include "../c++/utils.h"
using namespace XgsUtils;

#include <ctime>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <regex>
#include <vector>
#include <utility>
using namespace std;

#include <zlib.h>
#include <getopt.h>

// help message
int help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    haplo_cluster [OPTIONS] <fas>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    Output the clusters of the haplotype sequences in the given fasta file" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "    -t    specify the threshold of the distance measure (default:0.01) [FLT]" << endl;
    cerr << "    -v    print the runtime message" << endl;
    cerr << "    -h    print help message" << endl;
    return 0;
}

// parse the command line
int parseCommandLine(int argc, char *argv[], string &fas, double &thres, bool &verbose){
    thres = 0.015;
    verbose = false;

    int c;
    while((c=getopt(argc,argv,"t:vh"))!=-1){
        switch(c){
            case 'v':
                verbose=true;
                break;
            case 't':
                thres = stod(optarg);
                break;
            case 'h':
                help();
                exit(0);
                break;
            default:
                abort();
        }
    }
 
    fas = DirectoryUtils::AbsolutePath(argv[optind]);

    return 0;
}

KSEQ_INIT(gzFile,gzread)

int main(int argc, char *argv[]){
    // print help if no argument provided
    if (argc==1){
        help();
        exit(0);
    }

    string fas;
    double thres;
    bool verbose;
    
    // parse the command line
    parseCommandLine(argc,argv,fas,thres,verbose);

    // read the fasta sequences
    map<string,string> haplotypes;

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fas.c_str(),"r");
    seq = kseq_init(fp);
    while((l=kseq_read(seq))>=0){
        haplotypes[seq->name.s] = seq->seq.s;
    }
    gzclose(fp);
    kseq_destroy(seq);

    if (haplotypes.size()>1){
      // the current working directory
      string cwd = DirectoryUtils::CurrentWorkDir();
    
      // make a temporary directory
      string twd = cwd + "/temporary_" + to_string(clock()); 
      DirectoryUtils::MkDir(twd);
      DirectoryUtils::ChDir(twd);
    
      // the directory to the executing code
      string ewd = DirectoryUtils::GetExePath();

      string cmd;
      // compute the distance of pairs of haplotypes
      cmd = ewd + "/seq_ssw_sim " + fas + " 2>/dev/null | awk '{$3=1-$3;print}' > dist.txt";
      SystemUtils::Run(cmd);
   
      string result;
      cmd = "python " + ewd + "/hcluster.py -t " + to_string(thres) + " dist.txt 2>/dev/null";
      SystemUtils::Run(cmd,result);

      // parse the clustering result
      vector<pair<int,string>> clusters;
      stringstream ss(result);
      string line;
      while(getline(ss,line)){
          string item;
          vector<pair<int,string>> seqs;
          stringstream buf(line);
          buf >> item;
          while (!buf.eof()){
              buf >> item;
              smatch matches;
              regex_search(item,matches,regex("(\\d+)"));
              ssub_match ssm = matches[1];
              int key = stoi(ssm.str());
              seqs.emplace_back(make_pair(key,item));
          }
          sort(seqs.begin(), seqs.end(), [](pair<int,string> &a, pair<int,string> &b)->bool{
              return a.first<b.first;
          });
          clusters.emplace_back(*seqs.begin());
      }
    
      sort(clusters.begin(), clusters.end(), [](pair<int,string> &a, pair<int,string> &b)->bool{
          return a.first<b.first;
      });

      for (auto s : clusters){
          cout << ">" << s.second << endl;
          cout << haplotypes[s.second] << endl;
      }

      DirectoryUtils::ChDir(cwd);
      DirectoryUtils::RmDir(twd);

    }else{
        for (auto h : haplotypes){
            cout << ">" << h.first << endl;
            cout << h.second << endl;
        }
    }
    
    return 0;
}
