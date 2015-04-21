#include <limits>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <forward_list>
#include <map>
#include <set>
#include <tuple>
#include <regex>
#include <iterator>
using namespace std;

#include <getopt.h>


// xgsutils header
#include "../c++/utils.h"
using namespace XgsUtils;


// help message
void help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    long_read_cns_table [OPTION] <SAM>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    It reports the LongRead-Consensus occurrence matrix from the give SAM file" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "    -d,--diff        specify the NM threshold [INT]" << endl;
    cerr << "    -c,--complete    specify the complete list of consensus sequences" << endl;
    cerr << "    -t               specify threshold to filter long reads hit by equal to or less than t CNS [INT]" << endl;
    cerr << "    -T               specify threshold to filter long reads of average hamming distance equal to or larger than T [INT]" << endl;
    cerr << "    -h,--help        print help message" << endl;
}


// command line parser
void parseCommandLine(int argc, char *argv[], string& sam_file, int& diff, int& tt, int& tT, string& cns_file){
    diff=numeric_limits<int>::max();
    sam_file="";
    cns_file="";
    tt=0;
    tT=100;

    int c;
    while (1){
        static struct option long_options[] =
        {
            {"diff",    required_argument,0,'d'},
            {"complete",required_argument,0,'c'},
            {"help",    no_argument,      0,'h'},
            {0,0,0,0}
        };

        // getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "d:c:t:T:h", long_options, &option_index);
        
        // detect the end of the options
        if (c==-1) break;

        switch(c){
            case 'd':
                diff = stoi(optarg);
                break;
            case 'c':
                cns_file = optarg;
                break;
            case 't':
                tt = stoi(optarg);
                break;
            case 'T':
                tT = stoi(optarg);
                break;
            case 'h':
                help();
                exit(0);
                break;
            default:
                VerboseUtils::Verbose("long_read_cns_table","invalid option...abort");
                abort();
        }
    }

    sam_file = argv[optind];
}

// main function
int main(int argc, char *argv[]){
    // print help message, if no argument provided
    if (argc==1){
        help();
        exit(0);
    }

    // parse command line
    int diff, tt, tT;
    string cns_file;
    string sam_file;
    parseCommandLine(argc,argv,sam_file,diff,tt,tT,cns_file);

    // executing program directory
    string exd = DirectoryUtils::GetExePath();

    stringstream ssCMD, ssResult;
    string cmd, result;
    // collect long reads
    ssCMD << exd << "/../bamutils/xBamFilterByNM -n " << diff << " " << sam_file << " | cut -f3 | sort | uniq";
    cmd = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    SystemUtils::Run(cmd,result);

    ssResult.str(result);
    forward_list<string> long_reads;
    for(string line;getline(ssResult,line);){
        if (line.empty()) continue;
        long_reads.emplace_front(line);
    }   

    // collect consensus sequences
    map<string,forward_list<string>> cnss;
    if (!cns_file.empty()){
        ifstream input(cns_file);
        for (string cns; getline(input,cns);){
            size_t p = cns.find_last_of("_");
            string spot = cns.substr(0,p);
            cnss[spot].emplace_front(cns);
        }
    }else{
        ssCMD << exd << "/../bamutils/xBamFilterByNM -n " << diff << " " << sam_file << " | cut -f1 | sort | uniq";
        cmd = ssCMD.str();
        ssCMD.clear();
        ssCMD.str("");
        SystemUtils::Run(cmd,result);

        ssResult.clear();
        ssResult.str(result);
        for (string cns;getline(ssResult,cns);){
            if (cns.empty()) continue;
            size_t p = cns.find_last_of("_");
            string spot = cns.substr(0,p);
            cnss[spot].emplace_front(cns);
        }
    }

    // collect LongRead-Cns occurrence
    map<string,map<string,tuple<set<string>,double>>> lr_cns_table;
    
    ssCMD << exd << "/../bamutils/xBamFilterByNM -n " << diff << " " << sam_file;
    cmd = ssCMD.str();
    ssCMD.clear();
    ssCMD.str("");
    SystemUtils::Run(cmd,result);

    ssResult.clear();
    ssResult.str(result);
    for (string rec;getline(ssResult,rec);){
        if (rec.empty()) continue;
        string spot,cns,lr,cigar,dummy1,dummy2;
        int flag;
        stringstream buf;
        buf.str(rec);
        buf >> cns >> flag >> lr >> dummy1 >> dummy2 >> cigar;
        // get spot
        size_t p = cns.find_last_of("_");
        spot = cns.substr(0,p);
        // skip if no map
        if (flag&0x4) continue;
        // get nm
        int nm;
        regex nmex("NM:i:([0-9]+)");
        ssub_match ssm;
        smatch matches;
        if (regex_search(rec, matches, nmex)){
            ssm = matches[1];
            nm = stoi(ssm.str());
        }

        // get entire length
        regex enex("([0-9]+)[M=DI]");
        int en = 0;
        auto en_begin = sregex_iterator(cigar.begin(), cigar.end(), enex);
        auto en_end = sregex_iterator();
        for (sregex_iterator i = en_begin; i != en_end; i++){
            smatch match = *i;
            en += stoi(match.str());
        }
  
        // compute identity
        double iden = (en-nm) / (en+0.0);
        
        // save record
        if (lr_cns_table.count(lr)){
            if (lr_cns_table[lr].count(spot)){
                int nm0 = get<1>(lr_cns_table[lr][spot]);
                if (nm<nm0){
                    lr_cns_table[lr][spot] = make_tuple(set<string>{cns},nm);
                }else if (nm==nm0){
                    get<0>(lr_cns_table[lr][spot]).emplace(cns);
                }
            }else{
                lr_cns_table[lr][spot] = make_tuple(set<string>{cns},nm);
            }
        }else{
            lr_cns_table[lr][spot] = make_tuple(set<string>{cns},nm);
        }
    }

    // output
    // title line
    cout << "#";
    for (auto a : cnss)
        for (auto cns : a.second)
            cout << "\t" << cns;
    cout << endl;

    // long-read line
    for (auto lr : long_reads){
        // filtration
        int lrtt, lrtT=0, ntT=0;
        lrtt = lr_cns_table[lr].size();
        if (lrtt<=tt) continue;
        for (auto a : lr_cns_table[lr]){
            for ( auto b : get<0>(a.second)){
                lrtT += get<1>(a.second);
                ntT += 1;
            }
        }
        if (lrtT>=ntT*tT) continue;

        // report
        cout << lr;
        for (auto a : cnss){
            string spot = a.first;
            for (auto cns : a.second){
                auto b = get<0>(lr_cns_table[lr][spot]);
                auto c = get<1>(lr_cns_table[lr][spot]);
                if (b.count(cns)) cout << "\t" << c+1;
                else cout << "\t" << 0;
            }
        }
        cout << endl;
    }

    return 0;
}
