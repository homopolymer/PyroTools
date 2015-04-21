#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <regex>
#include <map>
#include <vector>
#include <iomanip>
using namespace std;

#include <getopt.h>

// help message
int help(){
    cerr << "SYNOPSIS" << endl;
    cerr << "    read_haplo_matrix [OPTIONS] <sam>" << endl;
    cerr << "" << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "    Print a matrix with the measures between reads and haplotypes" << endl;
    cerr << "" << endl;
    cerr << "OPTIONS" << endl;
    cerr << "    -m    specify the metrics to be reported" << endl;
    cerr << "              valid arguments:" << endl;
    cerr << "                  nm    the number of alignment difference" << endl;
    cerr << "                iden    the identity of alignment" << endl;
    cerr << "                like    the likelihood score" << endl;
    cerr << "                post    the posterior score" << endl;
    cerr << "    -v    print the runtime message" << endl;
    cerr << "    -h    print help message" << endl;
    return 0;
}

// parse the command line
int parseCommandLine(int argc, char *argv[], string &sam, string &metric, bool &verbose){
    metric = "nm";
    verbose = false;
    
    int c;
    while((c=getopt(argc,argv,"m:vh"))!=-1){
        switch(c){
            case 'v':
                verbose=true;
                break;
            case 'm':
                metric = optarg;
                break;
            case 'h':
                help();
                exit(0);
                break;
            default:
                abort();
        }
    }
  
    sam = argv[optind];
    return 0;
}

int main(int argc, char *argv[]){
    // print help if no argument provided
    if (argc == 1){
        help();
        exit(0);
    }

    bool verbose;
    string sam, metric;
    parseCommandLine(argc, argv, sam, metric, verbose);

    // parse the sam file
    ifstream input(sam);
    string line;

    int num_ref = 0;
    vector<string> ref_seqs;
    map<string,map<string,double>> matrix;
    int rw=0, fw=0;
    
    int cnt = 0;
    regex sqex("^@SQ\\s+SN:([\\w[:punct:]]+)");
    regex rcex("^@");
    regex nmex("NM:i:([0-9]+)");
    smatch matches;
    ssub_match ssm;
    while(getline(input,line)){
        if (regex_search(line, matches, sqex)){
            ssm = matches[1];
            ref_seqs.emplace_back(ssm.str());
            num_ref++;
        }else if (!regex_search(line,rcex)){
            if (verbose)
                cerr << "\r" << cnt++ << flush;            

            stringstream ss(line);

            string query,target;
            string cigar;
            int flag, mapq, pos;
            int nm, len;
            ss >> query >> flag >> target >> mapq >> pos >> cigar;
            
            if (regex_search(line, matches, nmex)){
                ssm = matches[1];
                nm = stoi(ssm.str());
            }

            // parse cigar
            len = 0;
            size_t prev=0, curr;
            while ((curr=cigar.find_first_of("HSM=XDI",prev))!=string::npos){
                if (curr>prev){
                    len += stoi(cigar.substr(prev,curr-prev));
                }
                prev = curr+1;
            }
    
            if (metric=="nm") matrix[query][target] = nm;
            if (metric=="iden") matrix[query][target] = 1-nm/(len+0.0);        
            if (metric=="like" || metric=="post") matrix[query][target] = exp((len-nm)*log(0.99)+nm*log(0.01));
            
            if (rw < query.length())
                rw = query.length();
        }
    }
    input.close();

    if (verbose)
        cerr << endl;

    if (metric=="post"){
        for (auto r : matrix){
            double z = 0;
            for (auto t : ref_seqs){
                z += matrix[r.first][t];
            }
            for (auto t : ref_seqs){
                matrix[r.first][t] /= z;
            }
        }
    }

    cout << setw(rw);
    for (auto t : ref_seqs){
        fw = t.length();
        cout << "\t" << setw(fw) << t;
    }
    cout << endl;

    for (auto r : matrix){
        cout << setw(rw) << r.first;
        for (auto t : ref_seqs){
            fw = t.length();
            cout << "\t" << setw(fw) << matrix[r.first][t];
        }
        cout << endl;
    }

    return 0;  
}
