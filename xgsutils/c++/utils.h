#include <cstdlib>
#include <climits>
#include <string>
#include <ctime>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <iostream>
using namespace std;

#include <sys/stat.h>
#include <unistd.h> // gnu library
#include <libgen.h> // for basename, http://linux.die.net/man/3/basename

#ifdef __APPLE__
    #include <mach-o/dyld.h>
#endif

// TimeUtils
// DirectoryUtils
// VerboseUtils
// SystemUtils
// TemporaryUtils

namespace XgsUtils{

// utilities for the time display
class TimeUtils{
    public:
        TimeUtils(){}
    public:
        static string CurrentTime(){
            time_t result = time(NULL);
            stringstream ss;
            ss << asctime(localtime(&result));
            return ss.str();
        }
}; // class TimeUtils


// utilities for the directory information
class DirectoryUtils{
    public:
        DirectoryUtils(){}
    public:
        // get the absolute path of a given path
        static string AbsolutePath(string p){
            int path_max;
            // redirect to http://linux.die.net/man/3/realpath
            #ifdef PATH_MAX
                path_max = PATH_MAX;
            #else
                path_max = pathconf((char*)p.c_str(), _PC_PATH_MAX);
                if (path_max < 0)
                    path_max = 4096;
            #endif
     
            // temporary result
            char result[path_max];
            realpath(p.c_str(), result);
            
            return string(result);
        }

        // get the name of the contain directory of a given file
        static string DirName(string f){
            // canonicalize the path
            string p2f = AbsolutePath(f);
            
            char *result = dirname((char*)p2f.c_str());
            return string(result);
        }

        // get the name of a given file in the form of the relative/absolute representation
        static string FileName(string f){
            // canonicalize the path
            string p2f = AbsolutePath(f);
            
            char *result = basename((char*)p2f.c_str());
            return string(result);
        }
     
        // get the current working directory
        static string CurrentWorkDir(){
            char result[PATH_MAX];
            getcwd(result, PATH_MAX);
            return string(result);
        }

        // get the path to the executing program
        static string GetExePath(){
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

        // make a temporary directory
        static void MkDir(string dir){
            mkdir(dir.c_str(),S_IRWXU);
        }

        // change the working directory
        static void ChDir(string dir){
            chdir(dir.c_str());
        }

        // remove the directory
        static void RmDir(string dir){
            string cmd = "rm -rf " + AbsolutePath(dir);
            system(cmd.c_str());
        }
}; // class DirectoryUtils

// utilities for verbose
class VerboseUtils{
    public:
        VerboseUtils(){}
    public:
        static void Verbose(string head, string msg){
            string t = TimeUtils::CurrentTime();
            stringstream ss;
            ss << "[" << t << " " << head << "]" << " " << msg;
            cerr << ss.str() << endl;
        }
}; // class VerboseUtils

// utilities for running the external command
class SystemUtils{
    public:
        SystemUtils(){}
    public:
        // run the external command and get the executation status
        static int Run(string &cmd){
            return system(cmd.c_str());
        }
        
        // run the external command and retrieve the output in string format
        static int Run(string &cmd, string &result){
            result = "";
            FILE* pipe = popen(cmd.c_str(), "r");
            if (!pipe) return -1;
            char buffer[128];
            while(!feof(pipe)) {
                if(fgets(buffer, 128, pipe) != NULL)
                    result += buffer;
            }
            pclose(pipe);
            return 0;
        }
}; // class SystemUtils

// utilities for temporary file and directory
class TemporaryUtils{
    public:
        TemporaryUtils(){}
    public:
        // get a temporary name
        static string TemporaryName(){
            // cpu clock
            auto t = clock();
            // random number
            srand(time(0));
            auto r = rand();
            // temporary name
            stringstream ss;
            ss << "temp_" << t << "_" << r;
            
            return ss.str();
        }
        
        // create a temporary directory
        static int CreateTemporaryDir(string dir_name){
            DirectoryUtils::MkDir(dir_name);
            return 0;
        }

        // remove a temporary directory
        static int RemoveTemporaryDir(string dir_name){
            DirectoryUtils::RmDir(dir_name);
            return 0;
        }

        // go to a temporary directory
        static int GotoTemporaryDir(string dir_name){
            DirectoryUtils::ChDir(dir_name);
            return 0;
        }
}; // class Temporary
}; // namespace
