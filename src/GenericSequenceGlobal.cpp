#include "GenericSequenceGlobal.h"
using namespace GenericSequenceTools;

#include <algorithm>
#include <vector>
#include <numeric>
#include <iostream>
#include <ctime>
#include <sstream>
#include <string>
using namespace std;

#ifdef __linux
    #include <unistd.h>
#elif __APPLE__
    #include <mach-o/dyld.h>
#endif

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <libgen.h>
#include <sys/stat.h>


void GenericSequenceTools::setVectorValue(VectorInteger &vec, int size, int value)
{
    vec = VectorInteger(size, value);
}

void GenericSequenceTools::setVectorValue(VectorDouble &vec, int size, long double value)
{
    vec = VectorDouble(size, value);
}

void GenericSequenceTools::setMatrixValue(Matrix2Integer &mat, int m, int n, int value)
{
    mat = Matrix2Integer(m, VectorInteger(n, value));
}

void GenericSequenceTools::setMatrixValue(Matrix2Double &mat, int m, int n, long double value)
{
    mat = Matrix2Double(m, VectorDouble(n, value));
}

void GenericSequenceTools::setMatrixValue(Matrix3Integer &mat, int m, int n, int k, int value)
{
    mat = Matrix3Integer(m, Matrix2Integer(n, VectorInteger(k, value)));
}

void GenericSequenceTools::setMatrixValue(Matrix3Double &mat, int m, int n, int k, long double value)
{
    mat = Matrix3Double(m, Matrix2Double(n, VectorDouble(k, value)));
}

void GenericSequenceTools::setMatrixValue(Matrix4Double &mat, int m, int n, int k, int l, long double value)
{
    mat = Matrix4Double(m, Matrix3Double(n, Matrix2Double(k, VectorDouble(l, value))));
}


double GenericRegressor::linearslop(vector<double> &x, vector<double> &y)
{
    const auto n    = x.size();
    const auto s_x  = std::accumulate(x.begin(), x.end(), 0.0);
    const auto s_y  = std::accumulate(y.begin(), y.end(), 0.0);
    const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);

    return a;
}

double GenericRegressor::linearfit(std::vector<double> &x, std::vector<double> &y, double &intercept, double &slop)
{
    slop = linearslop(x, y);

    const auto n   = x.size();
    const auto s_x = std::accumulate(x.begin(), x.end(), 0.0);
    const auto s_y = std::accumulate(y.begin(), y.end(), 0.0);
    const auto m_x = s_x / n;
    const auto m_y = s_y / n;

    intercept = m_y - m_x*slop;

    double residues = 0.0;
    for (int i=0; i<n; ++i)
    {
        residues += (y[i]-x[i]*slop)*(y[i]-x[i]*slop);
    }

    return residues;
}

string GenericSequenceTools::GetExecPath()
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

bool GenericSequenceTools::FileExist(const string &filename)
{
    struct stat buffer;
    return (stat (filename.c_str(), &buffer) == 0);
}


int GenericSequenceTools::RunCmdGetReturn(string &cmd, string &result)
{
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

string GenericSequenceTools::CurrentTime()
{
    string value;
    time_t result = time(NULL);
    stringstream ss;
    ss << asctime(localtime(&result));
    getline(ss, value);
    return value;
}
