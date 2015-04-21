#ifndef GENERICSEQUENCEGLOBAL_H
#define GENERICSEQUENCEGLOBAL_H

#include <map>
#include <vector>
#include <string>
#include <limits>
#include <sys/stat.h>

namespace GenericSequenceTools{

    // the allowed letters in DNA sequence
    static char         AlphabetChar[] = {'A','C','G','T','-','N'};
    static std::string  AlphabetString[] = {"A","C","G","T","-","N"};
    static int numAlphabet   = 6;
    static int numNucleotide = 4;

    typedef char Nucl;
    typedef char Alpha;

    static char Ade = 'A';
    static char Cyt = 'C';
    static char Gua = 'G';
    static char Thy = 'T';
    static char Spa = '-';
    static char Amb = 'N';
    static char Nil = '?';
    static char Zero = char(0);

typedef int SequenceBase;
#define ADENINE   0
#define CYTOSINE  1
#define GUANINE   2
#define THYMINE   3
#define SPACE     4
#define AMBIGUOUS 5

static std::map<Alpha, SequenceBase> Alpha2Base = {{'A',ADENINE}, {'C',CYTOSINE}, {'G', GUANINE},
                                                   {'T',THYMINE}, {'-',SPACE}, {'N',AMBIGUOUS}};

    // the mininum length of homopolymer being considered
    static int minHomopolymerSize = 1;
    // the middle length of homopolymer being considered
    static int midHomopolymerSize = 2;
    // the maximum length of homopolymer being considered
    static int maxHomopolymerSize = 6;
    // the number of homopolymer being consider
    static int numHomopolymerSize = maxHomopolymerSize+1;




    // generic types that are used
    typedef std::vector<int>            VectorInteger;
    typedef std::vector<long double>    VectorDouble;

    typedef std::vector<VectorInteger>  Matrix2Integer;
    typedef std::vector<VectorDouble>   Matrix2Double;

    typedef std::vector<Matrix2Integer> Matrix3Integer;
    typedef std::vector<Matrix2Double>  Matrix3Double;

    typedef std::vector<Matrix3Double>  Matrix4Double;

    VectorDouble&  createVector(int size, long double value);
    VectorInteger& createVector(int size, int value);

    Matrix2Double&  createMatrix(int m, int n, long double value);
    Matrix2Integer& createVector(int m, int n, int value);

    Matrix3Double&  createMatrix(int m, int n, int k, long double value);
    Matrix3Integer& createMatrix(int m, int n, int k, int value);

    Matrix4Double& createMatrix(int m, int n, int k, int l, long double value);

    void setVectorValue(VectorInteger& vec, int size, int value);
    void setVectorValue(VectorDouble&  vec, int size, long double value);

    void setMatrixValue(Matrix2Integer& mat, int m, int n, int value);
    void setMatrixValue(Matrix2Double&  mat, int m, int n, long double value);

    void setMatrixValue(Matrix3Integer& mat, int m, int n, int k, int value);
    void setMatrixValue(Matrix3Double&  mat, int m, int n, int k, long double value);

    void setMatrixValue(Matrix4Double&  mat, int m, int n, int k, int l, long double value);




#define DOUBLE_NEGATIVE_INFINITY std::numeric_limits<long double>::lowest()
#define DOUBLE_POSITIVE_INFINITY std::numeric_limits<long double>::max()



class GenericRegressor
{
    public:
        GenericRegressor(){}

    public:
        static double linearfit(std::vector<double>& x, std::vector<double>& y, double& intercept, double& slop);
        static double linearslop(std::vector<double>& x, std::vector<double>& y);
};


std::string GetExecPath();

bool FileExist(const std::string& filename);

int RunCmdGetReturn(std::string& cmd, std::string& result);

std::string CurrentTime();

}   // namespace

#endif // GENERICSEQUENCEGLOBAL_H
