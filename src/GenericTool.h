#ifndef GENERICTOOL_H
#define GENERICTOOL_H

namespace GenericSequenceTools
{

class GenericAbstractTool
{
    public:
        GenericAbstractTool(){}
        virtual ~GenericAbstractTool(){}

    public:
        virtual int Help() = 0;
        virtual int Run(int argc, char* argv[]) = 0;
};

}
#endif // GENERICTOOL_H
