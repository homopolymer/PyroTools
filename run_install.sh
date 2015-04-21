#!/bin/bash

#help message
help(){
  >&2 echo "SYNOPSIS"
  >&2 echo "    bash ./run_install.sh [OPTIONS]"
  >&2 echo ""
  >&2 echo "DESCRIPTION"
  >&2 echo "    Install PyroTools, xgsutils and the dependencies"
  >&2 echo ""
  >&2 echo "OPTIONS"
  >&2 echo "    --all        install the whole package"
  >&2 echo "    --pyrotools  build the PyroTools pakcage"
  >&2 echo "    --xgsutils   build the xgsutils package"
  >&2 echo "    --bamtools   build the bamtools package"
  >&2 echo "    --nlopt      build the nlopt package"
  >&2 echo "    --help       print help message"
  exit 0
}

# print help message if no arguments provided
if [ $# = 0 ];then help;fi

# parse the command line
PARSED_OPTIONS=$(getopt -n "$0" -o h --long all,pyrotools,xgsutils,bamtools,nlopt,help -- "$@")

# Bad arguments, something has gone wrong with the getopt command
if [ $? -ne 0 ]
then
  help
fi

# A little magic, necessary when using getopt
eval set -- "$PARSED_OPTIONS"

BUILD_PYROTOOLS=false
BUILD_XGSUTILS=false
BUILD_BAMTOOLS=false
BUILD_NLOPT=false
# Now goes through all the options with a case and using shift to analyse 1 argument at a time
# $1 identifies the first argument, and when we use shift we discard the first argument, 
# so $2 becomes $1 and goes again through the case
while true; do
  case "$1" in
    --help )
      help
      shift;;
    --all )
      BUILD_PYROTOOLS=true
      BUILD_XGSUTILS=true
      BUILD_BAMTOOLS=true
      BUILD_NLOPT=true
      shift;;
    --pyrotools )
      BUILD_PYROTOOLS=true
      shift;;
    --xgsutils )
      BUILD_XGSUTILS=true
      shift;;
    --bamtools )
      BUILD_BAMTOOLS=true
      shift;;
    --nlopt )
      BUILD_NLOPT=true
      shift;;
    -- )
      shift
      break;;
  esac
done

# build bamtools
if [ "$BUILD_BAMTOOLS" = true ]
then
  cd bamtools
  if [ -d "build" ];then rm -rf build;fi
  mkdir build && cd build && cmake .. && make && cd ../..
fi

# build nlopt
if [ "$BUILD_NLOPT" = true ]
then
  cd nlopt && ./configure --prefix="$PWD" && make && make install && cd ..
fi

# build PyroTools
if [ "$BUILD_PYROTOOLS" = true ]
then
  if [ -d "build" ];then rm -rf build;fi
  mkdir build && cd build && cmake .. && make && cd ..
  if [[ ! $PATH =~ $PWD/bin ]]
  then
    if [[ $OSTYPE =~ linux* ]];then echo -e "# add PyroTools to the PATH variable\nexport PATH=$PWD/bin:\$PATH\n" >> ~/.bashrc;source ~/.bashrc;fi
    if [[ $OSTYPE =~ darwin* ]];then echo -e "# add PyroTools to the PATH variable\nexport PATH=$PWD/bin:\$PATH\n" >> ~/.profile;source ~/.profile;fi
  fi
fi

# build xgsutils
if [ "$BUILD_XGSUTILS" = true ]
then
  cd xgsutils && g++ -std=c++11 -o xgsutils xgsutils.cpp && cd ..
#  if [[ ! $PATH =~ $PWD/xgsutils ]]
#  then
#    if [[ $OSTYPE =~ linux* ]];then echo -e "# add xgsutils to the PATH variable\nexport PATH=$PWD/xgsutils:\$PATH\n" >> ~/.bashrc;fi
#    if [[ $OSTYPE =~ darwin* ]];then echo -e "# add xgsutils to the PATH variable\nexport PATH=$PWD/xgsutils:\$PATH\n" >> ~/.profile;fi
#  fi
fi
