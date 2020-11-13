#include <iostream>
#include <sstream>
#include <unistd.h>
#include <getopt.h>

#include "TChain.h"

#include "Functions.hxx"
#include "Variables.hxx"

void Function::Help()
{
    std::cout<<"Usage: ./Lownu [OPTION]"<<std::endl;
    std::cout<<"  --help, -h                   : show this help message"<<std::endl;
    std::cout<<"  --inputfile, -i [file]       : input files"<<std::endl;
}

bool Function::Parsing(int argc, char * argv[])
{
    if(argc == 1)
    {
        Help();
        return false;
    }
    while((Option::option = getopt_long(argc, argv, "hi:", Option::options, &Option::option_index))!=EOF)
    {
        switch(Option::option)
        {
            case 'h' :
                {
                    Help();
                    return false;
                }
            case 'i' :
                {
                    Option::name_input_file = optarg;
                    break;
                }
        }
    }
    return true;
}

bool Function::InputFile(const std::unique_ptr<TChain>& tree)
{
    if(Option::name_input_file.size() != 0)
    {
        tree->Add(Option::name_input_file.c_str());
    }
    else
    {
        std::cout<<"No input file, quit"<<std::endl;
        return false;
    }
    return true;
}
