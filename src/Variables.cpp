#include <unistd.h>
#include <getopt.h>
#include <iostream>

#include "Variables.hxx"

struct option Option::options[] =
{
    {"help", 0, 0, 'h'},
    {"input", 1, 0, 'i'}
};
int Option::option_index = 0;
int Option::option = 0;
std::string Option::name_input_file = "";
