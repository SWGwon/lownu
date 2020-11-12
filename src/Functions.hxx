#ifndef FUNCTION_H
#define FUNCTION_H

#include "TChain.h"

namespace Function
{
    void Help();

    bool Parsing(int argc, char * argv[]);

    bool InputFile(TChain * tree);
}

#endif
