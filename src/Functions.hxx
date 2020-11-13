#ifndef FUNCTION_H
#define FUNCTION_H

#include "TChain.h"
#include <memory>

namespace Function
{
    void Help();

    bool Parsing(int argc, char * argv[]);

    bool InputFile(const std::unique_ptr<TChain>& tree);
}

#endif
