#include "Utils.hxx"

void PrintProgress(int i, int total) {
    if (i == 0)
        std::cout << "Processing..." << std::endl;
    float progress = (double)(i+1)/total;
    const int barWidth = 60 - 8;
    std::cout << "|-[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();

    if (int(progress * 100.0) >= 100) {
        std::cout << std::endl;
        std::cout << "|-done " << std::endl;
    }
}
