
#include <iostream>
#include "deBruijn/DeBruijnGraph.h"
int main(int argc, char** argv)
{
    std::cout << "Hello, World!" << std::endl;

    char test[] = "t";
    auto a = DeBruijnGraph("AGGCCCTGAAGC", 4);
    /*for(const auto& it : a.kmerToNode)
    {
        std::cout << "T  " << it.first << it.second.kmer << std::endl;
    }*/
    for(const auto& it : a.graph)
    {
        std::cout << "SRC " << it.first.kmer << std::endl;
        for(auto&  edge : it.second)
        {
            std::cout <<   "DEST:" <<  edge.kmer << std::endl;
        }

    }
    auto pls = a.hasEulerianWalkdOrCycle();
    std::cout << "---------------------------------------" << std::endl;
    for(const auto& it : pls.value())
    {
        std::cout << it.kmer << std::endl;
    }
    //TODO find out if g is multimap or just 1 to many
    return 0;
}
