
#include <iostream>
#include "deBruijn/DeBruijnGraph.h"
#include "lib/cxxopts.h"
#include "lib/bioio.hpp"
#include "fstream"
#include <sys/resource.h>
#include <stdio.h>

int main(int argc, char** argv)
{

    std::string line;
    std::string text;
    std::ifstream myfile ("../data/simulated/ecoli/ecoli1.fna");
    auto record = bioio::read_fasta(myfile , 20000);
    for(const auto& it : record)
    {
        text.append(it.sequence);
    }
    std::cout << "Text loaded" << " " << text.length() <<  std::endl;

    std::string fail = "AGGCCCTGAAGC";
    auto a = DeBruijnGraph(text, 31);
    a.toDot();
    std::cout << "Graph build" << std::endl;
    auto tour = a.hasEulerianWalkdOrCycle();
    //TODO find out if g is multimap or just 1 to many
    return 0;
}
