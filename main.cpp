
#include <iostream>
#include "deBruijn/DeBruijnGraph.h"
#include "lib/cxxopts.h"
#include "lib/bioio.hpp"
#include "fstream"
#include <sys/resource.h>
#include <stdio.h>
#include "lib/easylogging++.h"
#include <chrono>
INITIALIZE_EASYLOGGINGPP


int main(int argc, char** argv)
{

    std::string line;
    std::string text;
    auto load_start = std::chrono::system_clock::now();
    std::ifstream myfile ("../data/simulated/ecoli/ecoli1.fna");
    auto record = bioio::read_fasta(myfile );
    auto load_end= std::chrono::system_clock::now();
    auto append_start = std::chrono::system_clock::now();
    for(const auto& it : record)
    {
        text.append(it.sequence);
    }

    auto append_end = std::chrono::system_clock::now();
    LOG(INFO) << "Text appending took: " + std::to_string(std::chrono::duration<double>(append_end - append_start).count()) + "seconds";
    LOG(INFO) << "Text loaded: Length: " + std::to_string(text.length());
    LOG(INFO) << "Bytes loaded for the text string: " + std::to_string(text.size());
    LOG(INFO) << "KiloBytes loaded for the text string: " + std::to_string(text.size()/1024.0);
    LOG(INFO) << "MegaBytes loaded for the text string: " + std::to_string((float)text.size()/(1024.0*1024.0));
    LOG(INFO) << "GigaBytes loaded for the text string: " + std::to_string((float)text.size()/(1024.0*1024.0*1024.0));

    std::string fail = "AGGCCCTGAAGC";

    auto build_start = std::chrono::system_clock::now();
    auto a = DeBruijnGraph(text, 31);
    auto build_end = std::chrono::system_clock::now();
    LOG(INFO) << "Text building took: " + std::to_string(std::chrono::duration<double>((build_end - build_start)).count()/(60)) + " minutes";
    //a.toDot();
    std::cout << "Graph build" << std::endl;
    auto tour = a.hasEulerianWalkdOrCycle();
    //TODO add tour to_dot
    //TODO find out if g is multimap or just 1 to many
    return 0;
}
