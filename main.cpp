#include "deBruijn/DeBruijnGraph.h"
#include "deBruijn/DeBruijnGraphAlt.h"

#include "lib/bioio.hpp"
#include "lib/cxxopts.h"
#include "lib/easylogging++.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <sys/resource.h>
#include <thread>
#include <vector>
#include <cmath>
#include <cassert>

INITIALIZE_EASYLOGGINGPP

int main( int argc, char **argv ) {

    std::string line;
    std::string text;
    auto load_start = std::chrono::system_clock::now();
    std::ifstream myfile( "/home/td/dev/Bachelorarbeit/data/simulated/ecoli/ecoli1.fna" );
    auto record = bioio::read_fasta( myfile );
    std::cout << record.size() << std::endl;
    auto load_end = std::chrono::system_clock::now();
    auto append_start = std::chrono::system_clock::now();
    for ( const auto &it : record ) {
        text.append( it.sequence );
    }

    auto append_end = std::chrono::system_clock::now();
    LOG( INFO ) << "Text appending took: " +
                       std::to_string( std::chrono::duration<double>( append_end - append_start ).count() ) + "seconds";
    LOG( INFO ) << "Text loaded: Length: " + std::to_string( text.length() );
    LOG( INFO ) << "Bytes loaded for the text string: " + std::to_string( text.size() );
    LOG( INFO ) << "KiloBytes loaded for the text string: " + std::to_string( text.size() / 1024.0 );
    LOG( INFO ) << "MegaBytes loaded for the text string: " +
                       std::to_string( (float)text.size() / ( 1024.0 * 1024.0 ) );
    LOG( INFO ) << "GigaBytes loaded for the text string: " +
                       std::to_string( (float)text.size() / ( 1024.0 * 1024.0 * 1024.0 ) );
    std::string fail = "AGGCCCTGAAGC";
    std::string fail2 = "TAAGCTGATGTT"; // 4 good, 3bad
    std::string fail3 = "ATGCTGTAGCTAGATATCGTAGCTATGCTAGCTAATAGCTATTTCGATGCGGTAGCTAGTGCTAGCATGCGTATGCATGCGTACGGCTAGCTAG"
                        "TAGAGCTCGACTACGACGACGAGAGGGCATCGACGATTAGAGACTAGCGACTACGAGCTAGCGACT";

    //------------------------------------------------------------------------------------------------------------------
    auto build_start = std::chrono::system_clock::now();
    int kmerL = 31;
    unsigned int thread_count = 4;
    std::vector<std::thread> threads;
    int length = std::ceil(text.length()/thread_count);

    assert(length > kmerL);
    std::vector<std::unique_ptr<DeBruijnGraphAlt>> graphs;
    for(int i =0 ; i < text.length()-kmerL; i+= (length-kmerL))
    {
        auto seq = text.substr(i, length );
        graphs.push_back(std::make_unique<DeBruijnGraphAlt>(DeBruijnGraphAlt(seq, kmerL)));
        threads.emplace_back(&DeBruijnGraphAlt::build, graphs.back().get());
    }

    for(auto& it : threads)
    {
        it.join();
    }
    auto build_end = std::chrono::system_clock::now();
    //------------------------------------------------------------------------------------------------------------------
    //auto a = DeBruijnGraphAlt( fail3, 4 );
    LOG( INFO ) << "Text building took: " +
                       std::to_string( std::chrono::duration<double>( ( build_end - build_start ) ).count() / ( 60 ) ) +
                       " minutes";
    // std::cout << "Graph build" << a.kmerToNode.size() << std::endl;
    //auto tour = a.hasEulerianWalkdOrCycle();
    // LOG(INFO) << "HEAD:     " + a.head->kmer ;
    //a.toDot();
    // TODO add tour to_dot
    // TODO find out if g is multimap or just 1 to many
    return 0;
}
