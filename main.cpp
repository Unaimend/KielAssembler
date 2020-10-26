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

void do_stuff(std::unique_ptr<DeBruijnGraphAlt>& graph)
{
    std::cout << graph->m_sequence << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(10000));
}


int main( int argc, char **argv ) {

    std::string line;
    std::string text;
    auto load_start = std::chrono::system_clock::now();
    std::ifstream myfile( "../data/simulated/ecoli/ecoli1.fna" );
    auto record = bioio::read_fasta( myfile );
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
    std::string fail3 = "ATGCTGTAGCTAGATATCGTAGCTATGCTAGCTAATHGCTATTTCGATGCGGTAGCTAGTGCTAGCATGCGTATGCATGCGTACGGCTAGCTAG"
                        "TAGAGCTCGACTACGACGACGAGAGGGCATCGACGATTAGAGACTAGCGACTACGAGCTAGCGACT";

    //------------------------------------------------------------------------------------------------------------------
    int kmerL = 4;
    unsigned int thread_count = 32;
    std::vector<std::thread> threads;
    int length = std::ceil(fail3.length()/thread_count);

    assert(length > kmerL);
    std::cout << fail3.substr(0, length) << std::endl;
    std::vector<std::unique_ptr<DeBruijnGraphAlt>> graphs;
    for(int i = length; i < fail3.length()-kmerL; i+= (length-kmerL))
    {
        auto seq = fail3.substr(i, length );
        std::cout << seq << std::endl;
        graphs.emplace_back(std::make_unique<DeBruijnGraphAlt>(DeBruijnGraphAlt(seq, kmerL)));
        //threads.emplace_back(std::thread(do_stuff, std::ref(graphs.back())));
    }

    for(auto& it : threads)
    {
        it.join();
    }
    //------------------------------------------------------------------------------------------------------------------
    auto build_start = std::chrono::system_clock::now();
    auto a = DeBruijnGraphAlt( fail3, 4 );
    auto build_end = std::chrono::system_clock::now();
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
