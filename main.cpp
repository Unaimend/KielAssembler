#include "deBruijn/DeBruijnGraph.h"
#include "deBruijn/DeBruijnGraphAlt.h"

#include "lib/bioio.hpp"
#include "lib/cxxopts.h"
#include "lib/easylogging++.h"
#include "lib/CLI11.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <sys/resource.h>

INITIALIZE_EASYLOGGINGPP

int main( int argc, char **argv ) {
    CLI::App app{"App description"};
    //
    std::string filepath = "";
    int kmer_length = 0;
    app.add_option("-f,--file", filepath, "Filepath of the fasta file you want to assemble");
    app.add_option("-k,--kmer", kmer_length, "k-mer length with which the graph will be build");
    CLI11_PARSE(app, argc, argv);

    //-------------------------Parsed argument checking-------------
    if( filepath == "")
    {
        LOG(ERROR) << "You must supply a valid path to a fasta file";
        return -1;
    }

    if( kmer_length <= 0)
    {
        LOG(ERROR) << "You must supply a valid kmer length (length >= 0)";
        return -1;
    }
    LOG(INFO) << "File loaded from " + filepath << std::endl;


    //-------------------------FASTA LOADING -------------------------
    //Loades the supplied fasta file and concatenates it to a single string
    std::string text;
    auto load_start = std::chrono::system_clock::now();
    std::ifstream myfile( filepath );
    auto record = bioio::read_fasta( myfile );
    auto load_end = std::chrono::system_clock::now();
    auto append_start = std::chrono::system_clock::now();
    for ( const auto &it : record )
    {
        text.append( it.sequence );
    }

    auto append_end = std::chrono::system_clock::now();
    LOG( INFO ) << "Text appending took: " +
                       std::to_string( std::chrono::duration<double>( append_end - append_start ).count() ) + "seconds";
    LOG( INFO ) << "Text loaded: Length: " + std::to_string( text.length() );

    if(text.length()/kmer_length <= 2)
    {
        LOG(ERROR) << "Sequenced length divided by kmer length must be greater then 2";
        return -1;
    }

    LOG( INFO ) << "Bytes loaded for the text string: " + std::to_string( text.size() );
    LOG( INFO ) << "KiloBytes loaded for the text string: " + std::to_string( text.size() / 1024.0 );
    LOG( INFO ) << "MegaBytes loaded for the text string: " +
                       std::to_string( (float)text.size() / ( 1024.0 * 1024.0 ) );
    LOG( INFO ) << "GigaBytes loaded for the text string: " +
                       std::to_string( (float)text.size() / ( 1024.0 * 1024.0 * 1024.0 ) );

    auto build_start = std::chrono::system_clock::now();
    //Build a DeBruijn Graph from the supplied fast file
    auto graph = DeBruijnGraphAlt(text, kmer_length );
    auto build_end = std::chrono::system_clock::now();
    LOG( INFO ) << "Text building took: " +
                       std::to_string( std::chrono::duration<double>( ( build_end - build_start ) ).count() / ( 60 ) ) +
                       " minutes";

    graph.to_contigs(0);
    graph.toDot("test1");












    /*auto tour_start = std::chrono::system_clock::now();
    auto tour = graph.hasEulerianWalkdOrCycle();
    auto tour_end = std::chrono::system_clock::now();
    LOG( INFO ) << "Tour building took: " +
                   std::to_string( std::chrono::duration<double>( ( tour_start - tour_end) ).count() / ( 60 ) ) +
                   " minutes";
    for(auto& it : tour.value())
    {
        std::cout << it->kmer << std::endl;
    }*/
    graph.toDot("bravo");
    // TODO add tour to_dot
    // TODO find out if g is multimap or just 1 to many
    return 0;
}
