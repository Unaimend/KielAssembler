//
// Created by Benno Doerr on 10/21/2020.
//

#ifndef KIELASSEMBLER_GRAPH_H	
#define KIELASSEMBLER_GRAPH_H

#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

class DeBruijnGraphAlt {
  public:
    static constexpr size_t NONE = std::numeric_limits<size_t>::max();
    std::string m_sequence;
    std::vector<std::string_view> m_kmer;
    std::vector<std::vector<size_t>> m_edgesIn;
    std::vector<std::vector<size_t>> m_edgesOut;
    std::vector<size_t> m_mergedWith;
    std::vector<bool> m_isActive;
    std::unordered_map<std::string_view, size_t> m_kmerMap;

    size_t head;
    size_t tail;

    bool m_hasEulerianWalk;
    bool m_hasEulerianCycle;

  private:
    size_t find_or_create_node( std::string_view kmer ) {
        size_t index;
        const auto it = m_kmerMap.find( kmer );

        if ( it == m_kmerMap.end() ) {
            index = create_node( kmer );
        } else {
            index = it->second;
        }

        return index;
    }

    size_t create_node( std::string_view kmer ) {
        auto index = m_kmer.size();

        m_kmerMap.emplace( kmer, index );
        m_kmer.emplace_back( std::move( kmer ) );
        m_edgesIn.emplace_back();
        m_edgesOut.emplace_back();
        m_mergedWith.emplace_back( NONE );
        m_isActive.emplace_back( true );

        return index;
    }

    bool is_node_balanced( size_t index ) const {
        return m_edgesIn[index].size() == m_edgesOut[index].size();
    }

    bool is_node_semi_balanced( size_t index ) const {
        return std::abs( static_cast<long>( m_edgesIn[index].size() ) -
                         static_cast<long>( m_edgesOut[index].size() ) ) == 1;
    }

    size_t get_node_degree( size_t index ) {
        return get_node_in_degree( index ) + get_node_out_degree( index );
    }

    size_t get_node_in_degree( size_t index ) {
        return m_edgesIn[index].size();
    }

    size_t get_node_out_degree( size_t index ) {
        return m_edgesOut[index].size();
    }

  public:
    DeBruijnGraphAlt( const std::string &sequenceToAssemble, size_t kmerLength ) : m_sequence( sequenceToAssemble ) {
        for ( size_t i = 0; i < m_sequence.size() - ( kmerLength - 1 ); i++ ) {
            auto kmerL = std::string_view( m_sequence.data() + i, kmerLength - 1 );
            auto kmerR = std::string_view( m_sequence.data() + i + 1, kmerLength - 1 );
            auto iNodeL = find_or_create_node( kmerL );
            auto iNodeR = find_or_create_node( kmerR );

            m_edgesOut[iNodeL].push_back( iNodeR );
            m_edgesIn[iNodeR].push_back( iNodeL );
        }

        // check if graph has an eulerian walk/cycle
        // find head/tail
        {
            // size_t balanced = 0;
            size_t semiBalanced = 0;
            size_t neither = 0;

            for ( size_t i = 0; i < m_kmer.size(); i++ ) {
                if ( is_node_balanced( i ) ) {
                    // balanced++;
                } else if ( is_node_semi_balanced( i ) ) {
                    semiBalanced++;

                    if ( get_node_out_degree( i ) > get_node_in_degree( i ) ) {
                        head = i;
                    } else {
                        tail = i;
                    }
                } else {
                    neither++;
                }
            }

            m_hasEulerianWalk = ( neither == 0 && semiBalanced == 2 );
            m_hasEulerianCycle = ( neither == 0 && semiBalanced == 0 );
        }
    }

    bool is_eulerian() {
        return has_eulerian_cycle() || has_eulerian_walk();
    }

    bool has_eulerian_walk() {
        return m_hasEulerianWalk;
    }

    bool has_eulerian_cycle() {
        return m_hasEulerianCycle;
    }

    std::vector<size_t> get_euler_path() const {
        // stack St;
        std::vector<size_t> stack;
        std::vector<size_t> euler_path;
        auto edges = m_edgesOut;

        // put start vertex in St;
        stack.push_back( 0 );

        // until St is empty
        while ( !stack.empty() ) {

            // let V be the value at the top of St;
            const auto v = stack.back();

            auto &v_out_edges = edges[v];

            // if degree(V) = 0, then % probably meant to be outDegree
            if ( v_out_edges.empty() ) {
                // add V to the answer;
                euler_path.push_back( v );
                // remove V from the top of St;
                stack.pop_back();
            }
            // otherwise
            else {
                // find any edge coming out of V;
                auto tgt = v_out_edges.back();
                // remove it from the graph;
                v_out_edges.pop_back();
                // put the second end of this edge in St;
                stack.push_back( tgt );
            }
        }

        return euler_path;
    }

    void toDot( const std::string &filename ) {
        std::ofstream file;
        file.open( filename, std::ios::out );
        file << "digraph {\n";

        for ( size_t i = 0; i < m_edgesOut.size(); i++ ) {
            for ( const auto &j : m_edgesOut[i] ) {
                file << m_kmer[i] << "->" << m_kmer[j] << "\n";
            }
        }

        file << "}";
        file.close();
    }
};
#endif //KIELASSEMBLER_GRAPH_H 
