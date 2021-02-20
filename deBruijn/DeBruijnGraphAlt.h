//
// Created by td on 2/5/21.
//

#ifndef INC_4C5222FFA8CA43E1A84139B063DA915D
#define INC_4C5222FFA8CA43E1A84139B063DA915D
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
#include <stack>
class DeBruijnGraphAlt {
  public:
    static constexpr size_t NONE = std::numeric_limits<size_t>::max();
    using ID = size_t;
    std::string m_sequence;
    //Vector off k_mers
    std::vector<std::string_view> m_kmer;
    //All incoming edges for a
    std::vector<std::vector<ID>> m_edgesIn;
    std::vector<std::vector<ID>> m_edgesOut;
    std::vector<bool> m_isActive;
    //Connects k-mer to its ID, the in is just an incremented size_value
    //to which one is added for every new k_mer
    std::unordered_map<std::string_view, ID> m_kmerMap;
    //Only in DEBUG, stores all node ids which are safe to merge, should be added by to_contig
    std::vector<ID> safe_nodes;

    size_t head;
    size_t tail;

    bool m_hasEulerianWalk;
    bool m_hasEulerianCycle;

  private:
    size_t find_or_create_node( std::string_view kmer ) {
        size_t index;
        const auto it = m_kmerMap.find( kmer );

        if ( it == m_kmerMap.end() ) {
            //If no same k-mer is found we have to create a new one
            index = create_node( kmer );
            //and set the index
        } else {
            //Else set the id
            index = it->second;
        }
        //And return it
        return index;
    }

    size_t create_node( std::string_view kmer ) {
        auto index = m_kmer.size();

        m_kmerMap.emplace( kmer, index );
        m_kmer.emplace_back( kmer );
        m_edgesIn.emplace_back();
        m_edgesOut.emplace_back();
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

    void to_contigs(ID start_node)
    {
        std::vector<bool> marked(m_kmerMap.size(), false);
        std::stack<ID> nodes;
        nodes.push(start_node);
        while(!nodes.empty())
        {
            auto cur_node = nodes.top();
            nodes.pop();
            if(!marked[cur_node])
            {
                std::cout << cur_node << std::endl;
                auto outgoing = m_edgesOut[cur_node];
                auto incoming = m_edgesIn[cur_node];
                if((outgoing.size() == 1) && (incoming.size()  <= 1))
                {
                    //merge the current node with the next node
                    //todo this whe have to check if the next node is safe
                    //This is safe since we know that we have exactly one outgoing edge
                    //Todo Extera case for last node
                    auto outgoing_next = m_edgesOut[outgoing[0]];
                    auto incoming_next = m_edgesIn[outgoing[0]];
                    //check if next node is safe
                    //We have to check that our current_node is the only incoming node
                    //and also that the next_node does not have more then one outgoing edge
                    //Zero is also fine because then we have reached an end
                    if(outgoing_next.size() <= 1 && incoming_next.size() == 1)
                    {
                        //both of these will no longer exists

                        //Get both k_mers
                        auto& kmer1 = m_kmer[cur_node];
                        auto& kmer2 = m_kmer[outgoing[0]];
                        //merge kmer2 into kmer1
                        //TODO Make a new strinz view which is just larger
                        //kmer1 += kmer2;
                        //add outgoing edges from kmer2 to kmer1 so the graph stays connected
                        for(const auto it : outgoing_next)
                        {
                            m_edgesOut[cur_node].push_back(it);
                        }
                        //Remove both k_mers from the string_view to id map
                        //And new string_view to the map
                        //Remove all incoming and outgoing edges from k_mer2, or remove it from both vectors not sure atm
                        //We also have to update the incoming_edge from the node after next_node, i.e. add our new node and remove the old one
                        //Delete old k_mers from vector, we have to this since after the merge
                    }
                }
                marked[cur_node] = true;
            }

            for(auto it : m_edgesOut[cur_node])
            {
                if(!marked[it])
                {
                    nodes.push(it);
                }
            }
        }
    }
    void toDot( const std::string &filename ) {
        std::ofstream file;
        file.open( filename, std::ios::out );
        file << "digraph {\n";
        std::cout << "SAFE NODES" <<  safe_nodes.size() << std::endl;
        for ( size_t i = 0; i < m_edgesOut.size(); i++ ) {
            auto k = m_kmer[i];
            if(std::find(safe_nodes.begin(), safe_nodes.end(), i) != safe_nodes.end())
            {
                file << k << "[color = green]" << "\n";

            } else
            {
                file << k << "\n";

            }
        }
        file << '\n';
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

#endif // INC_4C5222FFA8CA43E1A84139B063DA915D
