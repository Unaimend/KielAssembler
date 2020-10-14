//
// Created by td on 10/7/20.
//

#ifndef KIELASSEMBLER_DEBRUIJNGRAPH_H
#define KIELASSEMBLER_DEBRUIJNGRAPH_H

#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <functional>
#include <memory>
struct Node
{

    Node() = default;
    explicit Node(std::string kmer_) : kmer(std::move(kmer_))
    {

    }
    //Text/Content of the node //TODO Test vs Char Array
    std::string kmer;
    // Number of incoming edges
    int inDegree{};
    // Number of outgoing edges
    int outDegree{};
    std::vector<Node*> connectedNodes;
    int degree() const
    {
        return inDegree + outDegree;
    }
    //Moving those somewhere else should save space
    [[nodiscard]] bool isSemiBalanced() const;
    [[nodiscard]] bool isBalanced() const;

    bool operator==(const Node& rhs)
    {
        return kmer == rhs.kmer;
    }
    bool operator==(const Node& rhs) const
    {
        return kmer == rhs.kmer;
    }

};


namespace std {
    template<> struct hash<Node>
    {
        std::size_t operator()(const Node & s) const noexcept
        {
            std::size_t h1;
            h1 = std::hash<std::string>{}(s.kmer);
            return h1;
        }
    };
}

class DeBruijnGraph {
public:
    explicit DeBruijnGraph(const std::string& sequenceToAssemble, int kmerLength_ = 30);

    bool hasEulerianWalk() const;

    bool hasEulerianCycle() const;

    bool isEulerian() const;

    std::optional<std::vector<Node>> hasEulerianWalkdOrCycle();

   //Maps Kmers to Node Objects
    std::unordered_map<std::string, std::unique_ptr<Node>> kmerToNode;
    void visit(std::vector<Node>& tour, Node* src);
    Node* head;
    Node* tail;

    int kmerLength;
    int semi = 0;
    int balanced = 0;
    int neither = 0;

};


#endif //KIELASSEMBLER_DEBRUIJNGRAPH_H

