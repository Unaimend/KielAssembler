//
// Created by td on 10/7/20.
//

#include "DeBruijnGraph.h"

bool Node::isBalanced() const {
    return inDegree == outDegree
;}

bool Node::isSemiBalanced() const {
    return abs(inDegree - outDegree) == 1;
}

bool DeBruijnGraph::hasEulerianWalk() const {
    return neither == 0 && semi == 2;
}

bool DeBruijnGraph::hasEulerianCycle() const {
    return neither == 0 && semi == 0;
}

bool DeBruijnGraph::isEulerian() const {
    return hasEulerianWalk() or hasEulerianCycle();
}

//TODO Das der letzte Param keine ref sein?
void visit(std::vector<Node>& tour, const Node& n, std::unordered_map<Node, std::vector<Node>>& graph)
{
    while(!graph.at(n).empty())
    {
        Node dst = graph[n].back();
        graph[n].pop_back();
        visit(tour, dst, graph);
    }
    tour.push_back(n);
};

std::optional<std::vector<Node>> DeBruijnGraph::hasEulerianWalkdOrCycle(){

    if (!isEulerian())
    {
        return std::optional<std::vector<Node>>{};
    }
    std::unordered_map<Node, std::vector<Node>> temp = graph;
    if(hasEulerianWalk())
    {
        try{
            // if graph contains l add NodeR as destination
            //TODO find out if g is multimap or just 1 to many
            temp.at(tail).push_back(head);
        } catch (std::out_of_range& exception) {
            temp[tail] = std::vector<Node>{head};
        }
    }
    //TODO CHECK IF CORRECT
    std::vector<Node> tour{};
    const Node& src = temp.begin()->first;

    /*std::function<void(const Node&)> visit;
    visit = [&temp, &visit, &tour](const Node& n) mutable
    {
        while(!temp[n].empty())
        {
            Node dst = temp[n].back();
            temp[n].pop_back();
            visit(dst);
        }
        tour.push_back(n);
    };
     */
    visit(tour, src, temp);

    std::reverse(tour.begin(), tour.end());
    //Remove last element
    tour.erase(tour.end()-1);

    if(hasEulerianWalk())
    {
        auto sti = std::find(tour.begin(), tour.end(), head);
        //TODO USELESS COPY
        std::vector<Node> t;
        for(auto it = sti; it < tour.end(); ++it)
        {
            t.push_back(*it);
        }
        for(auto it = tour.begin(); it < sti; ++it)
        {
            t.push_back(*it);
        }

        return std::optional<std::vector<Node>>{t};
    }
    return std::optional<std::vector<Node>>{tour};

}

DeBruijnGraph::DeBruijnGraph(const std::string &sequenceToAssemble, int kmerLength_) : kmerLength(kmerLength_) {
    //TODO StringViews?
    std::string kmer;
    std::string kmerL;
    std::string kmerR;
    for (int i = 0; i < sequenceToAssemble.length() - (kmerLength - 1); ++i) {
        //(st[i:i+k], st[i:i+k-1], st[i+1:i+k])
        Node* nodeL;
        Node* nodeR;
        kmer = sequenceToAssemble.substr(i, kmerLength);
        kmerL = kmer.substr(0, (kmerLength-1));
        kmerR = kmer.substr(1, kmerLength);

        try
        {
            nodeL = &kmerToNode.at(kmerL);
        } catch (std::out_of_range& exception)
        {
            nodeL = &(kmerToNode[kmerL] = Node(kmerL));
        }

        try
        {
            nodeR = &kmerToNode.at(kmerR);
        } catch (std::out_of_range& exception)
        {
            nodeR = &(kmerToNode[kmerR] = Node(kmerR));
        }

        nodeR->inDegree += 1;
        nodeL->outDegree += 1;

        try{
            // if graph contains l add NodeR as destination
            //TODO find out if g is multimap or just 1 to many
            graph.at(*nodeL).push_back(*nodeR);
        } catch (std::out_of_range& exception) {
            graph[*nodeL] = std::vector<Node>{*nodeR};
        }
    }

    for(const auto& it : kmerToNode)
    {
        const auto& node = it.second;
        if(node.isBalanced())
        {
            balanced += 1;
        }
        else if (node.isSemiBalanced())
        {
            if(node.inDegree == node.outDegree + 1)
            {
                tail = node;
            }
            if(node.inDegree == node.outDegree - 1)
            {
                head = node;
            }
            semi += 1;
        }
        else
        {
            neither += 1;
        }
    }

}
