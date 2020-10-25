//
// Created by td on 10/7/20.
//

#include "DeBruijnGraph.h"
#include <memory>
#include <stack>
#include <list>
#include <fstream>

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


std::optional<DeBruijnGraph::TourType> DeBruijnGraph::hasEulerianWalkdOrCycle(){

    if (!isEulerian())
    {
        return std::optional<TourType>{};
    }
    //BRAUCH MAN DAS
    //std::unordered_map<Node, std::vector<Node>> temp = graph;

    if(hasEulerianWalk())
    {
        //Make it to an eulerian cycle
        try{
            //TODO DAS WIRD NET AKTUALISIERT
            kmerToNode.at(tail->kmer)->connectedNodes.push_back(head);
            ++kmerToNode.at(tail->kmer)->outDegree;
            ++head->inDegree;
        } catch (std::out_of_range& exception) {
            //TODO DAS AUCH NET
            kmerToNode[tail->kmer]->connectedNodes = std::vector<Node*>{head};
            ++kmerToNode.at(tail->kmer)->outDegree;
            ++head->inDegree;
        }
    }
    merge();
    calculate_nodes();

    return std::optional<TourType>{};
    //TODO REBUILD kmerToNode
    /*TourType tour{};
    Node* src =  kmerToNode.begin()->second.get();

    std::stack<Node*> nodes;
    nodes.push(src);
    while(!nodes.empty())
    {
       auto& v = nodes.top();
       if(v->outDegree == 0)
       {
           //Useless copy
           tour.push_back(v);
           nodes.pop();
       }
       else
       {
            //Find outgoing edgr
            auto j = kmerToNode[v->kmer]->connectedNodes.back();
            //Remove it from graph
            kmerToNode[v->kmer]->connectedNodes.pop_back();
            //since we removed a edge we must decrease outgoing edges
            --v->outDegree;
            // TODO we probably have to remove in degree for this node too
            //--kmerToNode[j->kmer]->inDegree;
            nodes.push(j);
       }
    }


    //FIND OUT WHY
    std::reverse(tour.begin(), tour.end());
    //Remove last element
    tour.erase(tour.end()-1);

    if(hasEulerianWalk())
    {
        std::cout << "CALLED" << std::endl;
        //TODO this linke seem dangerours
        auto sti = std::find(tour.begin(), tour.end(), head);
        //TODO USELESS COPY
        TourType t;
        for(auto it = sti; it < tour.end(); ++it)
        {
            t.push_back(*it);
        }
        for(auto it = tour.begin(); it < sti; ++it)
        {
            t.push_back(*it);
        }

        return std::optional<TourType>{t};
    }
    return std::optional<TourType>{tour};
    */

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

        const auto it = kmerToNode.find(kmerL);
        if(it != kmerToNode.end())
        {
            nodeL =  it->second.get();
        }
        else
        {
            //if object with this kmers doenst exist make a one and store a pointer to it in the map
            kmerToNode[kmerL] = std::make_unique<Node>(kmerL);
            //set nodeL to this pointer
            nodeL = kmerToNode[kmerL].get();
        }

        const auto it2 = kmerToNode.find(kmerR);
        if(it2 != kmerToNode.end())
        {
            //check if object with this kmer exists if yes get it'
            nodeR = it2->second.get();
        }
        else
       {
            //if object with this kmers doenst exist make a one and store a pointer to it in the map
            kmerToNode[kmerR] = std::make_unique<Node>(kmerR);
            //set nodeL to this pointer
            nodeR = kmerToNode[kmerR].get();
        }
        nodeR->inDegree += 1;
        nodeL->outDegree += 1;
        nodeL->connectedNodes.push_back(nodeR);
    }

    // Build contigs
    calculate_nodes();

}

void DeBruijnGraph::calculate_nodes()
{
    balanced = 0;
    semi = 0;
    neither = 0;
    for(const auto& it : kmerToNode)
    {
        const auto& node = it.second;
        if(node->isBalanced())
        {
            balanced += 1;
        }
        else if (node->isSemiBalanced())
        {
            if(node->inDegree == node->outDegree + 1)
            {
                tail = node.get();
            }
            if(node->inDegree == node->outDegree - 1)
            {
                head = node.get();
            }
            semi += 1;
        }
        else
        {
            neither += 1;
        }
    }
}

void DeBruijnGraph::toDot() {
    std::ofstream myfile;
    myfile.open ("example.dot", std::ios::out );
    myfile << "digraph {\n";
    for(const auto& it : kmerToNode)
    {
        for(const auto& it2: it.second->connectedNodes)
        {
            myfile << it.second->kmer << "->" << it2->kmer << "\n";
        }

    }
    myfile << "}";
    myfile.close();
}


void DeBruijnGraph::dfs(Node* src_)
{
    src_->marked = true;
    auto beg = 0;
    auto end = src_->connectedNodes.size();
    while(beg < end)
    {
        //dont visited node twice
        if(!(src_->connectedNodes[beg]->marked))
        {
            //TO STUFF HERE
            //Safe to merge
            if((src_->connectedNodes[beg])->inDegree == 1 && src_->outDegree == 1) {
                //std::cout << beg << "CALL    " << src_->kmer  << std::endl;
                //if outdegree > 1 we need to fix sth
                //if ingegree > 1 we would fuck ptr
                //add kmer so we dont loose information
                //TODO MERGE CORRECTLY MUSS ICHC NIH NUR DEN LETZTEN RANPACKEN
                char test = src_->connectedNodes[beg]->kmer.back();
                src_->kmer.push_back(test);
                //Add outgoig nodes from the src to our outgoing nodese
                for (auto x : (src_->connectedNodes[beg])->connectedNodes) {
                    src_->connectedNodes.push_back(x);
                    src_->outDegree += 1;
                }

                //remove the currentcly watched at vertive, which muste be the beginn
                //since outDegree == 1
                //TODO IS POTENTIALLY DANGEROUS? Danglin ptr? da unique_ptr speicher freigibt
                kmerToNode.erase((src_->connectedNodes[beg])->kmer);
                src_->connectedNodes.erase(src_->connectedNodes.begin());
                //haben ja einen weg gemacht
                src_->outDegree--;
                //irgendwas in die richtung muss sein
                if(src_->outDegree == 1)
                {
                    //WR MUEDE KAOB RICHTIG, wenn wir uns nooch mal treffen
                    // und immernoch out==1 haben adden
                    src_->marked = 0;
                    //soll fail2 fixen
                }
                dfs(src_);
            }
            else
            {
                dfs((src_->connectedNodes[beg]));
            }
        }
        ++beg;
    }
}

void DeBruijnGraph::merge()
{
    std::cout << "HEAD" << head->kmer << std::endl;
    dfs(head);
}
