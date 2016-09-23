//
// File: Family.cpp
// Created by: Thomas Bigot
//

/*
 *   Copyright or © or Copr. Thomas Bigot 2012
 * 
 *   This software is a computer program whose purpose is to provide classes
 *   for phylogenetic data analysis.
 * 
 *   This software is governed by the CeCILL  license under French law and
 *   abiding by the rules of distribution of free software.  You can  use,
 *   modify and/ or redistribute the software under the terms of the CeCILL
 *   license as circulated by CEA, CNRS and INRIA at the following URL
 *   "http://www.cecill.info".
 * 
 *   As a counterpart to the access to the source code and  rights to copy,
 *   modify and redistribute granted by the license, users are provided only
 *   with a limited warranty  and the software's author,  the holder of the
 *   economic rights,  and the successive licensors  have only  limited
 *   liability.
 * 
 *   In this respect, the user's attention is drawn to the risks associated
 *   with loading,  using,  modifying and/or developing or reproducing the
 *   software by the user in light of its specific status of free software,
 *   that may mean  that it is complicated to manipulate,  and  that  also
 *   therefore means  that it is reserved for developers  and  experienced
 *   professionals having in-depth computer knowledge. Users are therefore
 *   encouraged to load and test the software's suitability as regards their
 *   requirements in conditions enabling the security of their systems and/or
 *   data to be ensured and,  more generally, to use and operate it in the
 *   same conditions as regards security.
 * 
 *   The fact that you are presently reading this means that you have had
 *   knowledge of the CeCILL license and that you accept its terms.
 */


#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <numeric>
#include <ios>
#include <math.h>

#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <boost/iterator/iterator_concepts.hpp>

#include "Family.hpp"
#include "Waiter.hpp"
#include "TreeTools.hpp"
#include "Taxon.hpp"
#include "DataBase.hpp"
#include "CandidateNode.hpp"
#include "Pattern.hpp"


using namespace std;
using namespace bpp;
using namespace tpms;


namespace tpms{
    Family::Family(stringstream* sIntro, string* sNewick, DataBase* dbp): db_(dbp), preamble_(sIntro), newick_(sNewick), containsUndefinedSequences_(false), doneMapping_NodeToTaxa_(false) {
        
    }
    
    
    void Family::initialize(){
        // extraction du nom de la famille
        string currLigne;
        getline(*preamble_,currLigne);
        istringstream sFamilleName(currLigne);
        sFamilleName >> name_;
        // extraction des synonymes
        getline(*preamble_, currLigne);
        if(currLigne.at(0) != '[') std::cerr << "Bracket (newick comments beginning) expected at family " << name_ << endl;
        
        // cas du fichier avec deux arbres en plus (concerne les arbres réconciliés
        getline(*preamble_, currLigne);
        if(currLigne.at(0) == '#') { getline(*preamble_, currLigne);getline(*preamble_, currLigne);getline(*preamble_, currLigne);} // on saute ces lignes
        
        bool mneFini;
        bool specFini;
        
        
        char lecar;
        
        
        while(currLigne.at(0) != ']') {
            stringstream espece;
            ostringstream mnemonique;
            mnemonique.str("");
            espece.str("");
            espece.seekg(0);
            mneFini = false;
            specFini = false;
            for(unsigned int i = 0; !specFini && i <currLigne.size(); i++) {
                lecar = currLigne.at(i);
                if(lecar == '"') mneFini=true;
                else if (lecar == ':') specFini=true;
                else {
                    if(mneFini) espece << lecar;
                    else mnemonique << lecar;
                }
                
            }
            
            // removing commas (that are forbidden in species names)
            ostringstream cleanSpeciesName;
            string spcPart;
            while(getline(espece,spcPart,',')){
                cleanSpeciesName << spcPart;
            }
            
            
            Taxon * currTaxon = db_->nameToTaxon(cleanSpeciesName.str());
            if(currTaxon != 00){
                mne2tax_.insert(pair<string,Taxon *>(mnemonique.str(),currTaxon ));
                taxa_.insert(currTaxon);
            } else {
                cout << "o Unable to find this species in the species tree:" << cleanSpeciesName.str() << endl;
            }
            
            getline(*preamble_, currLigne);
        }
        
        // on fabrique l'arbre
        //istringstream ssNewick(sNewick);
        //Newick newickReader(false);
        //tree = newickReader.read(ssNewick);
        
        tree_ = tpms::TreeTools::newickToTree(*newick_,false);
        
        /*if(tree->getNumberOfNodes() != tree2->getNumberOfNodes()) cout << "Incoherence famille " << _name << " bpp=" << tree->getNumberOfNodes() << " tpms=" <<tree2->getNumberOfNodes() << endl  ;
         */
        
        delete newick_;
        delete preamble_;
        
        highestID_ = tree_->getNumberOfNodes();
        
    }
    
    void Family::doMapping_LeavesToSpecies(){
        vector<Node *> leaves = tree_->getLeaves();
        this->leaves_.clear();
        this->leaves_.insert(leaves.begin(),leaves.end());
        mapping_NodesToTaxa_.clear();
        mapping_NodesToTaxa_.resize(tree_->getNumberOfNodes(),00);
        for(vector<Node *>::iterator leave = leaves.begin(); leave != leaves.end(); leave++){
            std::map<string,Taxon*>::iterator found = mne2tax_.find((*leave)->getName());
            if(found==mne2tax_.end()){
                containsUndefinedSequences_=true;
                cout << "Warning: The sequence " << (*leave)->getName() << " has not been associated to a species in family " << name_ << ". This familly will be ignored." << endl;
            }
            mapping_NodesToTaxa_.at((*leave)->getId())=found->second;
        }
    }
    
    void Family::doMapping_NodesToNatures(){
        vector<Node *> nodes = tree_->getNodes();
        mapping_NodesToNatures_.resize(tree_->getNumberOfNodes());
        for(vector<Node *>::iterator node = nodes.begin(); node != nodes.end(); node++){
            if((*node)->hasName() && (*node)->getName().at(0) == '#')
                mapping_NodesToNatures_.at((*node)->getId())=DUPLICATION;
            else mapping_NodesToNatures_.at((*node)->getId())=SPECIATION;
        }
    }
    
    
    void Family::doMapping_NodesToUnicityScores() {
        mapping_NodesToUnicityScores_.resize(tree_->getNumberOfNodes());
        compute_UnicityScoreOnNode_(mapping_NodesToUnicityScores_,tree_->getRootNode(),00); // 00 because no origin, we are at the real root
    }
    
    
    void Family::doRerooting_LessTransfers(){
        if(tree_->isRooted())
            tree_->unroot();
        
        // pass 1: best taxonomy
        vector<Node*> bestTaxoRoots = getTaxonomyBestRoots_(tree_->getNodes());
        
        vector<unsigned int> bestOGs = getLessTransfersBestRoots_(tpms::TreeTools::toNodesIDs(bestTaxoRoots));
        reRootAt_(bestOGs);
    }
    
    void Family::doRerooting_Daubin(){
        if(tree_->isRooted())
            tree_->unroot();
        vector<unsigned int> nodesIDofTheTree;
        vector<Node*> nodesOfTheTree = tree_->getNodes();
        
        for(vector<Node*>::iterator currNode = nodesOfTheTree.begin(); currNode != nodesOfTheTree.end(); currNode++)
            nodesIDofTheTree.push_back((*currNode)->getId());
        
        vector<unsigned int> bestOGs = getDaubinCriteriaBestRoots_(nodesIDofTheTree);
        reRootAt_(bestOGs);
    }
    
    void Family::reRootAt_(vector<unsigned int> bestOGs){
        vector<Node*> nodes;
        for(vector<unsigned int>::iterator currNodeID = bestOGs.begin(); currNodeID != bestOGs.end(); currNodeID++){
            nodes.push_back(tree_->getNode(*currNodeID));
        }
        reRootAt_(nodes);
    }
    
    
    vector<unsigned int> Family::getLessTransfersBestRoots_(vector<unsigned int> nodes){
        double bestGainPerTransfer;
        vector<unsigned int> bestOutgroups;
        
        
        TreeTemplate<Node> * savedTree = tree_->clone();
        
        for(vector<unsigned int>::iterator currOG = nodes.begin(); currOG != nodes.end(); currOG++){
            tree_->newOutGroup(tree_->getNode(*currOG));
            highestID_ = tree_->getNumberOfNodes();
            doMapping_NodesToTaxa();
            unsigned int gainSum = 0;
            unsigned int currNbOfTransfers = compute_detectTransfers_(false,gainSum);
            double gainPerTransfer = double(gainSum) / double(currNbOfTransfers);
            
            if(bestOutgroups.empty() || gainPerTransfer > bestGainPerTransfer){
                bestOutgroups.clear();
                bestOutgroups.push_back(*currOG);
                bestGainPerTransfer = gainPerTransfer;
            } else if(gainPerTransfer == bestGainPerTransfer) {
                bestOutgroups.push_back(*currOG);
            }
            delete(tree_);
            tree_ = savedTree->clone();
            doMapping_LeavesToSpecies();
        }
        delete savedTree;
        
        return(bestOutgroups);
    }
    
    vector<Node*> Family::getUnicityBestRoots_(vector<Node *> nodes){
        // aim of this function: trying all branches as roots, and return those which minimize the unicity score sum.
        vector<float> currScores;
        vector<Node*> bestRoots; // returned at the end
        float bestScoresSum;
        float currScoresSum;
        currScores.resize(tree_->getNumberOfNodes()+1);
        // we root at the branches leading to the node
        for(vector<Node *>::iterator currNode = nodes.begin(); currNode != nodes.end(); currNode++) {
            if(!(*currNode)->hasFather()) continue; // root has not father branch
            compute_UnicityScoreOnNode_(currScores,*currNode,00,true);
            currScoresSum = accumulate(currScores.begin(),currScores.end(),0);
            if(tree_->getRootNode()->getNumberOfSons()==2)
                currScoresSum -= currScores.at(tree_->getRootId());
            if(bestRoots.empty() || currScoresSum <  bestScoresSum){
                // we've found a first or better root candidate
                bestRoots.clear();
                bestRoots.push_back(*currNode);
                bestScoresSum = currScoresSum;
            } else if (currScoresSum == bestScoresSum){
                bestRoots.push_back(*currNode);
            }
        }
        return(bestRoots);
    }
    
    
    void Family::doRerooting_Unicity() {
        // aim of this function: choose the best branch to be a root according to minimum unicity scores sum
        // and re-root the tree according to this branch
        vector<Node *> nodes = tree_->getNodes();
        vector<Node*> bestRoots = getUnicityBestRoots_(nodes);
        reRootAt_(bestRoots);
        doMapping_NodesToUnicityScores();
    }
    
    unsigned int Family::getTaxonomicSum_(Node* node, Node* father, vector<Taxon*>* local_nodesToTaxa, set<Node*>* nodesToIgnore)
    {
        unsigned int sum = 0;
        vector<Node*> neighbors = node->getNeighbors();
        for(vector<Node*>::iterator n = neighbors.begin(); n != neighbors.end(); n++){
            if(*n == father || nodesToIgnore->find(*n) != nodesToIgnore->end()) continue;
            sum += getTaxonomicSum_(*n,node,local_nodesToTaxa,nodesToIgnore);
        }
        
        if(!node->isLeaf())
            sum += (local_nodesToTaxa->at(node->getId()) != 00? local_nodesToTaxa->at(node->getId())->getDepth() :0);
        return(sum);
    }
    
    
    
    vector<Node*> Family::getTaxonomyBestRoots_(vector<Node *> nodes){
        vector<Node*> bestRoots; // returned at the end
        float currDepthSum;
        float highestDepthSum = 0;
        // aim of this function: trying all branches as roots, and return those which maximize the taxo level sum
        // we’ll try to root at every branch and see taxonomic mapping
        for(vector<Node *>::iterator currPotentialRoot = nodes.begin(); currPotentialRoot != nodes.end(); currPotentialRoot++) {
            if(!(*currPotentialRoot)->hasFather()) continue; // there is no branch leading to the root node
            
            // considering the branch leading to currPotentialRoot as a root:
            vector<Taxon*> localMapping_NodesToTaxa(highestID_);
            mapNodeOnTaxon_(&localMapping_NodesToTaxa,&localMapping_NodesToTaxa,*currPotentialRoot,(*currPotentialRoot)->getFather());
            mapNodeOnTaxon_(&localMapping_NodesToTaxa,&localMapping_NodesToTaxa,(*currPotentialRoot)->getFather(),*currPotentialRoot);
            
            vector<Taxon*> localMapping_withoutGF(highestID_);
            vector<unsigned int> localMapping_perturbations(highestID_);
            set<Node*> localMapping_nodesInducingPerturbations;
            
            //computeTaxonomicShift(*currPotentialRoot,(*currPotentialRoot)->getFather(),00,00,&localMapping_NodesToTaxa,&localMapping_withoutGF,&localMapping_perturbations,&localMapping_nodesInducingPerturbations);
            //computeTaxonomicShift((*currPotentialRoot)->getFather(),*currPotentialRoot,00,00,&localMapping_NodesToTaxa,&localMapping_withoutGF,&localMapping_perturbations,&localMapping_nodesInducingPerturbations);
            
            mapNodeOnTaxon_(&localMapping_NodesToTaxa,&localMapping_NodesToTaxa,*currPotentialRoot,(*currPotentialRoot)->getFather(),&localMapping_nodesInducingPerturbations);
            mapNodeOnTaxon_(&localMapping_NodesToTaxa,&localMapping_NodesToTaxa,(*currPotentialRoot)->getFather(),*currPotentialRoot,&localMapping_nodesInducingPerturbations);
            
            currDepthSum = getTaxonomicSum_(*currPotentialRoot,(*currPotentialRoot)->getFather(),&localMapping_NodesToTaxa,&localMapping_nodesInducingPerturbations);
            currDepthSum += getTaxonomicSum_((*currPotentialRoot)->getFather(),*currPotentialRoot,&localMapping_NodesToTaxa,&localMapping_nodesInducingPerturbations);
            
            // putting Taxa into 
            
            
            // 	// DEBUG Display info
            // 	cout << "\n\nTOTAL= " << currDepthSum << ". Nodes inducing perturbations: " << localMapping_nodesInducingPerturbations.size() << endl;
            // 	//cout << "Br len :" << (*currPotentialRoot)->getDistanceToFather() <<  endl;
            // 	
            // 	    
            // 	
            // 	//DEBUG Print the tree
            // 	  print_tax_tree(*(*currPotentialRoot)->getFather(),0,*currPotentialRoot,&localMapping_NodesToTaxa,&localMapping_nodesInducingPerturbations,false);
            // 	  print_tax_tree(**currPotentialRoot,0,(*currPotentialRoot)->getFather(),&localMapping_NodesToTaxa,&localMapping_nodesInducingPerturbations,false);
            // 	    
            //removing the depth of the current root if rooted (it is twice in mapping_NodesToTaxa)
            //if(tree->isRooted())
            //    currDepthSum -= localMapping_NodesToTaxa.at(tree->getRootId())->getDepth();
            
            // storing the root if it’s better (replacing) or equal (adding)
            if(bestRoots.empty() || currDepthSum > highestDepthSum){
                // we've found a first or better root candidate
                bestRoots.clear();
                bestRoots.push_back(*currPotentialRoot);
                highestDepthSum = currDepthSum;
            } else if(currDepthSum == highestDepthSum)
                bestRoots.push_back(*currPotentialRoot);
        }
        
        return(bestRoots);
    }
    
    
    void Family::print_tax_tree_(Node& node, unsigned int depth, Node* origin, std::vector<tpms::Taxon*>* mapping, std::set<bpp::Node*>* ignoredNodes, bool subtreeIgnored){
        bool ignoredBegin=false;
        if(ignoredNodes->find(&node)!= ignoredNodes->end()) {subtreeIgnored=true; ignoredBegin=true;}
        for(unsigned int i=0 ; i < depth; i++) cout << " ";
        cout << (subtreeIgnored? (ignoredBegin? "***": "---") : "") << ((node.isLeaf() || mapping->at(node.getId())== 00 ? 0: mapping->at(node.getId())->getDepth())) << "|" << (node.isLeaf()? "*" : " " )<< (
            node.isLeaf() ? (mapping_NodesToTaxa_.at(node.getId()) == 0? "+++++" :mapping_NodesToTaxa_.at(node.getId())->getName()): (
                (mapping->at(node.getId()) == 00)? "": mapping->at(node.getId())->getName()));
        if(node.isLeaf()) cout << " (" << node.getName() << ")";
        cout << endl;
        
        vector<Node*> neighbors = node.getNeighbors();
        for(vector<Node*>::iterator n = neighbors.begin(); n != neighbors.end(); n++){
            if(*n == origin) continue;
            print_tax_tree_(**n,depth+1,&node,mapping,ignoredNodes,subtreeIgnored);
        }
    }
    
    void Family::doRerooting_Taxonomy() {
        vector<Node *> nodes = tree_->getNodes();
        vector<Node*> bestRoots = getTaxonomyBestRoots_(nodes);
        reRootAt_(bestRoots);
    }
    
    
    void Family::doRerooting_UnicityTaxonomy(){
        if(tree_->isRooted())
            tree_->unroot();
        vector<Node *> nodes = tree_->getNodes();
        
        // first using unicity criteria
        vector<Node*> bestRoots = getUnicityBestRoots_(nodes);
        // then, on ex-aequos, using taxonomi criteria
        bestRoots = getTaxonomyBestRoots_(bestRoots);
        reRootAt_(bestRoots);
    }
    
    
    void Family::reRootAt_(std::vector<Node*> bestRoots)
    {
        Node* bestRoot = *bestRoots.begin();
        double maxBranchSize = bestRoot->getDistanceToFather();
        for(vector<Node*>::iterator currRoot = bestRoots.begin()+1; currRoot != bestRoots.end(); currRoot++){
            if((*currRoot)->getDistanceToFather() > maxBranchSize){
                maxBranchSize = (*currRoot)->getDistanceToFather();
                bestRoot = *currRoot;
            }
        }
        
        // "rerooting" the tree according to this new outgroup :
        tree_->newOutGroup(bestRoot);
        //the number of nodes has potentially changed
        highestID_ = tree_->getNumberOfNodes();
        // the present taxa affectation changed: updating
        doMapping_NodesToTaxa();
    }
    
    
    
    map<Taxon*, unsigned int> Family::compute_UnicityScoreOnNode_(vector<float> &scores, Node * node, Node * origin, bool virtualRootOnTheBranch){
        unsigned int id = node->getId();
        map<Taxon*, unsigned int> thisNodeCount;
        
        // step 1: filling the species counts from this node
        vector<Node *> neighbors;
        neighbors = node->getNeighbors();
        map<Taxon*,unsigned int> currNeighborCount;
        if(!virtualRootOnTheBranch)
            for(vector<Node *>::iterator currNeighbor = neighbors.begin(); currNeighbor < neighbors.end(); currNeighbor ++) {
                if(*currNeighbor == origin) continue;
                currNeighborCount = compute_UnicityScoreOnNode_(scores,*currNeighbor,node,false);
                // adding this new map to the current
                for(map<Taxon*,unsigned int>::iterator currCount = currNeighborCount.begin(); currCount != currNeighborCount.end(); currCount++){	    
                    thisNodeCount[currCount->first] += currCount->second;
                }
            }
            else {
                currNeighborCount = compute_UnicityScoreOnNode_(scores,node,node->getFather(),false);
                for(map<Taxon*,unsigned int>::iterator currCount = currNeighborCount.begin(); currCount != currNeighborCount.end(); currCount++){	    
                    thisNodeCount[currCount->first] += currCount->second;
                }
                currNeighborCount = compute_UnicityScoreOnNode_(scores,node->getFather(),node,false);
                for(map<Taxon*,unsigned int>::iterator currCount = currNeighborCount.begin(); currCount != currNeighborCount.end(); currCount++){	    
                    thisNodeCount[currCount->first] += currCount->second;
                }
                
            }
            
            
            if(neighbors.size() == 1){ // leave case
                thisNodeCount.insert(pair<Taxon*,unsigned int>(mapping_NodesToTaxa_.at(id),1));
                
            }
            
            // step 2: Score computation on this node
            float score = 0;
            
            for(map<Taxon*,unsigned int>::iterator currCount = thisNodeCount.begin(); currCount != thisNodeCount.end(); currCount++){
                score += log(float(currCount->second));
            }    
            
            if(virtualRootOnTheBranch)
                *scores.rbegin() = score;
            else
                scores.at(id) = score;
            
            return(thisNodeCount);
    }
    
    void Family::writeTreeToStream(Node* root, ostream & sortie, unsigned int deep){
        for(unsigned int i=0; i< deep; i++) sortie << " ";
        if(root->hasName()) sortie << root->getName();
        else sortie << "<NONAME>";
        sortie << endl;
        for(unsigned int i = 0; i < root->getNumberOfSons() ; i++)
            writeTreeToStream(root->getSon(i), sortie, deep+1);
    }
    
    
    int Family::numberOfNodesBetween(Node * ancestor, Node * pnode){
        if(ancestor == pnode) return 0;
        else return(1+numberOfNodesBetween(ancestor, pnode->getFather()));
    }
    
    
    
    void Family::deleteFromLeavesToBif_(Node * pnode){
        if(pnode->hasFather()){
            Node * father = pnode->getFather();
            int nbFils = father->getNumberOfSons();
            if(nbFils == 1 && father->hasFather())
                deleteFromLeavesToBif_(father);
            else
                father->removeSon(pnode);
            delete pnode;
        } else cout << "deleteFromLeavesToBif : On essaye de supprimer la racine !!!" << endl;
    }
    
    
    Node * Family::removeUniqueSons_(Node * localRoot){
        // FONCTION RÉCURSIVE
        unsigned int nbFils = localRoot->getNumberOfSons();
        switch(nbFils){
            case 0: // Cas de base, on est sur une feuille, donc on retourne la feuille
                return(localRoot);
            case 1: // Cas récursif 1, le nœud a un seul fils
                return(removeUniqueSons_(localRoot->getSon(0)));
            default: // ce nœud n'a pas un fils unique, donc on va attribuer la prochaine bifurcation ou feuille aux fils de localRoot
                for(unsigned currSonIndex=0 ; currSonIndex < nbFils; currSonIndex++)
                    localRoot->setSon(currSonIndex, removeUniqueSons_(localRoot->getSon(currSonIndex)));
                return(localRoot);
        }
    }
    
    TreeTemplate<Node> * Family::getTree() {
        return(tree_);
    }
    
    bool Family::containsSpecie(Taxon* taxon) {
        return(taxa_.find(taxon)!=taxa_.end());
    }
    
    
    Taxon* Family::getTaxonOfNode(Node* node){
        return(mapping_NodesToTaxa_.at(node->getId()));
    }
    
    
    
    void Family::getLeavesFromNode(Node* pnode, std::vector< Node* >& leaves){
        if(pnode->isLeaf()) leaves.push_back(pnode);
        else
            for(unsigned int sonNr = 0; sonNr < pnode->getNumberOfSons() ; sonNr ++)
                getLeavesFromNode(pnode->getSon(sonNr),leaves);
    }
    
    
    void Family::getLeavesFromNode(Node* pnode, std::set< Node* >& leaves, int& leavesNumber){
        if(pnode->isLeaf()) {
            leaves.insert(pnode);
            leavesNumber++;
        }
        else
            for(unsigned int sonNr = 0; sonNr < pnode->getNumberOfSons() ; sonNr ++)
                getLeavesFromNode(pnode->getSon(sonNr),leaves, leavesNumber);
    }
    
    void Family::getLeavesFromNode(Node * pnode, std::set< Node* > & leaves){
        vector<Node *> tempVect;
        getLeavesFromNode(pnode, tempVect);
        // transformation en set
        for(vector<Node *>::iterator currNode = tempVect.begin(); currNode != tempVect.end(); currNode++)
            leaves.insert(*currNode);
    }
    
    set<string> Family::getLeavesNamesFromNode(Node * pnode){
        vector<Node *> leaves;
        getLeavesFromNode(pnode,leaves);
        set<string> result;
        for(vector<Node *>::iterator nit=leaves.begin(); nit != leaves.end(); nit++)
            result.insert((*nit)->getName());
        
        return(result);
    }
    
    
    std::string Family::getName(){
        return(name_);
    }
    
    unsigned int Family::getTaxonomicSum_(std::vector<Taxon*> &taxa){
        unsigned int sum = 0;
        for(std::vector<Taxon*>::iterator taxon = taxa.begin(); taxon != taxa.end(); taxon++){
            if (*taxon == 00)
                continue;
            sum += (*taxon)->getRelativeDepth(this->taxa_);
        }
        return sum;
    }
    
    Taxon* Family::mapNodeOnTaxon_(std::vector<Taxon*>* referenceMapping, vector<Taxon*>* recordResult,bpp::Node* node, bpp::Node* origin,set<Node*>* ignoredNodes,bool recursive, Node* ignoredNode){
        unsigned int currNodeID = node->getId();
        vector<Node*> neighbors;
        neighbors = node->getNeighbors();
        if(leaves_.find(node) != leaves_.end()) // BASE CASE: leaf
            return(referenceMapping->at(currNodeID));
        
        if(!recursive) // other base case: not recursive
            return(referenceMapping->at(currNodeID));
        
        // dealing with the case: topological leaf but not a real leaf (removed subtree)
        if(neighbors.size() == 1){
            //DEBUG
            cout << "\ntopological leaf but not a real leaf (removed subtree)" << endl;
            if(recordResult != 00) recordResult->at(currNodeID)=00;
            return 00;
        }
            
        set<Taxon*> virtualSonsTaxa;
        
        // Case of ignored node.
        // if the ignored node is in the sons no need to keep the funciton recursive: nothing has changed
        
        if(ignoredNode!=00 &&  find(neighbors.begin(),neighbors.end(),ignoredNode) != neighbors.end())
            recursive = false;
        
        
        //collecting taxa on sons (neighbors without the node "origin")
            
        for(vector<Node *>::iterator currNeighbor = neighbors.begin(); currNeighbor != neighbors.end(); currNeighbor++){
            if(*currNeighbor == origin || *currNeighbor == ignoredNode || (ignoredNodes != 00 && ignoredNodes->find(*currNeighbor) != ignoredNodes->end())) continue;
            virtualSonsTaxa.insert(mapNodeOnTaxon_(referenceMapping,recordResult,*currNeighbor,node,ignoredNodes,recursive,ignoredNode));
        }
        
        
        //removing null Taxa from list
        set<Taxon*>::iterator nullFound = virtualSonsTaxa.find(00);
        if(nullFound != virtualSonsTaxa.end()) virtualSonsTaxa.erase(nullFound);
        
        Taxon* currTaxon;
        if(!virtualSonsTaxa.empty()) currTaxon = Taxon::findLCA(virtualSonsTaxa);
        else currTaxon = 00;
        
        
        if(recordResult != 00) recordResult->at(currNodeID) = currTaxon ;
        return(currTaxon);
                
    }
    
    
    void Family::addSequencesNames(Node * currNode)
    {
        Node * seqNode = tree_->getNode(currNode->getId());
        string oldname;
        if(currNode->hasName()) oldname = currNode->getName();
        if(seqNode->hasName()&& currNode->isLeaf()) currNode->setName(oldname + "=" + seqNode->getName());
        else currNode->deleteName();
        vector<Node*> sons = currNode->getSons();
        for(vector<Node *>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++)
            addSequencesNames(*currSon);
    }
    
    
    void Family::labelWithSequencesNames(Node * currNode)
    {
        Node * seqNode = tree_->getNode(currNode->getId());
        if(seqNode->hasName()&& currNode->isLeaf()) currNode->setName(seqNode->getName());
        else currNode->deleteName();
        vector<Node*> sons = currNode->getSons();
        for(vector<Node *>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++)
            labelWithSequencesNames(*currSon);
    }
    
    std::vector<float>& Family::getUnicityScores()
    {
        return(mapping_NodesToUnicityScores_);
    }
    
    
    NodeNature Family::getNatureOf(Node* node)
    {
        return(mapping_NodesToNatures_.at(node->getId()));
    }
    
    
    std::set<Taxon*> & Family::getTaxaOnThisSubtree(Node* node)
    {
        if(cacheMapping_SubtreesToTaxalist_.at(node->getId()).empty()){
            vector<Node *> sons = node->getSons();
            for(vector<Node*>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++){
                set<Taxon*> listOnCurrSon = getTaxaOnThisSubtree(*currSon);
                cacheMapping_SubtreesToTaxalist_.at(node->getId()).insert(listOnCurrSon.begin(),listOnCurrSon.end());
            }
            if(sons.size() == 0) // leave case, return the species of this leave
                cacheMapping_SubtreesToTaxalist_.at(node->getId()).insert(mapping_NodesToTaxa_.at(node->getId()));
        }
        return(cacheMapping_SubtreesToTaxalist_.at(node->getId()));
    }
    
    unsigned int Family::getMaxDepthOfSubtree(Node* node){
        return(mapping_NodesToMaxDepths_.at(node->getId()));
    }
    
    void tpms::Family::clearCache()
    {
        cacheMapping_SubtreesToTaxalist_.resize(0);
    }
    
    void tpms::Family::initCache(){
        cacheMapping_SubtreesToTaxalist_.resize(highestID_);
    }
    
    
    set<Taxon *> &Family::getSpecies(){
        return taxa_;   
    }
    
    
    void Family::doMapping_NodesToTaxa(){
        mapping_NodesToTaxa_.resize(highestID_);
        mapNodeOnTaxon_(&mapping_NodesToTaxa_,&mapping_NodesToTaxa_,tree_->getRootNode());
        updateTaxa_();
        doneMapping_NodeToTaxa_ = true;
    }
    
    void Family::doMapping_NodesToMaxDepth(){
        mapping_NodesToMaxDepths_.resize(highestID_);
        mapNodeToMaxDepth_(tree_->getRootNode());
    }
    
    void Family::updateTaxa_(){
        taxa_.clear();
        for(vector<Taxon*>::iterator currTax = mapping_NodesToTaxa_.begin(); currTax != mapping_NodesToTaxa_.end(); currTax++){
            taxa_.insert(*currTax);}
    }
    
    
    void Family::doMapping_NodesToTaxonomicShift(){
        // (re)initializing the vector
        mapping_NodesToTaxonomicShift_.resize(0);
        mapping_NodesToTaxonomicShift_.resize(highestID_,0);
        
        computed_nodesInducingPerturbation_.clear();
        
        mapping_grandFatherWithoutThisNode_.resize(0);
        mapping_grandFatherWithoutThisNode_.resize(highestID_,0);
        
        mapNodeOnTaxon_(&mapping_NodesToTaxa_,&mapping_NodesToTaxa_,tree_->getRootNode());
        
        computeTaxonomicShift_(tree_->getRootNode(), 00, 00, 00, &mapping_NodesToTaxa_, &mapping_grandFatherWithoutThisNode_, &mapping_NodesToTaxonomicShift_, &computed_nodesInducingPerturbation_);
    }
    
    
    unsigned int Family::mapNodeToMaxDepth_(Node* node){
        unsigned int maxDepth = 0;
        if(!node->isLeaf()){
            vector<Node*> sons = node->getSons();
            for(vector<Node*>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++){
                unsigned int currMaxDepth = mapNodeToMaxDepth_(*currSon);
                if(currMaxDepth > maxDepth) maxDepth = currMaxDepth;
            }
        }
        mapping_NodesToMaxDepths_.at(node->getId()) = maxDepth;
        return(maxDepth);
    }
    
    void Family::computeTaxonomicShift_(Node* node, Node* father, Node* grandFather, Node* greatGrandFather, vector<Taxon*>* local_nodeToTaxon, vector<Taxon*>* local_GFmappingWithoutTheNode, vector<unsigned int>* local_perturbationsInduced, set<Node*>* local_nodesInducingPerturbations){
        
        vector<Node*> neighbors = node->getNeighbors();
        
        
        for(vector<Node*>::iterator currN = neighbors.begin(); currN != neighbors.end(); currN++){
            if(*currN != father)
                computeTaxonomicShift_(*currN,node,father,grandFather,local_nodeToTaxon, local_GFmappingWithoutTheNode, local_perturbationsInduced, local_nodesInducingPerturbations);
        }
        
        if(greatGrandFather != 00)
        {
            Taxon* initialGfTaxon = local_nodeToTaxon->at(grandFather->getId());
            Taxon* newTaxon = mapNodeOnTaxon_(local_nodeToTaxon,00,grandFather,greatGrandFather,00,true,node);
            //unsigned int currPerturbation = Taxon::computeRelativeDepthDifference(initialGfTaxon,newTaxon,&taxa);
            unsigned int currPerturbation =  newTaxon->getDepth() - initialGfTaxon->getDepth();
            // DBG
            // if(initialGfTaxon != newTaxon) cout << initialGfTaxon->getName() << " . " << newTaxon->getName() << " : " << currPerturbation << endl;
            local_perturbationsInduced->at(node->getId()) = currPerturbation;
            local_GFmappingWithoutTheNode->at(node->getId()) = newTaxon;
            if(initialGfTaxon != newTaxon) local_nodesInducingPerturbations->insert(node);
        } else{
            local_perturbationsInduced->at(node->getId()) = 0;
        }
        
    }
    
    bool Family::transfersRemaining_(){
        //     bool transfersRemaining = false;
        //     for(vector<unsigned int>::iterator currPerturbation = mapping_NodesToTaxonomicShift.begin(); !transfersRemaining && currPerturbation != mapping_NodesToTaxonomicShift.end(); currPerturbation++){
            // 	cout << *currPerturbation << " " << flush;
            // 	transfersRemaining |= (*currPerturbation > 0);
            //     }
            //     cout << (transfersRemaining? "true" : "false") << endl;
            return(!computed_nodesInducingPerturbation_.empty());
    }
    
    
    vector<unsigned int> Family::getDaubinCriteriaBestRoots_(vector<unsigned int> nodes){
        double bestGain = 0;
        vector<unsigned int> bestOutgroups;
        cout << "Trying on " << nodes.size() << " roots"<<endl;
        
        
        for(vector<unsigned int>::iterator currRoot = nodes.begin(); currRoot != nodes.end(); currRoot++){
            tree_->newOutGroup(*currRoot);
            highestID_ = tree_->getNumberOfNodes();
            doMapping_NodesToTaxa();
            doMapping_NodesToTaxonomicShift();
            double gainSum = 0;
            for(set<Node*>::iterator currPerturbator = computed_nodesInducingPerturbation_.begin(); currPerturbator != computed_nodesInducingPerturbation_.end(); currPerturbator++){
                vector<Taxon*> local_nodesToTaxa_perturbatorSubtree(highestID_,00);
                vector<Taxon*> local_nodesToTaxa_treeWithoutPetrurbator(highestID_,00);
                mapNodeOnTaxon_(&mapping_NodesToTaxa_,&local_nodesToTaxa_perturbatorSubtree,*currPerturbator,(*currPerturbator)->getFather(),00,true);
                mapNodeOnTaxon_(&mapping_NodesToTaxa_,&local_nodesToTaxa_treeWithoutPetrurbator,tree_->getRootNode(),00,00,true,*currPerturbator);
                unsigned int taxaSumWithoutPerturbator = getTaxonomicSum_(local_nodesToTaxa_treeWithoutPetrurbator);
                unsigned int taxaSumWithPerturbator = getTaxonomicSum_(mapping_NodesToTaxa_);
                unsigned int taxaSumOfPerturbator = getTaxonomicSum_(local_nodesToTaxa_perturbatorSubtree);
                gainSum += (taxaSumWithPerturbator - taxaSumOfPerturbator) - taxaSumWithoutPerturbator;
                
            }
            
            gainSum /= double(computed_nodesInducingPerturbation_.size());
            
            if (gainSum > bestGain){
                bestGain = gainSum;
                bestOutgroups.clear();
                bestOutgroups.push_back(*currRoot);
            }else if (gainSum == bestGain)
                bestOutgroups.push_back(*currRoot);
            cout<<gainSum << endl;
        }
        return(bestOutgroups);
    }
    
    
    void Family::compute_detectTransfers(){
        unsigned int gainSumFake; // will no be used
        compute_detectTransfers_(true,gainSumFake);
    }
    
    unsigned int Family::compute_detectTransfers_(bool writeResults, unsigned int &gainSum){
        computed_detectedTransfers_.clear();
        unsigned int foundTransfers = 0;
        // doing the first mapping:
        doMapping_NodesToTaxonomicShift();
        
        // then, we keep transfers that are not included in another
        //FIXME it would be better to confirm bootstrap are sufficient
        //for each potential transfer, (1) does this node include no transfer, (2) if yes, transfer confirmed
        // forbiddenNodes contains all the nodes that we know they can’t be transfers, because they include potential transfers
        while(transfersRemaining_()){
            set<Node*> forbiddenNodes;
            for(set<Node *>::iterator perturbator = computed_nodesInducingPerturbation_.begin(); perturbator != computed_nodesInducingPerturbation_.end(); perturbator++){
                // (1.0) it this node in forbidden nodes (it would mean that a descendant node has been found under this node)
                if(forbiddenNodes.find(*perturbator) != forbiddenNodes.end()) continue;
                // (1.1) it is here possible a descendant of this node is a potential transfer, but has not yet been explored
                bool transferFoundInDescendants = false;
                vector<Node*> descendants;
                tpms::TreeTools::getNodesOfTheSubtree(descendants,*perturbator);
                for(vector<Node*>::iterator currDescendant = descendants.begin()+1; !transferFoundInDescendants && currDescendant != descendants.end(); currDescendant++){
                    transferFoundInDescendants |= computed_nodesInducingPerturbation_.find(*currDescendant) != computed_nodesInducingPerturbation_.end();
                }
                if(transferFoundInDescendants) continue;
                
                
                // if we’re here, it means that “node” includes no transfer, we can consider it as a potential transfer, and so remove the subtree from the tree
                // first, we can blacklist all the ancestors (see 1.0)
                Node* nodeToBlacklist = *perturbator;
                while(nodeToBlacklist->hasFather()){
                    nodeToBlacklist = nodeToBlacklist->getFather();
                    forbiddenNodes.insert(nodeToBlacklist);
                }
                
                // then we need to insert the transfer into transfers
                // only if the bootstrap is enough
                
                bool transferAccepted = false;
                Taxon* donnor = mapping_grandFatherWithoutThisNode_.at((*perturbator)->getId());
                vector<string> donnorLeavesNames;
                vector<Node*> donnorLeaves;
                unsigned int perturbationIndex = mapping_NodesToTaxonomicShift_.at((*perturbator)->getId());
                
                double supportOfnodeGroupingIncongruentTaxa = 0;
                
                
                Node* nodeGroupingIncongruentTaxa = (*perturbator)->getFather();
                Node* incongruencyRoot = nodeGroupingIncongruentTaxa->getFather();
                Taxon* incongruencyRootBeforeRemoval =  mapping_NodesToTaxa_.at(incongruencyRoot->getId());

                if(incongruencyRoot == 00){
                    continue;}
                
                // here we now the peturbator node brings a perturbation, we have to test the bootstrap
                if(nodeGroupingIncongruentTaxa->hasBootstrapValue() && nodeGroupingIncongruentTaxa->getBootstrapValue() >= 90 && nodeGroupingIncongruentTaxa->getBootstrapValue() <=100)
                    transferAccepted = true;
                
                // if it’s not the case, we can go deeper in the tree to see if there is a sufficient bootstrap and the incongruency is still there
                
                
                while(!transferAccepted && incongruencyRoot->hasFather()) {
                    // shifting up the studyed nodes
                    nodeGroupingIncongruentTaxa = incongruencyRoot;
                    incongruencyRoot = incongruencyRoot->getFather();
                    
                    // 1st step: is there still an incongruency ? sees the changes of incongruencyRoot affectation
                    Taxon* incongruencyRootOriginalAffectation = mapNodeOnTaxon_(&mapping_NodesToTaxa_,00,incongruencyRoot,(incongruencyRoot->hasFather()?  incongruencyRoot->getFather():00),00,true);
                    incongruencyRootBeforeRemoval = incongruencyRootOriginalAffectation;
                    Taxon* incongruencyRootNewAffectation = mapNodeOnTaxon_(&mapping_NodesToTaxa_,00,incongruencyRoot,(incongruencyRoot->hasFather()?  incongruencyRoot->getFather():00),00,true,*perturbator);
                    donnor = incongruencyRootNewAffectation;
                    // perturbationIndex = tpms::Taxon::computeRelativeDepthDifference(incongruencyRootOriginalAffectation,incongruencyRootNewAffectation,&taxa);
                    perturbationIndex =  incongruencyRootNewAffectation->getDepth() - incongruencyRootOriginalAffectation->getDepth();
                    if(perturbationIndex == 0)
                        break; // because there is no incongruency anymore
                    
                    // 2nd step: checking bootstrap
                    if(nodeGroupingIncongruentTaxa->hasBootstrapValue() && nodeGroupingIncongruentTaxa->getBootstrapValue() >= 90 && nodeGroupingIncongruentTaxa->getBootstrapValue() <=100){
                        transferAccepted = true;
                    }
                    
                }
                
                supportOfnodeGroupingIncongruentTaxa = (nodeGroupingIncongruentTaxa->hasBootstrapValue()? nodeGroupingIncongruentTaxa->getBootstrapValue() : 0);

                
                // TRANSFER ACCEPTED -> TRANSFER RECORDED
                
                if(transferAccepted){
                    
                    // quantifying the gain
                    vector<Taxon*> local_nodesToTaxa_perturbatorSubtree(highestID_,00);
                    vector<Taxon*> local_nodesToTaxa_treeWithoutPetrurbator(highestID_,00);
                    mapNodeOnTaxon_(&mapping_NodesToTaxa_,&local_nodesToTaxa_perturbatorSubtree,*perturbator,(*perturbator)->getFather(),00,true);
                    mapNodeOnTaxon_(&mapping_NodesToTaxa_,&local_nodesToTaxa_treeWithoutPetrurbator,tree_->getRootNode(),00,00,true,*perturbator);
                    unsigned int taxaSumWithoutPerturbator = getTaxonomicSum_(local_nodesToTaxa_treeWithoutPetrurbator);
                    unsigned int taxaSumWithPerturbator = getTaxonomicSum_(mapping_NodesToTaxa_);
                    unsigned int taxaSumOfPerturbator = getTaxonomicSum_(local_nodesToTaxa_perturbatorSubtree);
                    gainSum += (taxaSumWithPerturbator - taxaSumOfPerturbator) - taxaSumWithoutPerturbator;
                    
                    foundTransfers++;
                    transfer currTransfer;
                    currTransfer.donnor = donnor;
                    currTransfer.perturbationIndex = perturbationIndex;
                    currTransfer.bootstrap = supportOfnodeGroupingIncongruentTaxa;
                    currTransfer.formerAncestorTaxon = incongruencyRootBeforeRemoval;
                    atomizeTaxon_(currTransfer.acceptors,currTransfer.acceptorsLeaves,currTransfer.donnor,*perturbator);
                    
                    // the donnor leaves are the leaves under the Grand-Father, ignoring the node *node which is the root of the acceptor group
                    tpms::TreeTools::getNodesOfTheSubtree(donnorLeaves,nodeGroupingIncongruentTaxa,true,*perturbator);
                    // getting names of the leaves
                    for(vector<Node*>::iterator currDonnorLeave = donnorLeaves.begin(); currDonnorLeave != donnorLeaves.end(); currDonnorLeave++)
                        if((*currDonnorLeave)->isLeaf())
                            donnorLeavesNames.push_back((*currDonnorLeave)->getName());
                        currTransfer.donnorLeaves = donnorLeavesNames;
                    
                    computed_detectedTransfers_.push_back(currTransfer);
                }
                
                // ANYWAY, THIS IS AN INCONGRUENCY, WE DELETE THE SUBTREE
                
                // then, we prune the subtree
                Node* P = (*perturbator)->getFather();
                P->removeSon(*perturbator);
                std::vector<Node*> nodesToDelete;
                tpms::TreeTools::getNodesOfTheSubtree(nodesToDelete,*perturbator);
                for(vector<Node*>::iterator currNode = nodesToDelete.begin(); currNode != nodesToDelete.end(); currNode++){
                    set<Taxon*>::iterator taxonToDelete = taxa_.find(mapping_NodesToTaxa_.at((*currNode)->getId()));
                    if(taxonToDelete != taxa_.end())
                        taxa_.erase(taxonToDelete);
                    delete *currNode;
                }
                if(P->getNumberOfSons() == 1){
                    Node* R = *(P->getSons().begin());
                    Node* A = P->getFather();
                    double newBS;
                    newBS = std::max(P->hasBootstrapValue()? P->getBootstrapValue(): 100,R->hasBootstrapValue()? R->getBootstrapValue():100);
                    A->removeSon(P);
                    A->addSon(R);
                    R->setDistanceToFather(P->getDistanceToFather()+R->getDistanceToFather());
                    R->setBranchProperty(bpp::TreeTools::BOOTSTRAP,bpp::Number<double>(newBS));
                    delete(P);
                }
                break;
                
            }
            doMapping_NodesToTaxonomicShift();
            // 	cout << computed_detectedTransfers.size() << ',' << flush;
        }
        
        // now, the gene tree should be congruent with the species tree.
        // outputing the result
        
        //     if(computed_detectedTransfers.size() > 1){
            //     cout << name << " : ";
            //     cout << computed_detectedTransfers.size() << " transfers detected. Original number of leaves : " << mne2tax.size() << ". Final number of leaves : " << tree->getNumberOfLeaves() << endl;}
            
            // NO TRANSFER REMAINING: WRITING ALL THE TRANSFERS TO THE RESULT STRING
            if(writeResults){
                        for(vector<transfer>::iterator currTransfer = computed_detectedTransfers_.begin(); currTransfer != computed_detectedTransfers_.end(); currTransfer++){
                            results_ << name_ << ',' << currTransfer->perturbationIndex << ',' << currTransfer->donnor->getName() << ",(" ;
                            string resultStr;
                            for(vector<Taxon*>::iterator currReceiver = currTransfer->acceptors.begin(); currReceiver != currTransfer->acceptors.end(); currReceiver++)
                                resultStr += (*currReceiver)->getName() + ",";
                            if(resultStr.empty())
                                resultStr = "<undetermined taxon> ";
                            resultStr.at(resultStr.size()-1) = ')';
                    results_ << resultStr << endl;
                    
                    // now giving details
                    // donnor leaves
                    results_ << "; brother-leaves:";
                    char separator = ' ';
                    for(vector<string>::iterator currDonnorLeave = currTransfer->donnorLeaves.begin(); currDonnorLeave != currTransfer->donnorLeaves.end(); currDonnorLeave++){
                        results_ << separator << *currDonnorLeave;
                        separator = ',';
                    }
                    results_ << endl;
                    
                    results_ << "; donnor-taxonomy: " + currTransfer->donnor->getTaxonomy() << endl;
                    
                    // acceptor(s) leaves
                    // for each acceptor
                    for(unsigned int currAcceptor = 0; currAcceptor != currTransfer->acceptorsLeaves.size(); currAcceptor++){
                        results_ << "; acceptor-leaves-" << currAcceptor << ":";
                        separator = ' ';
                    for(vector<string>::iterator currAcceptorLeave = currTransfer->acceptorsLeaves.at(currAcceptor).begin(); currAcceptorLeave != currTransfer->acceptorsLeaves.at(currAcceptor).end(); currAcceptorLeave++){
                        results_ << separator << *currAcceptorLeave;
                        separator = ',';
                    }
                    results_ << endl;
                    
                    results_ << "; acceptor-taxonomy-" << currAcceptor << ": " << currTransfer->acceptors.at(currAcceptor)->getTaxonomy() << endl;
                    results_ << "; bootstrap: " << currTransfer->bootstrap << endl;
                    results_ << "; former-ancestor-taxon: " << currTransfer->formerAncestorTaxon->getName() << endl;
                    }
                    
                    results_ << ";" << endl;
                        }
                        } 
            
        return(foundTransfers);
    }
    
    void Family::atomizeTaxon_(std::vector< Taxon* > &acceptors, vector<vector<string> > &acceptorsLeaves, Taxon* acceptor, Node* subtree)
    {
        
        Taxon* currTaxon = mapping_NodesToTaxa_.at(subtree->getId());
        // base case case
        if(acceptor != currTaxon && acceptor->getAncestors().find(currTaxon) == acceptor->getAncestors().end() && currTaxon != 00){
            // base case: here currTaxon is not included in ancestor’s ancestors
            vector<Node*> currAcceptorNodes;
            tpms::TreeTools::getNodesOfTheSubtree(currAcceptorNodes, subtree, true);
            vector<string> currAcceptorLeavesNames;
            for(vector<Node*>::iterator currN = currAcceptorNodes.begin(); currN != currAcceptorNodes.end(); currN++){
                if((*currN)->hasName())
                    currAcceptorLeavesNames.push_back((*currN)->getName());
            }
            acceptors.push_back(currTaxon);
            acceptorsLeaves.push_back(currAcceptorLeavesNames);
        }else{
            vector<Node*> sons = subtree->getSons();
            for(vector<Node*>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++)
                atomizeTaxon_(acceptors,acceptorsLeaves, acceptor, *currSon);
        } 
        
    }
    
    
    
    void Family::threadedWork_launchJobs(std::vector<Family *> families, void (Family::*function)(), unsigned int nbThreads, ostream *output){
        unsigned int nbFamilies = families.size();
        boost::mutex progressbarMutex;
        boost::mutex outputMutex;
        boost::thread_group tg;
        unsigned int blockSize = nbFamilies / nbThreads;
        cout << "\nMultithreaded operation. Number of threads: " << nbThreads << ". Lot size : " << blockSize << endl;Waiter progressbar(&cout, nbFamilies, '#');
        
        for(unsigned int currThreadIndex = 0; currThreadIndex < nbThreads; currThreadIndex++){
            vector<Family*>::iterator currPartBegin;
            vector<Family*>::const_iterator currPartEnd;
            currPartBegin = families.begin() + (blockSize*currThreadIndex);
            if(currThreadIndex+1 != nbThreads) currPartEnd = families.begin() + (blockSize*(currThreadIndex+1)); else currPartEnd = families.end();
                    boost::thread *currThread = new boost::thread(Family::threadedWork_oneThread,function,&progressbar,&progressbarMutex,output,&outputMutex,currPartBegin,currPartEnd);
            tg.add_thread(currThread);
        }
        tg.join_all();
        progressbar.drawFinal();
    }
    
    
    void Family::threadedWork_patternMatching(std::vector<Family *> & families, Pattern * pattern, vector<pair<Family*, CandidateNode*> > * matchingFamilies, unsigned int nbThreads){
        unsigned int nbFamilies = families.size();
        boost::mutex progressbarMutex;
        boost::mutex outputMutex;
        boost::thread_group tg;
        unsigned int blockSize = nbFamilies / nbThreads;
        cout << "\nMultithreaded operation. Number of threads: " << nbThreads << ". Lot size : " << blockSize << endl;Waiter progressbar(&cout, nbFamilies, '#');
        
        for(unsigned int currThreadIndex = 0; currThreadIndex < nbThreads; currThreadIndex++){
            vector<Family*>::iterator currPartBegin;
            vector<Family*>::const_iterator currPartEnd;
            currPartBegin = families.begin() + (blockSize*currThreadIndex);
            if(currThreadIndex+1 != nbThreads) currPartEnd = families.begin() + (blockSize*(currThreadIndex+1)); else currPartEnd = families.end();
                    boost::thread *currThread = new boost::thread(Family::threadedWork_onePMThread,pattern,matchingFamilies,&progressbar,&progressbarMutex,&outputMutex,currPartBegin,currPartEnd);
            tg.add_thread(currThread);
        }
        tg.join_all();
        progressbar.drawFinal();
    }
    
    void Family::threadedWork_oneThread(void(Family::*function)(),Waiter *progressbar, boost::mutex *progressbarMutex, ostream *output, boost::mutex *outputMutex , std::vector<tpms::Family*>::iterator &currPartBegin, std::vector<tpms::Family*>::const_iterator &currPartEnd){
        unsigned int waiterUpdate = 0;
        for(vector<Family*>::iterator currFamily = currPartBegin; currFamily != currPartEnd; currFamily++){
            if(!(*currFamily)->getContainsUndefinedSequences()){
                (*currFamily->*function)();
                if(output != 00){
                    outputMutex->lock();
                    *output << (*currFamily)->threadedWork_getResults().str();
                    outputMutex->unlock();
                }
            }
            waiterUpdate++;
            if(waiterUpdate == 50){
                waiterUpdate = 0;
                progressbarMutex->lock();
                progressbar->doStep(50);
                progressbarMutex->unlock();
            }
        }
    }
    
    void Family::threadedWork_onePMThread(Pattern * pattern, vector<pair<Family *, CandidateNode*> > *matchingFamilies, Waiter *progressbar, boost::mutex *progressbarMutex, boost::mutex *outputMutex , std::vector<tpms::Family*>::iterator &currPartBegin, std::vector<tpms::Family*>::const_iterator &currPartEnd){
        unsigned int waiterUpdate = 0;
        
        for(vector<Family*>::iterator currFamily = currPartBegin; currFamily != currPartEnd; currFamily++){
            if(!(*currFamily)->getContainsUndefinedSequences()){
                vector<pair<Family *, CandidateNode*> > resultVector;
                CandidateNode * candRoot = new CandidateNode();
                if(pattern->patternMatchInit(**currFamily,candRoot)){
                    resultVector.push_back(pair<Family *, CandidateNode*>(*currFamily,candRoot));
                }
                else
                    delete(candRoot);
                
                outputMutex->lock();
                matchingFamilies->insert(matchingFamilies->begin(),resultVector.begin(),resultVector.end());
                outputMutex->unlock();
            }
            waiterUpdate++;
            if(waiterUpdate == 50){
                waiterUpdate = 0;
                progressbarMutex->lock();
                progressbar->doStep(50);
                progressbarMutex->unlock();
            }
        }
    }
    
    bool Family::getContainsUndefinedSequences(){
        return(containsUndefinedSequences_);
    }
    
    
    
    ostringstream& Family::threadedWork_getResults()
    {
        return(results_);
    }
    
    
}
