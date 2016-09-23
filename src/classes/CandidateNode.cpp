//
// File: CandidateNode.cpp
// Created by: Thomas Bigot
//

/*
   Copyright or © or Copr. Thomas Bigot 2012

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "CandidateNode.hpp"
#include "TreeTools.hpp"

using namespace std;
using namespace tpms;
using namespace bpp;


namespace tpms{
CandidateNode::~CandidateNode(){
    // deleting sons
    for(map<bpp::Node*,vector<CandidateNode*> >::iterator currSons = sons_.begin(); currSons != sons_.end(); currSons++){
	for(vector<CandidateNode*>::iterator currSon = currSons->second.begin(); currSon != currSons->second.end(); currSon++){
	    delete(*currSon);
	}
    }
}

CandidateNode::CandidateNode(CandidateNode* father, Node* treeNode, Node* patternNode):
father_(father), treeNode_(treeNode), patternNode_(patternNode), isRoot_(false), addsons_(0)
{
}

Node * CandidateNode::getTreeNode(){ return(treeNode_);}

void CandidateNode::confirm(){
    if(sons_.size() == 2 || sons_.size() == 0)
	father_->addSon(this,patternNode_);
    else {
	cout << "Erreur Nb Fils : " << sons_.size() << " et addsons=" << addsons_ << ". " ;
	if(patternNode_->hasFather()) cout << "pas racine";
	cout << endl;
    }
    
    // setting the distance to the matching node father
    distanceToFather_ = tpms::TreeTools::getDistanceBetweenTwoNodes(father_->getTreeNode(),treeNode_);
}

CandidateNode::CandidateNode(): isRoot_(true), addsons_(0)
{
}

bool CandidateNode::isLeaf_()
{
    return(sons_.size() == 0);
}




void CandidateNode::addSon(CandidateNode* son, Node* sonPatternNode)
{
    addsons_++;
    sons_[sonPatternNode].push_back(son);
}


unsigned int CandidateNode::genTrees(std::vector< TreeTemplate< Node >* >& trees)
{
    vector<Node *> roots;
    if(isRoot_){
	//NOTE: the candidate root has no phylogenetics meaning. It's only designed to put found trees together.
		
	vector<CandidateNode*> candSons = sons_.begin()->second;
	for(vector<CandidateNode*>::iterator currSon = candSons.begin(); currSon != candSons.end(); currSon++){
	    vector<Node *> currRoots = (*currSon)->recGenTrees_();
	    roots.insert(roots.end(),currRoots.begin(),currRoots.end());
	}
	
	
	for(vector<Node *>::iterator cr = roots.begin(); cr != roots.end(); cr++){
	    trees.push_back(new TreeTemplate<Node>(*cr));
	}
    }
    return(roots.size());
}

unsigned int CandidateNode::getWholeMatchingSubtrees(std::vector< TreeTemplate< Node >* >& trees)
{
    vector<Node *> roots;
    if(isRoot_){
	//NOTE: the candidate root has no phylogenetics meaning. It's only designed to put found trees together.
		
	vector<CandidateNode*> candSons = sons_.begin()->second;
	for(vector<CandidateNode*>::iterator currSon = candSons.begin(); currSon != candSons.end(); currSon++){
	    roots.push_back(bpp::TreeTemplateTools::cloneSubtree<Node>(*(*currSon)->getTreeNode()));
	}
	
	
	for(vector<Node *>::iterator cr = roots.begin(); cr != roots.end(); cr++){
	    trees.push_back(new TreeTemplate<Node>(*cr));
	}
    }
    return(roots.size());
}


void CandidateNode::print(unsigned int lvl){
    for(unsigned int i = 0; i <= lvl; i++) cout << ' ';
    if(!isRoot_){
    cout << this << ' ' << patternNode_->getId();
    if(patternNode_->hasName()) cout <<" " << patternNode_->getName() << " : ";
    } else cout <<'\n'<< this << ' ';
    for(map<Node *,vector<CandidateNode *> >::iterator it = sons_.begin(); it != sons_.end(); it++){
	cout << it->first << ", " ;
	
    }
    if(isLeaf_()) cout << " <--------- " << endl; else { cout << endl;
    for(map<Node *,vector<CandidateNode *> >::iterator it = sons_.begin(); it != sons_.end(); it++){
	for(unsigned int i = 0; i < lvl; i++) cout << '+';
	cout << ">" << it->first << '-' << endl;
	
	for(vector<CandidateNode *>::iterator it2 = it->second.begin(); it2!= it->second.end(); it2++){
	    (*it2)->print(lvl+1);
	}
    }
    }
}

vector<Node *> CandidateNode::recGenTrees_(){
    vector<Node *> result;
    
    Node * currNode = new Node;
    currNode->setId(treeNode_->getId());
    currNode->setDistanceToFather(distanceToFather_);
    if(patternNode_->hasName()) currNode->setName(patternNode_->getName());
    
    result.push_back(currNode);
    
    if(!isLeaf_()){
	
	// generating subtrees
	vector<Node *> leftSubtrees;
	vector<Node *> rightSubtrees;
	
	vector<CandidateNode *> leftSons = sons_.begin()->second;
	vector<CandidateNode *> rightSons = sons_.rbegin()->second;
	
	for(vector<CandidateNode *>::iterator cc = leftSons.begin(); cc != leftSons.end(); cc++){
	    vector<Node*> currResult = (*cc)->recGenTrees_();
	    leftSubtrees.insert(leftSubtrees.end(),currResult.begin(),currResult.end());
	}
	
	for(vector<CandidateNode *>::iterator cc = rightSons.begin(); cc != rightSons.end(); cc++){
	    vector<Node*> currResult = (*cc)->recGenTrees_();
	    rightSubtrees.insert(rightSubtrees.end(),currResult.begin(),currResult.end());
	}
	
	// STEP 1: making all possible pairs of son
	vector<pair<Node*,Node*> > sonsCombinations;
	for(vector<Node *>::iterator leftSubtree = leftSubtrees.begin(); leftSubtree != leftSubtrees.end(); leftSubtree++) {
	    for(vector<Node *>::iterator rightSubtree = rightSubtrees.begin(); rightSubtree != rightSubtrees.end(); rightSubtree++) {
		sonsCombinations.push_back(pair<Node*,Node*>(*leftSubtree,*rightSubtree));
	    }
	}
		
	vector<pair<Node*,Node*> >::iterator currCombination;
	
	set<Node *> toClone;
	
	
	Node *leftSon,*rightSon;
	
	for(currCombination = sonsCombinations.begin(); currCombination != sonsCombinations.end(); currCombination++){
	    Node * newCurrNode;
	    if(currCombination != sonsCombinations.end()-1){
		newCurrNode = new Node(*currNode);
		result.push_back(newCurrNode);
	    } else newCurrNode = currNode;
	    leftSon = currCombination->first;
	    if(toClone.find(leftSon) != toClone.end()) leftSon = bpp::TreeTemplateTools::cloneSubtree<Node>(*leftSon);
	    else toClone.insert(leftSon);
	    rightSon = currCombination->second;
	    if(toClone.find(rightSon) != toClone.end()) rightSon = bpp::TreeTemplateTools::cloneSubtree<Node>(*rightSon);
	    else toClone.insert(rightSon);
	    newCurrNode->addSon(leftSon);
	    newCurrNode->addSon(rightSon);
	}
	
    }
    
    
    return(result);
}

/*
vector<Node *> CandidateNode::recGenTrees(){
    vector<Node *> result;
    
    Node * currNode = new Node;
    currNode->setId(treeNode->getId());
    if(patternNode->hasName()) currNode->setName(patternNode->getName());
    
    result.push_back(currNode);
    
    if(!isLeaf()){
	
	// generating subtrees
	vector<Node *> leftSubtrees;
	vector<Node *> rightSubtrees;
	
	vector<CandidateNode *> leftSons = sons.begin()->second;
	vector<CandidateNode *> rightSons = sons.rbegin()->second;
	
	for(vector<CandidateNode *>::iterator cc = leftSons.begin(); cc != leftSons.end(); cc++){
	    vector<Node*> currResult = (*cc)->recGenTrees();
	    leftSubtrees.insert(leftSubtrees.end(),currResult.begin(),currResult.end());
	}
	
	for(vector<CandidateNode *>::iterator cc = rightSons.begin(); cc != rightSons.end(); cc++){
	    vector<Node*> currResult = (*cc)->recGenTrees();
	    rightSubtrees.insert(rightSubtrees.end(),currResult.begin(),currResult.end());
	}
	
	// STEP 1: making all possible pairs of son
	vector<pair<Node*,Node*> > sonsCombinations;
	for(vector<Node *>::iterator leftSubtree = leftSubtrees.begin(); leftSubtree != leftSubtrees.end(); leftSubtree++) {
	    for(vector<Node *>::iterator rightSubtree = rightSubtrees.begin(); rightSubtree != rightSubtrees.end(); rightSubtree++) {
		sonsCombinations.push_back(pair<Node*,Node*>(*leftSubtree,*rightSubtree));
	    }
	}
		
	vector<pair<Node*,Node*> >::iterator currCombination;
	
	set<Node *> toClone;
	
	
	Node *leftSon,*rightSon;
	
	for(currCombination = sonsCombinations.begin(); currCombination != sonsCombinations.end(); currCombination++){
	    Node * newCurrNode;
	    if(currCombination != sonsCombinations.end()-1){
		newCurrNode = new Node(*currNode);
		result.push_back(newCurrNode);
	    } else newCurrNode = currNode;
	    leftSon = currCombination->first;
	    if(toClone.find(leftSon) != toClone.end()) leftSon = bpp::TreeTemplateTools::cloneSubtree<Node>(*leftSon);
	    else toClone.insert(leftSon);
	    rightSon = currCombination->second;
	    if(toClone.find(rightSon) != toClone.end()) rightSon = bpp::TreeTemplateTools::cloneSubtree<Node>(*rightSon);
	    else toClone.insert(rightSon);
	    newCurrNode->addSon(leftSon);
	    newCurrNode->addSon(rightSon);
	}
	
    }
    
    
    return(result);
}*/
}
