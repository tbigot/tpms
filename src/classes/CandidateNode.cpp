#include "CandidateNode.hpp"

using namespace std;
using namespace tpms;
using namespace bpp;


namespace tpms{
CandidateNode::~CandidateNode(){
    // deleting sons
    for(map<bpp::Node*,vector<CandidateNode*> >::iterator currSons = sons.begin(); currSons != sons.end(); currSons++){
	for(vector<CandidateNode*>::iterator currSon = currSons->second.begin(); currSon != currSons->second.end(); currSon++){
	    delete(*currSon);
	}
    }
}

CandidateNode::CandidateNode(CandidateNode* father, Node* treeNode, Node* patternNode):
father(father), treeNode(treeNode), patternNode(patternNode), isRoot(false), addsons(0)
{
}

Node * CandidateNode::getTreeNode(){ return(treeNode);}

void CandidateNode::confirm(){
    if(sons.size() == 2 || sons.size() == 0)
	father->addSon(this,patternNode);
    else {
	cout << "Erreur Nb Fils : " << sons.size() << " et addsons=" << addsons << ". " ;
	if(patternNode->hasFather()) cout << "pas racine";
	cout << endl;
    }
}

CandidateNode::CandidateNode(): isRoot(true), addsons(0)
{
}

bool CandidateNode::isLeaf()
{
    return(sons.size() == 0);
}




void CandidateNode::addSon(CandidateNode* son, Node* sonPatternNode)
{
    addsons++;
    sons[sonPatternNode].push_back(son);
}


unsigned int CandidateNode::genTrees(std::vector< TreeTemplate< Node >* >& trees)
{
    vector<Node *> roots;
    if(isRoot){
	//NOTE: the candidate root has no phylogenetics meaning. It's only designed to put found trees together.
		
	vector<CandidateNode*> candSons = sons.begin()->second;
	for(vector<CandidateNode*>::iterator currSon = candSons.begin(); currSon != candSons.end(); currSon++){
	    vector<Node *> currRoots = (*currSon)->recGenTrees();
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
    if(isRoot){
	//NOTE: the candidate root has no phylogenetics meaning. It's only designed to put found trees together.
		
	vector<CandidateNode*> candSons = sons.begin()->second;
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
    if(!isRoot){
    cout << this << ' ' << patternNode->getId();
    if(patternNode->hasName()) cout <<" " << patternNode->getName() << " : ";
    } else cout <<'\n'<< this << ' ';
    for(map<Node *,vector<CandidateNode *> >::iterator it = sons.begin(); it != sons.end(); it++){
	cout << it->first << ", " ;
	
    }
    if(isLeaf()) cout << " <--------- " << endl; else { cout << endl;
    for(map<Node *,vector<CandidateNode *> >::iterator it = sons.begin(); it != sons.end(); it++){
	for(unsigned int i = 0; i < lvl; i++) cout << '+';
	cout << ">" << it->first << '-' << endl;
	
	for(vector<CandidateNode *>::iterator it2 = it->second.begin(); it2!= it->second.end(); it2++){
	    (*it2)->print(lvl+1);
	}
    }
    }
}

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
