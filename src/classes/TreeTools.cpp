//
// File: TreeTools.cpp
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

#include <string>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTools.h>

#include "TreeTools.hpp"

using namespace std;
using namespace bpp;

namespace tpms{
bpp::TreeTemplate<bpp::Node> * TreeTools::newickToTree(string & in, bool distanceDelimitation){
	unsigned int jogger = 0;
	unsigned int nodeId = 0;
	Node * racine = newickToNode(in, &jogger, 00, nodeId,distanceDelimitation);
	TreeTemplate<Node> * retour = new TreeTemplate<Node>(racine);
	// cout << "Nous venons de créer un arbre a partir de Newick" << endl;
	//treePrint(racine,0);
	return(retour);
}

Node * TreeTools::newickToNode(string & in, unsigned int * jogger, Node * father, unsigned int& nodeId, bool distanceDelimitation){
	stringstream *ssNom = new stringstream();
	stringstream *ssDist = new stringstream();
	double distance;
	bool sharp = false; // duplication marker for reconciled trees
	Node * thisNode = 00;
	while(in.at(*jogger) != ')' && in.at(*jogger) != ';') { // énumération de tous les éléments d'une groupe parenthésé
		thisNode = new Node(nodeId);
		nodeId++;
		if(father != 00) father->addSon(thisNode);
		
		// case of DUPLICATION marker
		if(in.at(*jogger) == '#') {
		    sharp= true;
		    (*jogger)++;
		}
		
		if(in.at(*jogger) == '(') { //cas d'un groupe, on va faire un appel récursif pour parser le contenu de la parenthèse
				(*jogger)++;
				newickToNode(in,jogger,thisNode,nodeId,distanceDelimitation);
				// la fonction nous rend la main sur une fin de groupe : on est forcément sur une parenthèse fermante
				(*jogger)++; // donc on avance d'un cran pour se trouver sur le premier char du nom
		}
		
		// on se retrouve ici : soit un groupe qui est fini, soit on était sur une feuille. Dans tous les cas, c'est le nom du nœud créé dans cette itération qu'on va chercher
		
		ssDist->clear();
		ssDist->str("");
		ssNom->clear();
		ssNom->str("");
				
		// le nom n'est pas terminé tant qu'on ne tombe pas sur deux points
		// ou alors si /!DistanceDelimitation/, tout autre signe de fin de champ
		while(in.at(*jogger) != ':' && (distanceDelimitation || in.at(*jogger) != ',' && in.at(*jogger) != ')' && in.at(*jogger) != ';')){
			*ssNom << in.at(*jogger);
			(*jogger)++;
		}
		
		if(in.at(*jogger) == ':'){ // traitement éventuel de la distance
		    // maintenant, la distance, qui n'est pas terminée tant qu'on ne tombe pas sur une virgule ou une parenthèse fermante' ou un point-virgule pour le dernier groupe
		    // on est sur deux points, donc on avance d'un caractere
		    (*jogger)++;
		    while(in.at(*jogger) != ',' && in.at(*jogger) != ')' && in.at(*jogger) != ';'){
			    *ssDist << in.at(*jogger);
			    (*jogger)++;
		    }
		}
		if(!ssNom->str().empty()) {
		    unsigned int bootstrap=0;
		    thisNode->setName(ssNom->str());
		    *ssNom >> bootstrap;
		    if(bootstrap != 0) {
			thisNode->setBranchProperty(bpp::TreeTools::BOOTSTRAP, Number<double>(bootstrap));
		    }
		}
		else if(sharp) thisNode->setName("#");
		*ssDist >> distance;
		thisNode->setDistanceToFather(distance);
		
		if(in.at(*jogger) == ',') (*jogger)++;
		
	}
	delete(ssNom);
        delete(ssDist);
	return(thisNode);
	
}

void TreeTools::treePrint(Node * noeud, unsigned int cpt) {
	// cette fonction récursive affiche l'arbre sous forme hiérarchique espacée. cpt correspond au nombre d'espaces (profondeur en cours)
	for(unsigned int i=0; i<cpt; i++) cout << ".";
	cout << noeud->getId();
	if(noeud->hasName()) cout << '|' << noeud->getName(); else cout << "<no name>";
	cout << " (" << noeud->getNumberOfSons() << ")" << endl;
	for(unsigned int i=0; i< noeud->getNumberOfSons(); i++) treePrint(noeud->getSon(i),cpt+1);
}

bool TreeTools::isBinaryTree(Node * node){
    switch(node->getNumberOfSons()){
	case 0: return(true);
	case 2: return(isBinaryTree(node->getSon(0)) && isBinaryTree(node->getSon(1)));
	default: return(false);
    }
}

bool TreeTools::isAtLeastBinaryTree(Node * node){
    unsigned int numberOfSons = node->getNumberOfSons();
    if(numberOfSons ==0) return(true);
    if(numberOfSons >=2){
	bool result = true;
	for(unsigned int i = 0; result && i < numberOfSons; i++) result = isAtLeastBinaryTree(node->getSon(i));
	return(result);
    }
    return(false);
}

std::vector<unsigned int> TreeTools::toNodesIDs(std::vector<bpp::Node*> &nodes){
    vector<unsigned int> IDs;
    for(auto node:nodes)
        IDs.push_back(node->getId());
    return IDs;
}

std::vector<Node*> TreeTools::toNodes(bpp::TreeTemplate<Node> &tree, std::vector<unsigned int> &IDs){ 
    vector<Node*> nodes;
    for(auto ID:IDs)
        nodes.push_back(tree.getNode(ID));
    return nodes;
}


string TreeTools::nodeToNewick(Node* node)
{
    ostringstream result;
    vector<Node *> sons = node->getSons();
    if(sons.size() !=0){
	result << '(';
	bool first = true;
	for(vector<Node *>::iterator son = sons.begin(); son != sons.end(); son++){
	    if(!first) result << ','; else first = false;
	    result << nodeToNewick(*son);
	}
	result << ')';
    }
    if(node->hasName()) result << node->getName();
    if(node->hasDistanceToFather()) result << ":" << node->getDistanceToFather();
    return(result.str());
}

void TreeTools::multifurcated2binary(TreeTemplate< Node > * currTree, unsigned int parentID, vector< TreeTemplate< Node > *> &trees)
{
    Node * parent = currTree->getNode(parentID);
    vector<Node *> sons = parent->getSons();
    unsigned int nextID = currTree->getNextId();
    if(sons.size() == 2){
	multifurcated2binary(currTree, sons.at(0)->getId(), trees);
	multifurcated2binary(currTree, sons.at(1)->getId(), trees);
    } else if (sons.size() > 2){
	vector<Node *>::iterator currSon;
	
	// first sons : we have to create new trees
	for(currSon = sons.begin(); currSon != sons.end(); currSon++){
	    TreeTemplate<Node>* newTree;
	    if(currSon != sons.end()-1) {
		newTree = new TreeTemplate<Node>(*currTree);
		trees.push_back(newTree);
	    }
	    else newTree = currTree;
	    
	    Node * newParent = newTree->getNode(parentID);
	    Node * localCurrSon = newTree->getNode((*currSon)->getId());
	    Node * newSon = new Node(newTree->getNextId());
	    newSon->setId(nextID);
	    vector<Node *> brothers = newParent->getSons();
	    for(vector<Node *>::iterator brother = brothers.begin(); brother != brothers.end(); brother++)
		if(*brother != localCurrSon){
		    newParent->removeSon(*brother);
		    newSon->addSon(*brother);
		}
	    newParent->addSon(newSon);
	    multifurcated2binary(newTree, newSon->getId(), trees);
		
	}
	
    }
}


string TreeTools::extractNewickLineFromFile(ifstream& in)
{
    string currline;
    bool beforeComments = true;
    bool afterComments = false;
    while(getline(in,currline)){
	if(beforeComments){
	    if(currline.at(0) == '[') beforeComments=false;
	} else if(currline.at(0) == ']') { afterComments=true;}
	if(beforeComments || afterComments){
	    if(currline.find('(') != string::npos && currline.find(')') != string::npos && currline.find(';') != string::npos){
		return(currline);
	    }
	}
    }
    return(string());
}


double TreeTools::getDistanceBetweenTwoNodes(Node* ancestor, Node* descendant)
{
    // root case : ancestor has not been met
    if(!descendant->hasFather() || descendant == ancestor) return(0);
    return(descendant->getDistanceToFather() + getDistanceBetweenTwoNodes(ancestor, descendant->getFather()));
}


void TreeTools::getNodesOfTheSubtree(vector<Node*> &nodeListToFill,Node* node, bool leavesOnly, Node* ignoreSubtree)
{
    if(node == ignoreSubtree)
        return;
    if(!leavesOnly || node->isLeaf())
            nodeListToFill.push_back(node);
    vector<Node*> sons = node->getSons();
    for(vector<Node*>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++){
	getNodesOfTheSubtree(nodeListToFill,*currSon, leavesOnly, ignoreSubtree);
    }
}

void TreeTools::destroySubtree(Node* node)
{
    for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
    {
	Node* son = node->getSon(i);
	destroySubtree(son);
	delete son;
    }
}

unsigned int TreeTools::depthOfANode(Node* node)
{
    if(node->hasFather()) return(1+depthOfANode(node->getFather()));
    else return(0);
}


} //fin namespace tpms
