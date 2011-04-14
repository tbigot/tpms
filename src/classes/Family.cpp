#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <numeric>

#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>

#include "Family.hpp"
#include "DataBase.hpp"
#include "Waiter.hpp"
#include "TreeTools.hpp"


using namespace std;
using namespace bpp;

Family::Family(stringstream* sIntro, string sNewick, DataBase* dbp): db(dbp) {
    // extraction du nom de la famille
    string currLigne;
    getline(*sIntro,currLigne);
    istringstream sFamilleName(currLigne);
    sFamilleName >> _name;
    // extraction des synonymes
    getline(*sIntro, currLigne);
    if(currLigne.at(0) != '[') std::cerr << "Bracket (newick comments beginning) expected at family " << _name << endl;
    
    // cas du fichier avec deux arbres en plus (concerne les arbres réconciliés
    getline(*sIntro, currLigne);
    if(currLigne.at(0) == '#') { getline(*sIntro, currLigne);getline(*sIntro, currLigne);getline(*sIntro, currLigne);} // on saute ces lignes
	
	bool mneFini;
	bool specFini;
	ostringstream mnemonique;
	ostringstream espece;
	char lecar;
	
	while(currLigne.at(0) != ']') {
	    mnemonique.str("");
	    espece.str("");
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
	    
	    mne2spec.insert(pair<string,string>(mnemonique.str(),espece.str()));
	    species.insert(espece.str());
	    
	    //DEBUG: on vérifie que espece est bien dans la liste des espèces de la base !
	    if(db->getSpecies()->find(espece.str()) == db->getSpecies()->end())
		cout << "o Unable to find this specie in the species tree :" << espece.str() << endl;
	    
	    
	    getline(*sIntro, currLigne);
	}
	
	// on fabrique l'arbre
	//istringstream ssNewick(sNewick);
	//Newick newickReader(false);
	//tree = newickReader.read(ssNewick);
	
	tree = tpms::TreeTools::newickToTree(sNewick,false);
	
	/*if(tree->getNumberOfNodes() != tree2->getNumberOfNodes()) cout << "Incoherence famille " << _name << " bpp=" << tree->getNumberOfNodes() << " tpms=" <<tree2->getNumberOfNodes() << endl  ;
	*/
}

void Family::genSpTree(bool save, string path) {
    spTree = tree->clone();
    vector<Node *> noeuds = spTree->getLeaves();
    for(vector<Node *>::iterator it = noeuds.begin(); it < noeuds.end(); it++) {
	if((*it)->hasName())
	    (*it)->setName(mne2spec[(*it)->getName()]);
    }
    //DEBUG: on va sauvegarder temporairement ces fichiers
    //ofstream sauveRefTree(string("cache/"+_name+".sptree").c_str());
    //writeTreeToStream(spTree->getRootNode(),sauveRefTree,0);
    
    if (save) writeSpTreeToFile(path);
}

void Family::genUnicityScores() {
    unicityScores.resize(tree->getNumberOfNodes());
    computeUnicity(unicityScores,tree->getRootNode(),00); // 00 because no origin, we are at the real root
}

void Family::genBestUnicityScores() {
    // in this function, we'll try all nodes as potential roots, and keep the topology that gets the lower scores sum
    // then, the tree is rerooted
    
    vector<unsigned int> currScores;
    vector<unsigned int> bestScores;
    Node * bestRoot = 00;
    unsigned int bestScoresSum;
    unsigned int currScoresSum;
    currScores.resize(tree->getNumberOfNodes());
    vector<Node *> nodes = tree->getNodes();
    for(vector<Node *>::iterator currNode = nodes.begin(); currNode != nodes.end(); currNode++) {
	if((*currNode)->isLeaf()) continue; // leaves cannot be roots, skipping
	computeUnicity(currScores,*currNode,00);
	currScoresSum = accumulate(currScores.begin(),currScores.end(),0);
	if(bestRoot == 00 || currScoresSum <  bestScoresSum){
	    // we've found a first or better root candidate
	    bestRoot = *currNode;
	    bestScores = currScores;
	    bestScoresSum = currScoresSum;
	}
    }
    // here, we've the best root in bestRoot: we have to reRoot
    // indicating a new outgroup (not to have a trifurcated root)
    // FIXME:tree->rootAt(bestRoot);
    
    vector<Node *> sons = bestRoot->getNeighbors();
    Node* bestOutgroup = 00;
    unsigned int minScore = 0;
    for(vector<Node *>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++){
	if(bestOutgroup == 00 || bestScores.at((*currSon)->getId())<minScore ){
	    minScore = bestScores.at((*currSon)->getId());
	    bestOutgroup = *currSon;
	}
    }
    
    // "rerooting" the tree according to this new outgroup :
    tree->newOutGroup(bestOutgroup);
    
    // after reroot, nodes id might have changed, we must recalculate scores
    genUnicityScores();
    
}

map<string, unsigned int> Family::computeUnicity(vector<unsigned int> &scores, Node * node, Node * origin){
    unsigned int id = node->getId();
    map<string, unsigned int> thisNodeCount;
    
    // step 1 : this node count
        
    vector<Node *> neighbors = node->getNeighbors();
    for(vector<Node *>::iterator currSon = neighbors.begin(); currSon < neighbors.end(); currSon ++) {
	if(*currSon == origin) continue;
	map<string,unsigned int> currSonCount = computeUnicity(scores,*currSon,node);
	// adding this new map to the current
	for(map<string,unsigned int>::iterator currCount = currSonCount.begin(); currCount != currSonCount.end(); currCount++){	    
	    thisNodeCount[currCount->first] += currCount->second;
	}
    }
    
    if(neighbors.size() == 1){ // leave case
	    thisNodeCount.insert(pair<string,unsigned int>(spTree->getNodeName(id),1));
	    
    }
    
    // step 2 : this node score computation
    unsigned int score = 1;
    
    for(map<string,unsigned int>::iterator currCount = thisNodeCount.begin(); currCount != thisNodeCount.end(); currCount++){
	 score *= currCount->second;
    }    
    
    scores[id] = score;
        
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

void Family::genRefTree(bool save, string path) {
    // l'arbre est initialement celui de référence de la base de données
    refTree = db->getSpeciesTree()->clone();
    
    //ensuite, pour chaque feuille, on regarde si elle est dans la liste des espèces de la famille
    vector<Node*> leaves = refTree->getLeaves();
    
    for(vector<Node*>::iterator currLeaf = leaves.begin(); currLeaf != leaves.end(); currLeaf++){
	if(species.find((*currLeaf)->getName()) == species.end())
	    //si elle ne l'est pas, on supprime le nœud feuille et on remonte jusqu'à la prochaine bifurcation
	    deleteFromLeavesToBif(*currLeaf);
	
	    
    }
    
   // on va ensuite parcourir récursivement l'arbre et éliminer les fils uniques. 
   removeUniqueSons(refTree->getRootNode());
   
   //maintenant qu'on a généré l'arbre, on l'écrit dans le fichier cache
   if (save) writeRefTreeToFile(path);
   
}

int Family::numberOfNodesBetween(Node * ancestor, Node * pnode){
    if(ancestor == pnode) return 0;
    else return(1+numberOfNodesBetween(ancestor, pnode->getFather()));
}



void Family::deleteFromLeavesToBif(Node * pnode){
    if(pnode->hasFather()){
	Node * father = pnode->getFather();
	int nbFils = father->getNumberOfSons();
	if(nbFils == 1 && father->hasFather())
	    deleteFromLeavesToBif(father);
	else
	    father->removeSon(pnode);
        delete pnode;
    } else cout << "deleteFromLeavesToBif : On essaye de supprimer la racine !!!" << endl;
}


Node * Family::removeUniqueSons(Node * localRoot){
    // FONCTION RÉCURSIVE
    unsigned int nbFils = localRoot->getNumberOfSons();
    switch(nbFils){
	case 0: // Cas de base, on est sur une feuille, donc on retourne la feuille
	    return(localRoot);
	case 1: // Cas récursif 1, le nœud a un seul fils
	    return(removeUniqueSons(localRoot->getSon(0)));
	default: // ce nœud n'a pas un fils unique, donc on va attribuer la prochaine bifurcation ou feuille aux fils de localRoot
	    for(unsigned currSonIndex=0 ; currSonIndex < nbFils; currSonIndex++)
		localRoot->setSon(currSonIndex, removeUniqueSons(localRoot->getSon(currSonIndex)));
	    return(localRoot);
    }
}

TreeTemplate<Node> * Family::getSpTree() {
    return(spTree);
}

TreeTemplate<Node> * Family::getTree() {
    return(tree);
}

bool Family::containsSpecie(string specie) {
    return(!(species.find(specie)==species.end()));
}

bool Family::containsSpecies(set<string> speciesList) {
    bool answer = true;
    for(set<string>::iterator specie = speciesList.begin(); answer && specie != speciesList.end(); specie++)
	answer = answer && containsSpecie(*specie);
    
    return(answer);
}

set<string> * Family::getSpecies() {return(&species);}

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


std::set<std::string> Family::nodesToNames(std::set<bpp::Node *> &nodes){
    set<string> answer;
    for(set<Node *>::iterator node = nodes.begin(); node != nodes.end(); node++)
	answer.insert((*node)->getName());
    return(answer);
}

std::string Family::getName(){
    return(_name);
}

std::string Family::mapNodeOnTaxon(bpp::Node & pnode){
    // cette fonction va se servir de l'abre vrai de cette famille pour attribuer un taxon à ce nœud
    // Première étape, obtenir toutes les feuilles de ce nœud
    set<string> leavesNames = getLeavesNamesFromNode(&pnode);
    
    // maintenant, on prend un nom de feuille quelconque dans le set, et on essaye de remonter l'abre vrai de cette famille
    Node * currNode = refTree->getNode(*(leavesNames.begin()));
    
    // il ne peut pas y avoir de doublons dans un set. Donc pour savoir si un taxon mappe sur un nœud, il suffit de tester l'égalité
    while(getLeavesNamesFromNode(currNode)!=leavesNames)
	currNode = currNode->getFather();
    
    //à ce point (sauf échec), currNode est le taxon que nous recherchons. On retourne son nom :
    return(currNode->getName());
}

void Family::loadRefTreeFromFile(string path){
    Newick reader;
    refTree = reader.read(path+"/refTrees/"+_name+".refTree");
}

void Family::loadSpTreeFromFile(string path){
    Newick reader;
    spTree = reader.read(path+"/spTrees/"+_name+".spTree");
}

void Family::writeRefTreeToFile(string path){
    Newick writer;
    writer.write(*refTree,path+"/refTrees/"+_name+".refTree",true);
}

void Family::writeSpTreeToFile(string path){
    Newick writer;
    writer.write(*spTree,path+"/spTrees/"+_name+".spTree",true);
}


void Family::addSequencesNames(Node * currNode)
{
    Node * seqNode = tree->getNode(currNode->getId());
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
    Node * seqNode = tree->getNode(currNode->getId());
    if(seqNode->hasName()&& currNode->isLeaf()) currNode->setName(seqNode->getName());
    else currNode->deleteName();
    vector<Node*> sons = currNode->getSons();
    for(vector<Node *>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++)
	labelWithSequencesNames(*currSon);
}

std::vector< unsigned int >& Family::getUnicityScores()
{
    return(unicityScores);
}

