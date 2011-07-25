#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <numeric>
#include <ios>

#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>

#include "Family.hpp"
#include "Waiter.hpp"
#include "TreeTools.hpp"
#include "Taxon.hpp"
#include "DataBase.hpp"


using namespace std;
using namespace bpp;
using namespace tpms;


namespace tpms{
    Family::Family(stringstream* sIntro, string* sNewick, DataBase* dbp): db(dbp), preamble(sIntro), newick(sNewick), containsUndefinedSequences(false) {
	
    }
    
    
    void Family::initialize(){
	// extraction du nom de la famille
	string currLigne;
	getline(*preamble,currLigne);
	istringstream sFamilleName(currLigne);
	sFamilleName >> name;
	// extraction des synonymes
	getline(*preamble, currLigne);
	if(currLigne.at(0) != '[') std::cerr << "Bracket (newick comments beginning) expected at family " << name << endl;
	   
	   // cas du fichier avec deux arbres en plus (concerne les arbres réconciliés
	getline(*preamble, currLigne);
	if(currLigne.at(0) == '#') { getline(*preamble, currLigne);getline(*preamble, currLigne);getline(*preamble, currLigne);} // on saute ces lignes
	
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
	    
	    
	    Taxon * currTaxon = db->nameToTaxon(cleanSpeciesName.str());
	    if(currTaxon != 00){
		mne2tax.insert(pair<string,Taxon *>(mnemonique.str(),currTaxon ));
		species.insert(currTaxon);
	    } else {
		cout << "o Unable to find this specie in the species tree:" << cleanSpeciesName.str() << endl;
	    }
	    
	    getline(*preamble, currLigne);
	}
	
	// on fabrique l'arbre
	//istringstream ssNewick(sNewick);
	//Newick newickReader(false);
	//tree = newickReader.read(ssNewick);
	
	tree = tpms::TreeTools::newickToTree(*newick,false);
	
	/*if(tree->getNumberOfNodes() != tree2->getNumberOfNodes()) cout << "Incoherence famille " << _name << " bpp=" << tree->getNumberOfNodes() << " tpms=" <<tree2->getNumberOfNodes() << endl  ;
	 */
	
	delete newick;
	delete preamble;
	
}

void Family::doMapping_LeavesToSpecies(){
    vector<Node *> leaves = tree->getLeaves();
    mapping_NodesToTaxa.resize(tree->getNodes().size());
    for(vector<Node *>::iterator leave = leaves.begin(); leave != leaves.end(); leave++){
	std::map<string,Taxon*>::iterator found = mne2tax.find((*leave)->getName());
	if(found==mne2tax.end()){
	    containsUndefinedSequences=true;
	    cout << "Danger: The sequence " << (*leave)->getName() << " has not been associated to a species in family " << name << "." << endl;
	}
	mapping_NodesToTaxa.at((*leave)->getId())=found->second;
    }
}

void Family::doMapping_NodesToNatures(){
    vector<Node *> nodes = tree->getNodes();
    mapping_NodesToNatures.resize(nodes.size());
    for(vector<Node *>::iterator node = nodes.begin(); node != nodes.end(); node++){
	if((*node)->hasName() && (*node)->getName().at(0) == '#')
	    mapping_NodesToNatures.at((*node)->getId())=DUPLICATION;
	else mapping_NodesToNatures.at((*node)->getId())=SPECIATION;
    }
}


void Family::doMapping_NodesToUnicityScores() {
    mapping_NodesToUnicityScores.resize(tree->getNumberOfNodes());
    compute_UnicityScoreOnNode(mapping_NodesToUnicityScores,tree->getRootNode(),00); // 00 because no origin, we are at the real root
}

void Family::doMapping_NodesToBestUnicityScores() {
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
	compute_UnicityScoreOnNode(currScores,*currNode,00);
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
    doMapping_NodesToUnicityScores();
    
}


void Family::doMapping_NodesToLowestTaxa() {
    unsigned int initDepthSum, currDepthSum;
    unsigned int highestDepthSum = 0;
    Node * bestRoot = 00;
    vector<Node *> nodes = tree->getNodes();
    
    // we’ll try each node as the root, and see taxonomic mapping
    for(vector<Node *>::iterator currNode = nodes.begin(); currNode != nodes.end(); currNode++) {
	if((*currNode)->isLeaf()) continue; // we don’t want to root by a leaf
	mapNodeOnTaxon(true,**currNode);
	currDepthSum = 0;
	for(vector<Taxon*>::iterator currTaxon = mapping_NodesToTaxa.begin(); currTaxon != mapping_NodesToTaxa.end(); currTaxon++)
	    currDepthSum += (*currTaxon)->getDepth();
	if(bestRoot == 00 || currDepthSum > highestDepthSum){
	    // we've found a first or better root candidate
	    bestRoot = *currNode;
	    highestDepthSum = currDepthSum;
	}
    }
    
    tree->rootAt(bestRoot);
    doMapping_NodesToTaxa();
    // DEBUG to see depth sum medification
    // cout << initDepth << "->" << highestDepthSum << endl;
    
}

map<Taxon*, unsigned int> Family::compute_UnicityScoreOnNode(vector<unsigned int> &scores, Node * node, Node * origin){
    unsigned int id = node->getId();
    map<Taxon*, unsigned int> thisNodeCount;
    
    // step 1 : this node count
    
    vector<Node *> neighbors = node->getNeighbors();
    for(vector<Node *>::iterator currSon = neighbors.begin(); currSon < neighbors.end(); currSon ++) {
	if(*currSon == origin) continue;
	   map<Taxon*,unsigned int> currSonCount = compute_UnicityScoreOnNode(scores,*currSon,node);
	// adding this new map to the current
	   for(map<Taxon*,unsigned int>::iterator currCount = currSonCount.begin(); currCount != currSonCount.end(); currCount++){	    
	       thisNodeCount[currCount->first] += currCount->second;
	   }
    }
    
    if(neighbors.size() == 1){ // leave case
	    thisNodeCount.insert(pair<Taxon*,unsigned int>(mapping_NodesToTaxa.at(id),1));
	    
    }
    
    // step 2 : this node score computation
    unsigned int score = 1;
    
    for(map<Taxon*,unsigned int>::iterator currCount = thisNodeCount.begin(); currCount != thisNodeCount.end(); currCount++){
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

TreeTemplate<Node> * Family::getTree() {
    return(tree);
}

bool Family::containsSpecie(Taxon* taxon) {
    return(species.find(taxon)!=species.end());
}


Taxon* Family::getSpeciesOfNode(Node* node){
    return(mapping_NodesToTaxa.at(node->getId()));
}

// bool Family::containsSpecies(set<string> speciesList) {
	    //     bool answer = true;
	    //     for(set<string>::iterator specie = speciesList.begin(); answer && specie != speciesList.end(); specie++)
// 	answer = answer && containsSpecie(*specie);
//     
//     return(answer);
// }


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
    return(name);
}

Taxon* Family::mapNodeOnTaxon(bool recordResult,bpp::Node & node, bpp::Node* origin,bool recursive, bpp::Node* ignoredNode){
    unsigned int currNodeID = node.getId();
    vector<Node*> neighbors = node.getNeighbors();
    if(neighbors.size() == 1 || !recursive) // BASE CASE: leaf, or don’t continue if asked
    {
	return(mapping_NodesToTaxa.at(currNodeID));
    }
    set<Taxon*> taxaListOnSons;
 
    // if the ignored node is in the sons no need to keep the funciton recursive: nothing has changed
    if(ignoredNode!=00 &&  find(neighbors.begin(),neighbors.end(),ignoredNode) != neighbors.end())
	recursive = false;
    
    //collecting taxa on sons (neighbors without origin)
    for(vector<Node *>::iterator currNeighbor = neighbors.begin(); currNeighbor != neighbors.end(); currNeighbor++){
	if(*currNeighbor == origin || *currNeighbor == ignoredNode) continue;
	taxaListOnSons.insert(mapNodeOnTaxon(recordResult,**currNeighbor,&node,recursive,ignoredNode));
    }
    Taxon* currTaxon = Taxon::findSmallestCommonTaxon(taxaListOnSons);
    if(recordResult) mapping_NodesToTaxa.at(currNodeID) = currTaxon ;
    return(currTaxon);
    
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
    return(mapping_NodesToUnicityScores);
}


Family::NodeNature Family::getNatureOf(Node* node)
{
    return(mapping_NodesToNatures.at(node->getId()));
}

void Family::getTaxaOnThisSubtree(Node* node, std::vector< Taxon* >& speciesList)
{
    vector<Node *> sons = node->getSons();
    for(vector<Node*>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++){
	getTaxaOnThisSubtree(*currSon,speciesList);
    }
    if(sons.size() == 0) // leave case, return the species of this leave
	speciesList.push_back(mapping_NodesToTaxa.at(node->getId()));
}

set<Taxon *> &Family::getSpecies(){
    return species;   
}

void Family::doMapping_NodesToTaxa(){
    mapNodeOnTaxon(true,*tree->getRootNode());
}

void Family::doMapping_NodesToTaxonomicShift(){
    // initializing the vector
    vector<Node*> nodes = tree->getNodes();
    mapping_NodesToTaxonomicShift.resize(nodes.size());
    for(vector<Node*>::iterator currNode = nodes.begin(); currNode != nodes.end(); currNode++){
	computeMappingShiftWithoutTheNode(*currNode);
    }
    
}


void Family::threadedWork_launchJobs(std::vector<Family *> families, void (Family::*function)(), unsigned int nbThreads){
    unsigned int nbFamilies = families.size();
    Waiter progressbar(&cout, nbFamilies, '#');
    boost::mutex progressbarMutex;
    boost::thread_group tg;
    unsigned int blockSize = nbFamilies / nbThreads;
    cout << "Multithreaded operation. Number of threads: " << nbThreads << ". Lot size : " << blockSize << endl;
    for(unsigned int currThreadIndex = 0; currThreadIndex < nbThreads; currThreadIndex++){
	vector<Family*>::iterator currPartBegin;
	vector<Family*>::const_iterator currPartEnd;
	currPartBegin = families.begin() + (blockSize*currThreadIndex);
	if(currThreadIndex+1 != nbThreads) currPartEnd = families.begin() + (blockSize*(currThreadIndex+1)); else currPartEnd = families.end();
	   boost::thread *currThread = new boost::thread(Family::threadedWork_oneThread,function,&progressbar,&progressbarMutex,currPartBegin,currPartEnd);
	tg.add_thread(currThread);
    }
    tg.join_all();
    progressbar.drawFinal();
}

void Family::threadedWork_oneThread(void(Family::*function)(),Waiter *progressbar, boost::mutex *progressbarMutex, std::vector<tpms::Family*>::iterator &currPartBegin, std::vector<tpms::Family*>::const_iterator &currPartEnd){
    unsigned int waiterUpdate = 0;
    for(vector<Family*>::iterator currFamily = currPartBegin; currFamily != currPartEnd; currFamily++){
	if(!(*currFamily)->getContainsUndefinedSequences())
	    (*currFamily->*function)();
	waiterUpdate++;
	if(waiterUpdate == 50){
	    waiterUpdate = 0;
	    progressbarMutex->lock();
	    progressbar->step(50);
	    progressbarMutex->unlock();
	}
    }
}


bool Family::getContainsUndefinedSequences(){
    return(containsUndefinedSequences);
}

unsigned int Family::computeMappingShiftWithoutTheNode(Node* node){
    if(!(node->hasFather() && node->getFather()->hasFather()))
	return(0);
    Node* grandFatherNode = node->getFather()->getFather();
    Node* grandFatherNodeFather = 00;
    if(grandFatherNode->hasFather()) grandFatherNodeFather = grandFatherNode->getFather();
    Taxon* initialGfTaxon = mapping_NodesToTaxa[grandFatherNode->getId()];
    Taxon* newTaxon = mapNodeOnTaxon(false,*grandFatherNode,grandFatherNodeFather,true,node);
    
    mapping_NodesToTaxonomicShift[node->getId()] = newTaxon->getDepth() - initialGfTaxon->getDepth();
    
    //DEBUG: printing result
    //cout << "Initial Taxon Assignation : " << initialGfTaxon->getName() << ". Final : " << newTaxon->getName() << ". Difference : " << mapping_NodesToTaxonomicShift[node->getId()] << endl;
    
}


}
