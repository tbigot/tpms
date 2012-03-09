#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <numeric>
#include <ios>
#include <math.h>

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
		taxa.insert(currTaxon);
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
    this->leaves.insert(leaves.begin(),leaves.end());
    mapping_NodesToTaxa.resize(tree->getNumberOfNodes(),00);
    for(vector<Node *>::iterator leave = leaves.begin(); leave != leaves.end(); leave++){
	std::map<string,Taxon*>::iterator found = mne2tax.find((*leave)->getName());
	if(found==mne2tax.end()){
	    containsUndefinedSequences=true;
	    cout << "Warning: The sequence " << (*leave)->getName() << " has not been associated to a species in family " << name << ". This familly will be ignored." << endl;
	}
	mapping_NodesToTaxa.at((*leave)->getId())=found->second;
    }
}

void Family::doMapping_NodesToNatures(){
    vector<Node *> nodes = tree->getNodes();
    mapping_NodesToNatures.resize(tree->getNumberOfNodes());
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



vector<Node*> Family::getUnicityBestRoots(vector<Node *> nodes){
    // aim of this function: trying all branches as roots, and return those which minimize the unicity score sum.
    vector<float> currScores;
    vector<Node*> bestRoots; // returned at the end
    float bestScoresSum;
    float currScoresSum;
    currScores.resize(tree->getNumberOfNodes()+1);
    // we root at the branches leading to the node
    for(vector<Node *>::iterator currNode = nodes.begin(); currNode != nodes.end(); currNode++) {
        if(!(*currNode)->hasFather()) continue; // root has not father branch
        compute_UnicityScoreOnNode(currScores,*currNode,00,true);
        currScoresSum = accumulate(currScores.begin(),currScores.end(),0);
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
    vector<Node *> nodes = tree->getNodes();
    vector<Node*> bestRoots = getUnicityBestRoots(nodes);
    
    // we have to choose between the possible roots, choosing the longer branch
    Node* bestRoot = *bestRoots.begin();
    double maxBranchSize = bestRoot->getDistanceToFather();
    for(vector<Node*>::iterator currRoot = bestRoots.begin()+1; currRoot != bestRoots.end(); currRoot++){
        if((*currRoot)->getDistanceToFather() > maxBranchSize)
            bestRoot = *currRoot;
    }
    
    // "rerooting" the tree according to this new outgroup :
    tree->newOutGroup(bestRoot);
    
    // after reroot, nodes id might have changed, we must recalculate scores
    doMapping_NodesToUnicityScores();
    
}

vector<Node*> Family::getTaxonomyBestRoots(vector<Node *> nodes){
    if(!doneMapping_NodeToTaxa) doMapping_NodesToTaxa();
    
    vector<Node*> bestRoots; // returned at the end
    unsigned int initDepthSum, currDepthSum;
    unsigned int highestDepthSum = 0;
    // aim of this function: trying all branches as roots, and return those which maximize the taxo level sum
    // we’ll try each node as the root, and see taxonomic mapping
    for(vector<Node *>::iterator currNode = nodes.begin(); currNode != nodes.end(); currNode++) {
	if(!(*currNode)->hasFather()) continue;
	Taxon* thisRoot = mapNodeOnTaxon(true,true,**currNode);
	currDepthSum = thisRoot->getDepth();
	for(vector<Taxon*>::iterator currTaxon = mapping_NodesToTaxa.begin(); currTaxon != mapping_NodesToTaxa.end(); currTaxon++)
	    currDepthSum += (*currTaxon)->getDepth();
	if(bestRoots.empty() || currDepthSum > highestDepthSum){
	    // we've found a first or better root candidate
	    bestRoots.clear();
	    bestRoots.push_back(*currNode);
	    highestDepthSum = currDepthSum;
	} else if(currDepthSum == highestDepthSum)
	    bestRoots.push_back(*currNode);
    }
    
    return(bestRoots);
}

void Family::doRerooting_Taxonomy() {
    vector<Node *> nodes = tree->getNodes();
    
    vector<Node*> bestRoots = getTaxonomyBestRoots(nodes);
    
    // Step 1: we have to choose between the candidate roots, choosing the longer branch
    Node* bestRoot = *bestRoots.begin();
    double maxBranchSize = bestRoot->getDistanceToFather();
    for(vector<Node*>::iterator currRoot = bestRoots.begin()+1; currRoot != bestRoots.end(); currRoot++){
        if((*currRoot)->getDistanceToFather() > maxBranchSize)
            bestRoot = *currRoot;
    }
    
    // Step 2: rerooting at the best place (testing sons, find the best outgroup)
    vector<float> currScores;
    compute_UnicityScoreOnNode(currScores,bestRoot,00);
    vector<Node *> sons = bestRoot->getNeighbors();
    Node* bestOutgroup = 00;
    unsigned int minScore = 0;
    for(vector<Node *>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++){
	if(bestOutgroup == 00 || currScores.at((*currSon)->getId())<minScore ){
	    minScore = currScores.at((*currSon)->getId());
	    bestOutgroup = *currSon;
	}
    }
    
    // "rerooting" the tree according to this new outgroup :
    tree->newOutGroup(bestOutgroup);
    doMapping_NodesToTaxa();
    // DEBUG to see depth sum medification
    // cout << initDepth << "->" << highestDepthSum << endl;
    
}

void Family::doRerooting_UnicityTaxonomy(){
    vector<Node *> nodes = tree->getNodes();
    
    // first using unicity criteria
    vector<Node*> bestRoots = getTaxonomyBestRoots(nodes);
    
    // then, on ex-aequos, using taxonomi criteria
    bestRoots = getTaxonomyBestRoots(bestRoots);
    
    // then, keep the longest branch.
    // we have to choose between the possible roots, choosing the longer branch
    Node* bestRoot = *bestRoots.begin();
    double maxBranchSize = bestRoot->getDistanceToFather();
    for(vector<Node*>::iterator currRoot = bestRoots.begin()+1; currRoot != bestRoots.end(); currRoot++){
        if((*currRoot)->getDistanceToFather() > maxBranchSize)
            bestRoot = *currRoot;
    }
    
    
    // "rerooting" the tree according to this new outgroup :
    tree->newOutGroup(bestRoot);
}

map<Taxon*, unsigned int> Family::compute_UnicityScoreOnNode(vector<float> &scores, Node * node, Node * origin, bool virtualRootOnTheBranch){
    unsigned int id = node->getId();
    map<Taxon*, unsigned int> thisNodeCount;
    
    // step 1: filling the species counts from this node
    vector<Node *> neighbors;
    if(virtualRootOnTheBranch){
	neighbors.push_back(node);
	neighbors.push_back(node->getFather());
    } else
	neighbors = node->getNeighbors();
    for(vector<Node *>::iterator currNeighbor = neighbors.begin(); currNeighbor < neighbors.end(); currNeighbor ++) {
	if(*currNeighbor == origin) continue;
	map<Taxon*,unsigned int> currNeighborCount = compute_UnicityScoreOnNode(scores,*currNeighbor,node,false);
	// adding this new map to the current
	for(map<Taxon*,unsigned int>::iterator currCount = currNeighborCount.begin(); currCount != currNeighborCount.end(); currCount++){	    
	    thisNodeCount[currCount->first] += currCount->second;
	}
    }
    
    if(neighbors.size() == 1){ // leave case
	    thisNodeCount.insert(pair<Taxon*,unsigned int>(mapping_NodesToTaxa.at(id),1));
	    
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
    return(taxa.find(taxon)!=taxa.end());
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

Taxon* Family::mapNodeOnTaxon(bool virtualRootOnTheBranch, bool recordResult,bpp::Node & node, bpp::Node* origin,bool recursive, bpp::Node* ignoredNode){
    unsigned int currNodeID = node.getId();
    vector<Node*> neighbors = node.getNeighbors();
    if(leaves.find(&node) != leaves.end() || !recursive) // BASE CASE: leaf, or don’t continue if asked
    {
	return(mapping_NodesToTaxa.at(currNodeID));
    }
    // dealing with the case: topological leaf but not a real leaf (removed subtree)
    if(neighbors.size() == 1){
	if(recordResult) mapping_NodesToTaxa.at(currNodeID)=00;
	return 00;
    }
    
    set<Taxon*> taxaListOnSons;
    
    // Case of ignored node.
    // if the ignored node is in the sons no need to keep the funciton recursive: nothing has changed
    if(ignoredNode!=00 &&  find(neighbors.begin(),neighbors.end(),ignoredNode) != neighbors.end())
	recursive = false;
    
    
    //collecting taxa on sons (neighbors without the node "origin")
    for(vector<Node *>::iterator currNeighbor = neighbors.begin(); currNeighbor != neighbors.end(); currNeighbor++){
	if(*currNeighbor == origin || *currNeighbor == ignoredNode) continue;
	taxaListOnSons.insert(mapNodeOnTaxon(false,recordResult,**currNeighbor,&node,recursive,ignoredNode));
    }
    Taxon* currTaxon = Taxon::findSmallestCommonTaxon(taxaListOnSons);
    if(recordResult && !virtualRootOnTheBranch) mapping_NodesToTaxa.at(currNodeID) = currTaxon ;
    taxa.insert(currTaxon);
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

std::vector<float>& Family::getUnicityScores()
{
    return(mapping_NodesToUnicityScores);
}


NodeNature Family::getNatureOf(Node* node)
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
    return taxa;   
}

void Family::doMapping_NodesToTaxa(){
    mapNodeOnTaxon(false,true,*tree->getRootNode());
    doneMapping_NodeToTaxa = true;
}

void Family::doMapping_NodesToTaxonomicShift(){
    // initializing the vector
    vector<Node*> nodes = tree->getNodes();
    if(mapping_NodesToTaxonomicShift.empty()) // first using? Initialization
	mapping_NodesToTaxonomicShift.resize(tree->getNumberOfNodes(),0);
    computed_nodesInducingPerturbation.clear();
    // we must fill the vector with zeros, or else, some parts of the algorithm could belive putative transfers still exist!
    fill(mapping_NodesToTaxonomicShift.begin(),mapping_NodesToTaxonomicShift.end(),0);
    
    if(mapping_grandFatherWithoutThisNode.empty())
	mapping_grandFatherWithoutThisNode.resize(tree->getNumberOfNodes());
    for(vector<Node*>::iterator currNode = nodes.begin(); currNode != nodes.end(); currNode++){
	if(computeMappingShiftWithoutTheNode(*currNode) > 0)
	    computed_nodesInducingPerturbation.insert(*currNode);
	
    }
    
}

bool Family::transfersRemaining(){
    bool transfersRemaining = false;
    for(vector<unsigned int>::iterator currPerturbation = mapping_NodesToTaxonomicShift.begin(); !transfersRemaining && currPerturbation != mapping_NodesToTaxonomicShift.end(); currPerturbation++){
	transfersRemaining |= (*currPerturbation > 1);
    }
    return transfersRemaining;
}

void Family::compute_detectTransfers(){
    // doing the first mapping:
    doMapping_NodesToTaxonomicShift();
    // then, we keep transfers that are not included in another
    //FIXME it would be better to confirm bootstrap are sufficient
    //for each potential transfer, (1) does this node include no transfer, (2) if yes, transfer confirmed
    // forbiddenNodes contains all the nodes that we know they can’t be transfers, because they include potential transfers
    while(transfersRemaining()){
	set<Node*> forbiddenNodes;
	for(set<Node *>::iterator node = computed_nodesInducingPerturbation.begin(); node != computed_nodesInducingPerturbation.end(); node++){
	    // (1.0) it this node in forbidden nodes (it would mean that a descendant node has been found under this node)
	    if(forbiddenNodes.find(*node) != forbiddenNodes.end()) continue;
	    // (1.1) it is here possible a descendant of this node is a potential transfer, but has not yet been explored
	    bool transferFoundInDescendants = false;
	    vector<Node*> descendants;
	    tpms::TreeTools::getNodesOfTheSubtree(descendants,*node);
	    for(vector<Node*>::iterator currDescendant = descendants.begin(); !transferFoundInDescendants && currDescendant != descendants.end(); currDescendant++){
		transferFoundInDescendants |= computed_nodesInducingPerturbation.find(*currDescendant) != computed_nodesInducingPerturbation.end();
	    }
	    if(transferFoundInDescendants) continue;
	    
	    
	    // if we’re here, it means that “node” includes no transfer, we can consider it as a potential transfer, and so remove the subtree from the tree
	    // first, we can blacklist all the ancestors (see 1.0)
	    Node* nodeToBlacklist = *node;
	    while(nodeToBlacklist->hasFather()){
		nodeToBlacklist = nodeToBlacklist->getFather();
		forbiddenNodes.insert(nodeToBlacklist);
	    }
	    
	    // then we need to insert the transfer into transfers
	    
	    // only if the bootstrap is enough
	    if((*node)->getBootstrapValue() >= .90 && (*node)->getFather()->getBootstrapValue() >= .90){
		transfer currTransfer;
		currTransfer.donnor = mapping_NodesToTaxa.at((*node)->getId());
		currTransfer.receiver = mapping_grandFatherWithoutThisNode.at((*node)->getId());
		currTransfer.perturbationIndex = mapping_NodesToTaxonomicShift.at((*node)->getId());
		computed_detectedTransfers.push_back(currTransfer);
	    }
	    
	    // then, we prune the subtree
	    (*node)->getFather()->removeSon(*node);
	    std::vector<Node*> nodesToDelete;
	    tpms::TreeTools::getNodesOfTheSubtree(nodesToDelete,*node);
	    for(vector<Node*>::iterator currNode = nodesToDelete.begin(); currNode != nodesToDelete.end(); currNode++){
		set<Taxon*>::iterator taxonToDelete = taxa.find(mapping_NodesToTaxa.at((*currNode)->getId()));
		if(taxonToDelete != taxa.end())
		    taxa.erase(taxonToDelete);
		delete *currNode;
	    }
	    delete *node;
	    
	}
	doMapping_NodesToTaxa();
	doMapping_NodesToTaxonomicShift();
// 	cout << computed_detectedTransfers.size() << ',' << flush;
    }
    
    // now, the gene tree should be congruent with the species tree.
    // outputing the result
    
//     if(computed_detectedTransfers.size() > 1){
//     cout << name << " : ";
//     cout << computed_detectedTransfers.size() << " transfers detected. Original number of leaves : " << mne2tax.size() << ". Final number of leaves : " << tree->getNumberOfLeaves() << endl;}
    
    for(vector<transfer>::iterator currTransfer = computed_detectedTransfers.begin(); currTransfer != computed_detectedTransfers.end(); currTransfer++){
	results << name << ',' << currTransfer->perturbationIndex << ',' << currTransfer->receiver->getName() << ',' << currTransfer->donnor->getName() << endl;
    }
    
}


void Family::threadedWork_launchJobs(std::vector<Family *> families, void (Family::*function)(), unsigned int nbThreads, ostream *output){
    unsigned int nbFamilies = families.size();
    Waiter progressbar(&cout, nbFamilies, '#');
    boost::mutex progressbarMutex;
    boost::mutex outputMutex;
    boost::thread_group tg;
    unsigned int blockSize = nbFamilies / nbThreads;
    cout << "Multithreaded operation. Number of threads: " << nbThreads << ". Lot size : " << blockSize << endl;
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
	    progressbar->step(50);
	    progressbarMutex->unlock();
	}
    }
}


bool Family::getContainsUndefinedSequences(){
    return(containsUndefinedSequences);
}

unsigned int Family::computeMappingShiftWithoutTheNode(Node* node){
    if(!(node->hasFather() && node->getFather()->hasFather()) || node->getNeighbors().size() <= 2)
	return(0);
    Node* grandFatherNode = node->getFather()->getFather();
    Node* grandFatherNodeFather = 00;
    if(grandFatherNode->hasFather()) grandFatherNodeFather = grandFatherNode->getFather();
    Taxon* initialGfTaxon = mapping_NodesToTaxa.at(grandFatherNode->getId());
    Taxon* newTaxon = mapNodeOnTaxon(false,false,*grandFatherNode,grandFatherNodeFather,true,node);
    
    if(newTaxon == 00 || initialGfTaxon == 00) return 0;
    
    mapping_NodesToTaxonomicShift[node->getId()] = Taxon::computeRelativeDepthDifference(initialGfTaxon,newTaxon,&taxa);
    mapping_grandFatherWithoutThisNode[node->getId()] = newTaxon;
    return(mapping_NodesToTaxonomicShift[node->getId()]);
    
    //DEBUG: printing result
//     unsigned int gain = mapping_NodesToTaxonomicShift[node->getId()];
//     if(gain > 1)
// 	cout << newTaxon->getName() << "->" << (mapping_NodesToTaxa[node->getId()])->getName() << " GAIN: " << gain << endl;
//     
}


ostringstream& Family::threadedWork_getResults()
{
    return(results);
}


}
