#include "Taxon.hpp"
#include "DataBase.hpp"

using namespace std;
using namespace bpp;
using namespace tpms;


namespace tpms{
Taxon::Taxon(string name, Node* nodeInSpTree,DataBase &database): name(name), db(database), nodeInSpTree(nodeInSpTree)
{
   
}

void Taxon::genRelations(){
    genDescendantsList(nodeInSpTree);
    genAncestorsList(00);
    depth = ancestors.size();
}


void Taxon::genAncestorsList(Node* localNode)
{
    bool firstTime = false;
    
    if(localNode == 00){
	firstTime = true;
	localNode = nodeInSpTree;
    }
    
    if(localNode->hasName())
	ancestors.insert(db.nameToTaxon(localNode->getName()));
    if(localNode->hasFather())
        genAncestorsList(localNode->getFather());
    
    if(firstTime){
	Node* currNode = localNode;
	do{
	if(currNode->hasFather() && currNode->getFather()->hasName()){
	    directAncestor = db.nameToTaxon(currNode->getFather()->getName());
	    currNode = currNode->getFather();
	} else
	    directAncestor = 00;
	} while (directAncestor != 00 && directAncestor->getName() == name);
    }
}

void Taxon::genDescendantsList(Node* localNode)
{
    if(localNode->hasName())
	descendants.insert(db.nameToTaxon(localNode->getName()));

    vector<Node*> sons = localNode->getSons();
    for(vector<Node*>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++) {
	genDescendantsList(*currSon);
    }
}

string Taxon::getName()
{
    return(name);
}

std::set< Taxon* >& Taxon::getAncestors()
{
    cout << name << " has " << ancestors.size() << " ancestors." << endl;
    return(ancestors);
}

std::set< Taxon* >& Taxon::getDescendants()
{
    cout << name << " has " << ancestors.size() << " descendants." << endl;
    return(descendants);
}



bool Taxon::belongsTo(Taxon* ancestor)
{
    return(ancestors.find(ancestor)!=ancestors.end());
}

bool Taxon::contains(Taxon* descendant)
{
    return(descendants.find(descendant)!=descendants.end());
}

bool Taxon::hasAncestor()
{
    return(!ancestors.empty());
}

Taxon* Taxon::getDirectAncestor()
{
    return(directAncestor);
}

bool Taxon::hasDirectAncestor(){
    return(directAncestor != 00);
}

bool Taxon::containsAllTheseSpecies(std::set< Taxon* > species)
{
    bool allSpeciesAreInThisTaxon = true;
    for(std::set<Taxon *>::iterator currSp = species.begin(); allSpeciesAreInThisTaxon && currSp != species.end(); currSp++){
	allSpeciesAreInThisTaxon &= contains(*currSp);
    }
    return(allSpeciesAreInThisTaxon);
}


tpms::Taxon* Taxon::findSmallestCommonTaxon(set<Taxon*> taxa){
    //starting with the first taxon
    set<Taxon*>::iterator currTaxonIt = taxa.begin();
    // to manage unresolved nodes : propagates unresolution to nodes
    while(*currTaxonIt == 00 && currTaxonIt != taxa.end()) currTaxonIt++;
    if(*currTaxonIt == 00 || currTaxonIt==taxa.end()) return 00;
    
    Taxon *currTaxon = *currTaxonIt;
   
    while(currTaxon->hasDirectAncestor() && !currTaxon->containsAllTheseSpecies(taxa)){
	currTaxon = currTaxon->getDirectAncestor();
    }
    if(currTaxon->containsAllTheseSpecies(taxa)) return(currTaxon);
    else return(00);
}


unsigned int Taxon::getDepth(){
   return(depth);    
}

unsigned int Taxon::computeRelativeDepthDifference(Taxon* ancestor, Taxon* descendant, set<Taxon*>* taxaList)
{
    if(ancestor == descendant)
	return 0;
    if(descendant->hasDirectAncestor()){
	unsigned int toAdd = 0;
	if(taxaList->find(descendant) != taxaList->end()) toAdd = 1;
	return(toAdd + computeRelativeDepthDifference(ancestor,descendant->getDirectAncestor(), taxaList));
	
    }
    return 0;
}



}