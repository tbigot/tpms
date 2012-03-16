#include "Taxon.hpp"
#include "DataBase.hpp"

using namespace std;
using namespace bpp;
using namespace tpms;


namespace tpms{
Taxon::Taxon(string name, Node* nodeInSpTree,DataBase &database): name(name), db(database), nodeInSpTree(nodeInSpTree)
{
   if(name.empty()) {
       ostringstream ssid;
       ssid << nodeInSpTree->getId();
       this->name="<unnamed taxon " + ssid.str() + ">";
   } 
}

void Taxon::genRelations(){
    genDescendantsList(nodeInSpTree);
    genAncestorsList(00);
    depth = ancestors.size()-1;
}


void Taxon::genAncestorsList(Node* localNode)
{
    bool firstTime = false;
    
    if(localNode == 00){
	firstTime = true;
	localNode = nodeInSpTree;
    }
    
    ancestors.insert(db.nodeIdToTaxon(localNode->getId())); // base case
    if(localNode->hasFather()) // recursive case
        genAncestorsList(localNode->getFather());
    
    if(firstTime){ // initialization of direct ancestor
	if(localNode->hasFather())
	    directAncestor = db.nodeIdToTaxon(localNode->getFather()->getId());
	else
	    directAncestor = 00;
    }
}

void Taxon::genDescendantsList(Node* localNode)
{

    descendants.insert(db.nodeIdToTaxon(localNode->getId()));

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
    // cout << name << " has " << ancestors.size() << " ancestors." << endl;
    return(ancestors);
}

std::set< Taxon* >& Taxon::getDescendants()
{
    // cout << name << " has " << descendants.size() << " descendants:" << endl;
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
	allSpeciesAreInThisTaxon &= (contains(*currSp) || *currSp == 00);
    }
    return(allSpeciesAreInThisTaxon);
}


tpms::Taxon* Taxon::findLCA(set<Taxon*> taxa){
    set<Taxon*>::iterator nullFound = taxa.find(00);
    if(nullFound != taxa.end())
	cerr << "Cannot find LCA from null species." << endl;
    
    //starting with the first taxon
    Taxon *currTaxon = *taxa.begin();
   
    while(currTaxon->hasDirectAncestor() && !currTaxon->containsAllTheseSpecies(taxa)){
	currTaxon = currTaxon->getDirectAncestor();
    }
    if(currTaxon->containsAllTheseSpecies(taxa)) return(currTaxon);
    else {
	cout << "Impossible de trouver LCA pour " ;
	for(set<Taxon*>::iterator currTx = taxa.begin(); currTx != taxa.end(); currTx++){
	    cout << " " << (*currTx)->getName() << " :\nD= ";
	    for(set<Taxon*>::iterator currDesc = (*currTx)->getDescendants().begin(); currDesc != (*currTx)->getDescendants().end(); currDesc++){
		cout << "- " << (*currDesc)->getName(); 
	    }
	    cout << "\nA= " ;
	    for(set<Taxon*>::iterator currDesc = (*currTx)->getAncestors().begin(); currDesc != (*currTx)->getAncestors().end(); currDesc++){
		cout << "- '" << (*currDesc)->getName() << "' ***"; 
	    }
	    cout << endl;
	}
	cout << endl;
	return(00);}
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