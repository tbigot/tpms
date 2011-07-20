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
	if(localNode->hasFather() && localNode->getFather()->hasName())
	    directAncestor = db.nameToTaxon(localNode->getFather()->getName());
	else
	    directAncestor = 00;
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
    cout << "Nous avons " << ancestors.size() << " ancètres." << endl;
    return(ancestors);
}

std::set< Taxon* >& Taxon::getDescendants()
{
    cout << "Nous avons " << ancestors.size() << " descendants ." << endl;
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
    Taxon *currTaxon = *(taxa.begin());
    while(currTaxon->hasAncestor() && !currTaxon->containsAllTheseSpecies(taxa))
	currTaxon = currTaxon->getDirectAncestor();
    if(currTaxon->containsAllTheseSpecies(taxa)) return(currTaxon);
    else return(00);
}

}