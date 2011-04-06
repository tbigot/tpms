#include "Taxon.hpp"
#include "DataBase.hpp"

using namespace std;
using namespace bpp;
using namespace tpms;

Taxon::Taxon(string name, Node* nodeInSpTree,DataBase &database): name(name)
{
    genDescendantsList(nodeInSpTree,database);
    genAncestorsList(nodeInSpTree,database);
}

void Taxon::genAncestorsList(Node* localNode, DataBase &database)
{
    if(localNode->hasName())
	ancestors.insert(database.nameToTaxon(localNode->getName()));
    if(localNode->hasFather())
        genAncestorsList(localNode->getFather(),database);
}

void Taxon::genDescendantsList(Node* localNode, DataBase &database)
{
    if(localNode->hasName())
	descendants.insert(database.nameToTaxon(localNode->getName()));

    vector<Node*> sons = localNode->getSons();
    for(vector<Node*>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++) {
	genDescendantsList(*currSon,database);
    }
}

string Taxon::getName()
{
    return(name);
}

std::set< Taxon* >& Taxon::getAncestors()
{
    return(ancestors);
}

std::set< Taxon* >& Taxon::getDescendants()
{
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
