#include "Taxon.hpp"
#include "DataBase.hpp"

using namespace std;
using namespace bpp;
using namespace tpms;

Taxon::Taxon(string name, Node* nodeInSpTree,DataBase &database): _name(name)
{
    genDescendantsList(nodeInSpTree,database);
    genAncestorsList(nodeInSpTree,database);
}

void Taxon::genAncestorsList(Node* localNode, DataBase &database)
{
    
}

void Taxon::genDescendantsList(Node* localNode, DataBase &database)
{
    if(localNode->hasName())
	_descendants.insert(database.nameToTaxon(localNode->getName()));

    vector<Node*> sons = localNode->getSons();
    for(vector<Node*>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++) {
	genDescendantsList(*currSon,database);
    }
}

string Taxon::getName()
{
    return(_name);
}
