#include "Taxon.hpp"
#include "DataBase.hpp"

using namespace std;
using namespace bpp;
using namespace tpms;

Taxon::Taxon(string name, Node* nodeInSpTree,DataBase &database): name(name), db(database), nodeInSpTree(nodeInSpTree)
{
   
}

void Taxon::genRelations(){
    genDescendantsList(nodeInSpTree);
    genAncestorsList(nodeInSpTree);
}

void Taxon::genAncestorsList(Node* localNode)
{
    if(localNode->hasName())
	ancestors.insert(db.nameToTaxon(localNode->getName()));
    if(localNode->hasFather())
        genAncestorsList(localNode->getFather());
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
