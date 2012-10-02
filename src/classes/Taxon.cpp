//
// File: Taxon.cpp
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

string Taxon::getTaxonomy(){
    ostringstream result;
    bool sep=false;
    for(set<Taxon*>::reverse_iterator currAncestor = ancestors.rbegin(); currAncestor != ancestors.rend(); currAncestor++){
        result << (sep? "/": "") << (*currAncestor)->getName();
        sep=true;
    }
    return(result.str());
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
	if(taxaList->find(descendant->getDirectAncestor()) != taxaList->end()) toAdd = 1;
	return(toAdd + computeRelativeDepthDifference(ancestor,descendant->getDirectAncestor(), taxaList));
	
    }
    return 0;
}



}