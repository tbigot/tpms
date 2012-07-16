//
// File: NodeConstraints.cpp
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


#include "NodeConstraints.hpp"
#include "Family.hpp"
#include "Taxon.hpp"
#include "DataBase.hpp"


#include <sstream>
#include <iostream>
#include <cmath>



using namespace std;
using namespace tpms;


namespace tpms{

NodeConstraints::NodeConstraints(DataBase &pRefDB): refDB(pRefDB), ok(true) {
    nature = ANY;
    speciesRestrictionsAsSon = false;
    speciesRestrictions = false;
    type = NODE;
    direct = false;
}

void NodeConstraints::setConstraints(DataBase &pRefDB, string constraintsString, NodeType type){    
    // the params are in the form :
    // D!/+SUBTREE-SPECIES/+LEAVE+SPECIES
    
    this->type = type;
    
    istringstream initSString(constraintsString);
    string paramsString;
    getline(initSString,paramsString,'/');
    
    getline(initSString,subtreeConstraintsString,'/');
    buildAllowedSpecies(allowedSpeciesOnSubtree,subtreeConstraintsString,pRefDB);
    
    getline(initSString,constraintsOnNodeString,'/');
    buildAllowedSpecies(allowedSpeciesOnNode,constraintsOnNodeString,pRefDB);
  


    // dealing with nodeNature
    if(paramsString.find('D') != string::npos) nature = DUPLICATION;
    else if (paramsString.find('S') != string::npos) nature = SPECIATION;
    else nature = ANY;
    
    // dealing with branch type
    if(paramsString.find('!') != string::npos) direct = true;
    
    size_t colonPos = paramsString.find('$');
    if(colonPos != string::npos) {
	string minBootstrapS;
	minBootstrapS =  paramsString.substr(colonPos+1);
		
	istringstream minBootstrapSs(minBootstrapS);
	minBootstrapSs >> minBootstrap;
	
	
    } else minBootstrap = 0;
    
//     readConstraintsString();
//     if(!initString.empty()){
// 	// if authorized species start with a minus (or a ! - compat reasons), implicit +ALL before
// 	if(initString.at(0) == '!')
// 	    initString.at(0) = '-';
// 	if(initString.at(0) == '-')
// 	    initString = "+" + pRefDB.getSpeciesTree()->getRootNode()->getName() + initString;
// 	// case of single species
// 	if(initString.find('+')== string::npos && initString.find('-')== string::npos)
// 	    initString = '+' + initString;
// 	buildAllowedSpecies(allowedSpeciesOnSubtree,initString);
// 	speciesRestrictionsAsSon=true;
//     }
}


void NodeConstraints::buildAllowedSpecies(set<Taxon*>& spset,string spstr, DataBase &pRefDB){
    if(spstr.empty()) return;
    //cout << "Building species list with this string : " << spstr << endl;
    
    // if authorized species start with a minus (or a ! - compat reasons), implicit +ALL before
	if(spstr.at(0) == '!')
	    spstr.at(0) = '-';
	if(spstr.at(0) == '-' && spset.empty())
	    spstr = "+" + pRefDB.getSpeciesTree()->getRootNode()->getName() + spstr;
	// case of single species
	if(spstr.find('+')== string::npos && spstr.find('-')== string::npos)
	    spstr = '+' + spstr;
    
	// cout << "    this string will be used : " << spstr << endl;
    
    // fonction récursive qui analyse la liste des espèces à autoriser.
    // le premier caractère est un + ou un -, il détermine si le taxon est à ajouter ou soustraire
    if(spstr.size() != 0){
	bool toAdd;
	if(spstr.at(0) == '+') toAdd = true; else toAdd = false; //FIXME: il faudrait gérer les erreurs
	unsigned int cmpt = 1;
	while(cmpt < spstr.size() && spstr.at(cmpt) != '+' && spstr.at(cmpt) != '-') cmpt++;
	// quoi qu'il en soit, on vient d'extraire un taxon, qu'on doit ajouter ou supprimer à la liste
	if(toAdd)
	    addTaxon(spset,spstr.substr(1,cmpt-1));
	else
            deleteTaxon(spset,spstr.substr(1,cmpt-1));
	if(cmpt < spstr.size()) buildAllowedSpecies(spset,spstr.substr(cmpt),pRefDB);
    }
}

void NodeConstraints::addTaxon(set<Taxon*>& spset,string taxon)
{
    // très simple, on transforme ce taxon en une liste d'espèces
    // trying to find the taxon in the db
    Taxon* foundTaxon = refDB.nameToTaxon(taxon);
    if (foundTaxon == 00){
	cout << "The taxon “" << taxon << "” in your pattern does not seem to exist in the species tree of your collection." << endl;
        ok = false;
	return;
    }
    set<Taxon*> spList = foundTaxon->getDescendants();
    // et on l'ajoute intégralement à la liste des taxons autorisés, le set ne disposant pas de doublons
    spset.insert(spList.begin(), spList.end());
    
// DEBUG:    
//      cout << "Voici les espèces ajoutées :" << endl;
//      for(set<Taxon*>::iterator ct = spList.begin(); ct != spList.end(); ct++){
// 	cout << (*ct)->getName() << ' ';
//     }
//     cout << endl;
}

bool NodeConstraints::isOk(){
    return(ok);
}

void NodeConstraints::deleteTaxon(set<Taxon*>& spset,string taxon){
    // trying to find the taxon in the db
    Taxon* foundTaxon = refDB.nameToTaxon(taxon);
    if (foundTaxon == 00){
	cout << "The taxon " << taxon << " in your pattern does not seem to appear in the species tree of your collection" << endl;
	return;
    }
    set<Taxon*> spList = foundTaxon->getDescendants();
    for(set<Taxon*>::iterator taxonToDelete = spList.begin(); taxonToDelete != spList.end(); taxonToDelete++){
	set<Taxon*>::iterator found = spset.find(*taxonToDelete);
	if(found != spset.end())
	    spset.erase(found);
    }
}

bool NodeConstraints::allows(Family& family, bpp::Node* node)
{
    if(type == LEAF){
	Taxon* nodeSpecies = family.getSpeciesOfNode(node);
	return(allowedSpeciesOnNode.find(nodeSpecies)!= allowedSpeciesOnNode.end());
    } else {
	// a node is accepted if the nature matches & sufficient bootstrap
	double nodeBootstrap = 0;
	if(node->hasBootstrapValue()) nodeBootstrap = node->getBootstrapValue();
	return((nature==ANY || family.getNatureOf(node) == nature) && nodeBootstrap >= minBootstrap);
    }
}

bool NodeConstraints::allowsAsSon(Family& family, bpp::Node* node)
{
    if(speciesRestrictionsAsSon) // we allow only if the subtree contains only allowed species from father
    {
	// 1st step: getting all the species on the gene tree subtree
	vector<Taxon*> speciesList;
	family.getTaxaOnThisSubtree(node,speciesList);
	// seeing, for each species of the taxaList
	bool allSpAreAllowed = true;
	for(vector<Taxon*>::iterator currSp = speciesList.begin(); currSp != speciesList.end(); currSp++){
	    allSpAreAllowed &= allowedSpeciesOnSubtree.find(*currSp) != allowedSpeciesOnSubtree.end();
	}
	return(allSpAreAllowed);
    }
    return(true);
}

bool NodeConstraints::isDirect(){
    return(direct);
}

std::set< Taxon* >& NodeConstraints::getAllowedSpecies()
{
    return(allowedSpeciesOnNode);
}

bool NodeConstraints::isLeaf()
{
    return(type == LEAF);
}

string NodeConstraints::getStr()
{
    string result;
    switch(nature){
	case DUPLICATION: result= "<DUP>"; break;
	case SPECIATION: result = "<SPE>"; break;
	case ANY: result = "<ANY>"; break;
    }
    if(direct){
        result += " -DIRECT- ";
    }
    ostringstream minBootstrapSs;
    minBootstrapSs << minBootstrap;
    result += "{"+ subtreeConstraintsString +"} " + constraintsOnNodeString + " :" + minBootstrapSs.str();
    return(result);
}

}
