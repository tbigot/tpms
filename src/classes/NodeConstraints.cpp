#include "NodeConstraints.hpp"
#include "Family.hpp"
#include "Taxon.hpp"
#include "DataBase.hpp"


#include <sstream>
#include <iostream>



using namespace std;
using namespace tpms;


NodeConstraints::NodeConstraints(DataBase &pRefDB): refDB(pRefDB) {
    nature = ANY;
    speciesRestrictionsAsSon = false;
    speciesRestrictions = false;
    type = NODE;
}

void NodeConstraints::setConstraints(DataBase &pRefDB, string constraintsString, Type type){
    initString = constraintsString;
    readConstraintsString();
    if(!initString.empty()){
	// if authorized species start with a minus (or a ! - compat reasons), implicit +ALL before
	if(initString.at(0) == '!')
	    initString.at(0) = '-';
	if(initString.at(0) == '-')
	    initString = "+" + pRefDB.getSpeciesTree()->getRootNode()->getName() + initString;
	buildAllowedSpecies(allowedSpeciesOnSubtree,initString);
	speciesRestrictionsAsSon=true;
    }
}


void NodeConstraints::buildAllowedSpecies(set<Taxon*>& spset,string spstr){
    // fonction récursive qui analyse la liste des espèces à autoriser.
    // le premier caractère est un + ou un -, il détermine si le taxon est à ajouter ou soustraire
    
    bool toAdd;
    if(spstr.at(0) == '+') toAdd = true; else toAdd= false; //FIXME: il faudrait gérer les erreurs
    unsigned int cmpt = 1;
    while(cmpt < spstr.size() && spstr.at(cmpt) != '+' && spstr.at(cmpt) != '-') cmpt++;
    // quoi qu'il en soit, on vient d'extraire un taxon, qu'on doit ajouter ou supprimer à la liste
    if(toAdd)
	addTaxon(spset,spstr.substr(1,cmpt-1));
    else
	deleteTaxon(spset,spstr.substr(1,cmpt-1));
    // si on n'est pas arrivé à la fin, on passe la fin de la chaîne à la fonction
    if(cmpt < spstr.size()) buildAllowedSpecies(spset,spstr.substr(cmpt));
}

void NodeConstraints::addTaxon(set<Taxon*>& spset,string taxon)
{
    // très simple, on transforme ce taxon en une liste d'espèces
    set<Taxon*> spList = refDB.nameToTaxon(taxon)->getDescendants() ;
    // et on l'ajoute intégralement à la liste des taxons autorisés, le set ne disposant pas de doublons
    spset.insert(spList.begin(), spList.end());
}

void NodeConstraints::deleteTaxon(set<Taxon*>& spset,string taxon){
    set<Taxon*> spList = refDB.nameToTaxon(taxon)->getDescendants();
    for(set<Taxon*>::iterator taxonToDelete = spList.begin(); taxonToDelete != spList.end(); taxonToDelete++){
	set<Taxon*>::iterator found = spset.find(*taxonToDelete);
	if(found != spset.end())
	    spset.erase(found);
    }
}

void NodeConstraints::readConstraintsString()
{
    // the params are in the form :
    // D!/+SUBTREE-SPECIES/+LEAVE+SPECIES
    
    istringstream initSString(initString);
    string paramsString;
    getline(initSString,paramsString,'/');
    
    getline(initSString,subtreeConstraintsString,'/');
    buildAllowedSpecies(allowedSpeciesOnSubtree,subtreeConstraintsString);
    
    getline(initSString,constraintsOnNodeString,'/');
    buildAllowedSpecies(allowedSpeciesOnNode,constraintsOnNodeString);


    // dealing with nodeNature
    if(paramsString.find('D') != string::npos) nature = DUPLICATION;
    else if (paramsString.find('S') != string::npos) nature = SPECIATION;
    else nature = ANY;
    
    // dealing with branch type
    if(paramsString.find('!') != string::npos) direct = true;
    
}


bool NodeConstraints::allows(Family& family, bpp::Node* node)
{
    if(type == LEAF){
	Taxon* nodeSpecies = family.getSpeciesOfNode(node);
	return(allowedSpeciesOnNode.find(nodeSpecies)!= allowedSpeciesOnNode.end());
    } else {
	// a node is accepted if the nature matches
	return(nature==ANY || family.getNatureOf(node) == nature);
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
    result += "{"+ subtreeConstraintsString +"} " + constraintsOnNodeString ;
    return(result);
}


