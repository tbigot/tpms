#include "NodeConstraints.hpp"
#include "DataBase.hpp"


#include <sstream>
#include <iostream>

using namespace std;

NodeConstraints::NodeConstraints(DataBase &pRefDB): refDB(pRefDB) {
    nature = ANY;
    speciesRestrictions = false;
}

NodeConstraints::NodeConstraints(DataBase &pRefDB, string constraintsString): refDB(pRefDB), initString(constraintsString){
    nature = ANY;
    setConstraints(pRefDB,constraintsString);
}

void NodeConstraints::setConstraints(DataBase &pRefDB, string constraintsString){
    initString = constraintsString;
    cachParams();
    if(!initString.empty()){
	// if authorized species start with a minus, implicit +ALL before
	if(initString.at(0) == '-')
	    initString = "+" + pRefDB.getSpeciesTree()->getRootNode()->getName() + initString;
	buildAuthorisedSpecies(initString);
	speciesRestrictions=true;   
    }
}


void NodeConstraints::buildAuthorisedSpecies(string speciesList){
    // fonction récursive qui analyse la liste des espèces à autoriser.
    // le premier caractère est un + ou un -, il détermine si le taxon est à ajouter ou soustraire
    
    bool toAdd;
    if(speciesList.at(0) == '+') toAdd = true; else toAdd= false; //FIXME: il faudrait gérer les erreurs
    unsigned int cmpt = 1;
    while(cmpt < speciesList.size() && speciesList.at(cmpt) != '+' && speciesList.at(cmpt) != '-') cmpt++;
    // quoi qu'il en soit, on vient d'extraire un taxon, qu'on doit ajouter ou supprimer à la liste
    if(toAdd)
	addTaxon(speciesList.substr(1,cmpt-1));
    else
	deleteTaxon(speciesList.substr(1,cmpt-1));
    // si on n'est pas arrivé à la fin, on passe la fin de la chaîne à la fonction
    if(cmpt < speciesList.size()) buildAuthorisedSpecies(speciesList.substr(cmpt));
}

void NodeConstraints::addTaxon(string taxon)
{
    // très simple, on transforme ce taxon en une liste d'espèces
    set<string> spList = refDB.getDescendants(taxon);
    // et on l'ajoute intégralement à la liste des taxons autorisés, le set ne disposant pas de doublons
    authorisedSpecies.insert(spList.begin(), spList.end());
}

void NodeConstraints::deleteTaxon(string taxon){
    set<string> spList = refDB.getDescendants(taxon);
    for(set<string>::iterator taxonToDelete = spList.begin(); taxonToDelete != spList.end(); taxonToDelete++){
	set<string>::iterator found = authorisedSpecies.find(*taxonToDelete);
	if(found != authorisedSpecies.end())
	    authorisedSpecies.erase(found);
    }
}

bool NodeConstraints::hasSpeciesRestrictions(){
    return(speciesRestrictions);
}

set<string> & NodeConstraints::getAuthorisedSpecies(){
    return(authorisedSpecies);
}

std::string NodeConstraints::getString(){
    return(initString);
}

void NodeConstraints::cachParams()
{
    // we have to find where starts the species list
    size_t speStart = initString.find('+');
    size_t speStart2 = initString.find('-');
    if(speStart2 != string::npos && speStart2 < speStart) speStart = speStart2;
    string paramsString = initString.substr(0,speStart);
    if(speStart != string::npos) initString = initString.substr(speStart);
    else initString = "";
    
    if(paramsString.find('D') != string::npos) nature = DUPLICATION;
    else if (paramsString.find('S') != string::npos) nature = SPECIATION;
    else nature = ANY;
    
}

NodeConstraints::NodeNature NodeConstraints::getNature()
{
    return(nature);
}

bool NodeConstraints::isAuthorizedNature(NodeConstraints::NodeNature nature)
{
    if(this->nature == ANY || this->nature == nature) return(true);
    else return(false);
}

