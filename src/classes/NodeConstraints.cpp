#include "NodeConstraints.hpp"
#include "Taxon.hpp"
#include "DataBase.hpp"


#include <sstream>
#include <iostream>


using namespace std;
using namespace tpms;


NodeConstraints::NodeConstraints(DataBase &pRefDB): refDB(pRefDB) {
    nature = ANY;
    speciesRestrictionsFromFather = false;
}

NodeConstraints::NodeConstraints(DataBase &pRefDB, string constraintsString): refDB(pRefDB), initString(constraintsString){
    nature = ANY;
    setConstraints(pRefDB,constraintsString);
}

void NodeConstraints::setConstraints(DataBase &pRefDB, string constraintsString){
    initString = constraintsString;
    catchParams();
    if(!initString.empty()){
	// if authorized species start with a minus (or a ! - compat reasons), implicit +ALL before
	if(initString.at(0) == '!')
	    initString.at(0) = '-';
	if(initString.at(0) == '-')
	    initString = "+" + pRefDB.getSpeciesTree()->getRootNode()->getName() + initString;
	buildAllowedSpecies(fromFatherAllowedSpecies,initString);
	speciesRestrictionsFromFather=true;
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
    if(cmpt < spstr.size()) buildAllowedSpecies(spstr.substr(cmpt));
}

void NodeConstraints::addTaxon(set<Taxon*>& spset,string taxon)
{
    // très simple, on transforme ce taxon en une liste d'espèces
    set<string> spList = refDB.nameToTaxon(taxon)->getDescendants() ;
    // et on l'ajoute intégralement à la liste des taxons autorisés, le set ne disposant pas de doublons
    spset.insert(spList.begin(), spList.end());
}

void NodeConstraints::deleteTaxon(set<Taxon*>& spset,string taxon){
    set<string> spList = refDB.nameToTaxon(taxon)->getDescendants();
    for(set<Taxon*>::iterator taxonToDelete = spList.begin(); taxonToDelete != spList.end(); taxonToDelete++){
	set<Taxon*>::iterator found = spset.find(*taxonToDelete);
	if(found != spset.end())
	    spset.erase(found);
    }
}

bool NodeConstraints::hasSpeciesRestrictionsFromFather(){
    return(speciesRestrictionsFromFather);
}

set<string> & NodeConstraints::getAuthorisedSpecies(){
    return(fromFatherAllowedSpecies);
}

std::string NodeConstraints::getString(){
    return(initString);
}

void NodeConstraints::catchParams()
{
    // we have to find where starts the species list
    size_t speStart = initString.find('+');
    size_t speStart2 = initString.find('-');
    if(speStart2 != string::npos && speStart2 < speStart) speStart = speStart2;
    string paramsString = initString.substr(0,speStart);
    if(speStart != string::npos) initString = initString.substr(speStart);
    else initString = "";
    
    
    // dealing with nodeNature
    if(paramsString.find('D') != string::npos) nature = DUPLICATION;
    else if (paramsString.find('S') != string::npos) nature = SPECIATION;
    else nature = ANY;
    
    // dealing with branch type
    //TODO: handle branch type
    
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

