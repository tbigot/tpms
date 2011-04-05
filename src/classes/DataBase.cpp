#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>

#include "TreeTools.hpp"
#include "DataBase.hpp"
#include "Waiter.hpp"
#include "Family.hpp"

#include <Bpp/Phyl/Tree.h>

using namespace std;
using namespace bpp;

DataBase::DataBase(string path): refTreesBuilded(false), speciesTreesBuilded(false) {
	
	// ouverture du fichier
	ifstream RAPfile(path.c_str(), ifstream::in);
	
	// chargement à partir du fichier
	loadFromFile(RAPfile);
	
	// on ferme le fichier proprement :
	RAPfile.close();
	
	//DEBUG OUT
	cout << "\nLoad sucessful!\n" << endl;
	
	std::cout << "Species tree loaded: " << speciesTree->getNumberOfLeaves() << " species." << std::endl;
	
	cout << "Loaded " << families.size() << " families / " << nbFamilies << " expected." << endl;
	if(families.size() != nbFamilies) {
		cout << "WARNING : The file " << path << " is supposed to contain more trees." << endl;
	}
	if(reconciled) cout << "Trees are reconciled."; else cout << "Trees are not reconciled.";
	cout << '\n' << endl;
		
	/*
	cout << "\n---------\nTests sur le premier arbre :" << endl;
	families.begin()->second->genSpTree();
	Newick * afficheur = new Newick(true);
	afficheur->write(*(families.begin()->second->getSpTree()), cout);
	
	cout << "est-ce que STREPTOMYCES COELICOLOR est dans les especes ?" << endl;
	if(families.begin()->second->isInSpecies("STREPTOMYCES COELICOLOR")) cout << "Oui"; else cout << "Non";
	cout << endl;
	
	cout << "est-ce que XXX est dans les especes ?" << endl;
	if(families.begin()->second->isInSpecies("XXX")) cout << "Oui"; else cout << "Non";
	cout << endl;
	*/
	
	
}

bpp::TreeTemplate<bpp::Node> * DataBase::getSpeciesTree() { return(speciesTree); }

void DataBase::iNeedSpeciesTrees(bool verbose, string path, bool generate) {
    if(!speciesTreesBuilded){
	cout << "Generating trees with species names as labels:" << endl;
	Waiter patienteur(&cout, nbFamilies, '#');
	for(vector<Family*>::iterator currFamily = families.begin(); currFamily != families.end(); currFamily++){
	    if(generate) (*currFamily)->genSpTree(false,path);
	    else (*currFamily)->loadSpTreeFromFile(path);
	    patienteur.step();
	} 
    speciesTreesBuilded = true;
    patienteur.drawFinal();
    }
}

void DataBase::iNeedMapping(bool verbose, string path, bool generate) {
    iNeedSpeciesTrees(verbose,path,generate);
    if(!refTreesBuilded){
	Waiter patienteur(&cout, nbFamilies, '#');
	for(vector<Family*>::iterator currFamily = families.begin(); currFamily != families.end(); currFamily++){
	    if(generate) (*currFamily)->genRefTree(false,path);
	    else (*currFamily)->loadRefTreeFromFile(path);
	    patienteur.step();
	}
    refTreesBuilded = true;
    patienteur.drawFinal();
    }
	
}

void DataBase::genUnicityScores() {
    if(!unicityScoresComputed) {
	Waiter patienteur(&cout, nbFamilies, 'o');
	for(vector<Family*>::iterator currFamily = families.begin(); currFamily != families.end(); currFamily++){
	    (*currFamily)->genUnicityScores();
	    patienteur.step();
	}
	unicityScoresComputed = true;
	patienteur.drawFinal();
    }
}



void DataBase::loadFromFile(ifstream & RAPfile) {
	// on extrait les lignes
	
	cout << "Loading preamble..." << flush;
	
	string currLigne;
	// ===PREMIERE LIGNE=== nous donne le nombre d'arbres contenus et leur réconciliation
	getline(RAPfile, currLigne);
	
	// Nombre d'arbres
	int currEspace = currLigne.find(" ");
	stringstream ssNbFam;
	ssNbFam << currLigne.substr(0,currEspace);
	ssNbFam >> nbFamilies;
	
	// Réconciliés ?
	int espacePrecedent = currEspace;
	currEspace = currLigne.find(" ",espacePrecedent+1);
	string sReconciled = currLigne.substr(espacePrecedent+1,currEspace-espacePrecedent-1);
	
	if(sReconciled=="unreconciled") reconciled = false; else reconciled = true;
	
	cout << "    [DONE]" << endl;
	// la seconde ligne contient l'arbre exhaustif des espèces
	cout << "Loading species tree..." << flush;
	
	getline(RAPfile, currLigne);
	loadSpeciesTree(currLigne);
	
	cout << "    [DONE]" << endl;

	// maintenant, prise en charge des familles : blocs suivants :
	cout << "Loading family trees:" << endl;
	
	unsigned int famillesChargees = 0;
	stringstream curPreambule;
	string curNewick;
	bool familleSuivante;
	bool crochetFermant;
	Family * laFamille;
	Waiter patienteur(&cout, nbFamilies, '#');
	while(famillesChargees < nbFamilies) {
		// la première ligne est le nom de la famille, on l'ajoute
		curPreambule.str("");
		familleSuivante = false;
		crochetFermant = false;
		while(!familleSuivante && getline(RAPfile,currLigne)){
			if(crochetFermant) { // on est sur la dernière ligne
				familleSuivante = true;
				curNewick = currLigne + '\n';
			} else {
				curPreambule << currLigne << endl;
				if(currLigne.at(0) == ']') crochetFermant = true;
			}
		}
		laFamille = new Family(&curPreambule,curNewick,this);
		families.push_back(laFamille);
		famillesChargees++;
		patienteur.step();
		
	}
	
	
}

void DataBase::loadSpeciesTree(string newickLine)
{
	// nettoyage des synonymes ET indexation des espèces de la base
	stringstream sNewickPropre;
	char currChar;
	
	bool premierGuillemet = false;
	bool deuxiemeGuillemet = false;
	
	for(unsigned int jogger=0; jogger < newickLine.size(); jogger++) {
		currChar = newickLine.at(jogger);
		
		if(currChar == '"')
			if(premierGuillemet) deuxiemeGuillemet = true;
			else premierGuillemet = true;
				
		else if(currChar == ':') { premierGuillemet = false; deuxiemeGuillemet = false; sNewickPropre << currChar;}
			
		else if((premierGuillemet && !deuxiemeGuillemet) || !premierGuillemet)
			sNewickPropre << currChar;
	}
	
	string snp = sNewickPropre.str();
	speciesTree = tpms::TreeTools::newickToTree(snp,true);
	
	
	//Debug : toutes les espèces :
	//cout << "\n\n\n### Voici la liste de toutes les espèces de l'abre des espèces :" << endl;	
	
	vector<Node *> listeNoeuds = speciesTree->getNodes();
	
	for(vector<Node *>::iterator it = listeNoeuds.begin(); it != listeNoeuds.end(); it++) {
		if((*it)->hasName()) species.insert((*it)->getName());
	}
	
}

vector<Family *> & DataBase::getFamilies() {return(families);}

unsigned int DataBase::getNbFamilies() {
	return(families.size());
}

string DataBase::getParentTaxon(string pTaxon, unsigned int level) {
	Node * pNode = speciesTree->getNode(pTaxon);
	for(unsigned int i =0 ; i < level; i++) pNode = pNode->getFather();
	return(pNode->getName());
}

bool DataBase::taxonExists(string ptax) {
	try{
		speciesTree->getNode(ptax);
		return(true);
	} catch (bpp::NodeNotFoundException e) {
		return(false);
	}
}

int DataBase::nbFamiliesContaining(string pTax){
	// première possibilité : on a déjà fait la recherche
	map<string,int>::iterator nbIt = nbFC.find(pTax);
	if(nbIt != nbFC.end()) return(nbIt->second);
	
	// deuxième cas, il faut faire la recherche et lister les familles
	int result = 0;
	
	// on développe le taxon
	set<string> tMembers = getDescendants(pTax); // contient la liste de toutes les taxons appartenant au taxon pTax
	
	for(vector<Family *>::iterator oneFam = families.begin(); oneFam != families.end(); oneFam++) {
		bool auMoinsUnTaxonPourCetteFamille = false;
		for(set<string>::iterator unTaxon = tMembers.begin(); !auMoinsUnTaxonPourCetteFamille && unTaxon != tMembers.end(); unTaxon++)
			auMoinsUnTaxonPourCetteFamille = auMoinsUnTaxonPourCetteFamille || (*oneFam)->getSpecies()->find(*unTaxon) != (*oneFam)->getSpecies()->end();
		if(auMoinsUnTaxonPourCetteFamille) result++;
	}
	nbFC.insert(pair<string,int>(pTax,result));
	return(result);
}

std::set<std::string> DataBase::getDescendants(string taxon, bool nodesWanted){
	return(getAllNodes(speciesTree->getNode(taxon),nodesWanted));
}

std::set<std::string> DataBase::getDescendants(vector<string> taxaList, bool nodesWanted){
	set<string> speciesList;
	for(vector<string>::iterator taxon = taxaList.begin(); taxon != taxaList.end(); taxon++){
	  set<string> currSpecies = getDescendants(*taxon,nodesWanted);
	  speciesList.insert(currSpecies.begin(), currSpecies.end());
	}
	return(speciesList);
}

// fonction récursive qui retourne le noms de tous les noeuds sous un noeud

std::set<std::string> DataBase::getAllNodes(Node * localRoot, bool nodesWanted){
	std::set<std::string> sons;
	if(localRoot->getNumberOfSons() == 0) { // cas de base : la feuille
		sons.insert(localRoot->getName());
	} else { // cas récursif : on explore les fils
		// mais on entre quand même le nom du nœud
		if(nodesWanted)
		    if(localRoot->hasName())
			sons.insert(localRoot->getName());
		    
		std::set<std::string> desc;
		for(unsigned int i=0; i < localRoot->getNumberOfSons(); i++) {
			desc = getAllNodes(localRoot->getSon(i),nodesWanted);
			sons.insert(desc.begin(),desc.end());
		}
	}
	return(sons);
}

set<string> * DataBase::getSpecies(){
    return(&species);
}



