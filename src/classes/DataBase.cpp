#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string.hpp>


#include "TreeTools.hpp"
#include "DataBase.hpp"
#include "Waiter.hpp"
#include "Family.hpp"
#include "Taxon.hpp"

#include <Bpp/Phyl/Tree.h>


using namespace std;
using namespace bpp;
using namespace tpms;

namespace fs = boost::filesystem;


namespace tpms{

DataBase::DataBase(string path, unsigned int nbThreads): mappingDone_NodesToTaxa(false), mappingDone_LeavesToSpecies(false), nbThreads(nbThreads) {
	
	// checking this file exists
	fs::path dbFile(path);
	if(!fs::exists(dbFile)) {
	    cout << "The file "<< path << " does not exists, exiting."<< endl;
	    exit(1);
	}
    
    
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

// void DataBase::doFamiliesMapping_LeavesToSpecies() {
//     if(!mappingDone_LeavesToSpecies){
// 	cout << "Mapping leaves to Species:" << endl;
// 	Waiter patienteur(&cout, nbFamilies, '#');
// 	for(vector<Family*>::iterator currFamily = families.begin(); currFamily != families.end(); currFamily++){
// 	    (*currFamily)->doMapping_LeavesToSpecies();
// 	    patienteur.step();
// 	} 
//     mappingDone_LeavesToSpecies = true;
//     patienteur.drawFinal();
//     }
// }

// 
void DataBase::doFamiliesMapping_LeavesToSpecies() {
    // multithreading version
    if(!mappingDone_LeavesToSpecies){
	cout << "Mapping leaves to Species:" << endl;
	Waiter patienteur(&cout, nbFamilies, '#');
	boost::thread_group tg;
	unsigned int blockSize = families.size() / nbThreads;
	cout << "Multithreaded operation. Lot size : " << blockSize << endl;
	for(unsigned int currThreadIndex = 0; currThreadIndex < nbThreads; currThreadIndex++){
	    vector<Family*>::iterator familyBegin, familyEnd;
	    familyBegin = families.begin() + (blockSize*currThreadIndex);
	    if(currThreadIndex+1 != nbThreads) familyEnd = families.begin() + (blockSize*(currThreadIndex+1)); else familyEnd = families.end();
	    boost::thread *currThread = new boost::thread(doFamiliesMapping_LeavesToSpecies_oneThread,&patienteur,&waiterMutex,familyBegin,familyEnd);
	    tg.add_thread(currThread);
	}
	tg.join_all();
	mappingDone_LeavesToSpecies = true;
	patienteur.drawFinal();
    }
}

void DataBase::doFamiliesMapping_LeavesToSpecies_oneThread(Waiter *waiter, boost::mutex *waiterMutex, vector<Family*>::iterator &familiesBegin, vector<Family*>::iterator &familiesEnd) {
	for(vector<Family*>::iterator currFamily = familiesBegin; currFamily != familiesEnd; currFamily++){
	    (*currFamily)->doMapping_LeavesToSpecies();
	    waiterMutex->lock();
	    waiter->step();
	    waiterMutex->unlock();
	} 
    
	
}


void DataBase::doFamiliesMapping_NodesToTaxa() {
    // to use the mapping “node on taxon”, we must ensure the mapping “species on leave” has been performed
    doFamiliesMapping_LeavesToSpecies();
    
    if(!mappingDone_NodesToTaxa){
	cout << "Mapping nodes to Taxa:" << endl;
	Waiter patienteur(&cout, nbFamilies, '#');
	for(vector<Family*>::iterator currFamily = families.begin(); currFamily != families.end(); currFamily++){
	    (*currFamily)->doMapping_NodesToTaxa();
	    patienteur.step();
	}
    mappingDone_NodesToTaxa = true;
    patienteur.drawFinal();
    }
	
}

void DataBase::doFamiliesMapping_NodesToUnicityScores() {
    if(!mappingDone_NodesToUnicityScores) {
	Waiter patienteur(&cout, nbFamilies, 'o');
	for(vector<Family*>::iterator currFamily = families.begin(); currFamily != families.end(); currFamily++){
	    (*currFamily)->doMapping_NodesToUnicityScores();
	    patienteur.step();
	}
	mappingDone_NodesToUnicityScores = true;
	patienteur.drawFinal();
    }
}

void DataBase::doFamiliesMapping_NodesToBestUnicityScores() {
    if(!mappingDone_NodesToUnicityScores) {
	Waiter patienteur(&cout, nbFamilies, 'o');
	for(vector<Family*>::iterator currFamily = families.begin(); currFamily != families.end(); currFamily++){
	    (*currFamily)->doMapping_NodesToBestUnicityScores();
	    patienteur.step();
	}
	mappingDone_NodesToUnicityScores = true;
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
	stringstream * curPreambule;
	string *curNewick;
	bool familleSuivante;
	bool crochetFermant;
	Family * laFamille;
	Waiter patienteur(&cout, nbFamilies, '#');
	while(famillesChargees < nbFamilies) {
		// la première ligne est le nom de la famille, on l'ajoute
		// dynamic allocation: will be destroyed during family creation
		curPreambule = new stringstream;
		curNewick = new string;
		familleSuivante = false;
		crochetFermant = false;
		while(!familleSuivante && getline(RAPfile,currLigne)){
			if(crochetFermant) { // on est sur la dernière ligne
				familleSuivante = true;
				*curNewick = currLigne + '\n';
			} else {
				*curPreambule << currLigne << endl;
				if(currLigne.at(0) == ']') crochetFermant = true;
			}
		}
		laFamille = new Family(curPreambule,curNewick,this);
		families.push_back(laFamille);
		famillesChargees++;
		patienteur.step();
		
	}
	patienteur.drawFinal();
	
	cout << "Initializing families:" << endl;
	Waiter patienteur2(&cout, nbFamilies, '#');
	boost::thread_group tg;
	unsigned int blockSize = families.size() / nbThreads;
	cout << "Multithreaded operation. Lot size : " << blockSize << endl;
	for(unsigned int currThreadIndex = 0; currThreadIndex < nbThreads; currThreadIndex++){
	    vector<Family*>::iterator familyBegin, familyEnd;
	    familyBegin = families.begin() + (blockSize*currThreadIndex);
	    if(currThreadIndex+1 != nbThreads) familyEnd = families.begin() + (blockSize*(currThreadIndex+1)); else familyEnd = families.end();
	    boost::thread *currThread = new boost::thread(Family::threadWork_initialize,&patienteur2,&waiterMutex,familyBegin,familyEnd);
	    tg.add_thread(currThread);
	}
	tg.join_all();
	patienteur2.drawFinal();
	
	
}

void DataBase::loadSpeciesTree(string newickLine)
{
	// nettoyage des synonymes ET indexation des espèces de la base
	stringstream sNewickPropre;
	char currChar;
	
	
	bool insideQuote = false;
	bool firstQuoteSeen = false;
	
	for(unsigned int jogger=0; jogger < newickLine.size(); jogger++) {
		currChar = newickLine.at(jogger);
		
		if(currChar == '"')
			if(insideQuote) {insideQuote = false; firstQuoteSeen = true;}
			else insideQuote = true;
				
		else if(firstQuoteSeen && currChar == ':')
		    //reinitialization
		    { firstQuoteSeen = false; sNewickPropre << currChar;}
			
		else if(!firstQuoteSeen) if(!(insideQuote && (currChar == ',' || currChar == ':')))
			sNewickPropre << currChar;
	}
	
	string snp = sNewickPropre.str();
	speciesTree = tpms::TreeTools::newickToTree(snp,true);
	
	
	//Debug : toutes les espèces :
	//cout << "\n\n\n### Voici la liste de toutes les espèces de l'abre des espèces :" << endl;	
	
	vector<Node *> listeNoeuds = speciesTree->getNodes();
		
	for(vector<Node *>::iterator it = listeNoeuds.begin(); it != listeNoeuds.end(); it++) {
		if((*it)->hasName()) {
		    taxa.insert(pair<string,Taxon*>((*it)->getName(),new tpms::Taxon((*it)->getName(),(*it),*this)));
		}
	}
	
	for(map<string,Taxon*>::iterator ct = taxa.begin(); ct != taxa.end(); ct++){
	    ct->second->genRelations();
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

// int DataBase::nbFamiliesContaining(string pTax){
// 	// première possibilité : on a déjà fait la recherche
// 	map<string,int>::iterator nbIt = nbFC.find(pTax);
// 	if(nbIt != nbFC.end()) return(nbIt->second);
// 	
// 	// deuxième cas, il faut faire la recherche et lister les familles
// 	int result = 0;
// 	
// 	// on développe le taxon
// 	set<string> tMembers = getDescendants(pTax); // contient la liste de toutes les taxons appartenant au taxon pTax
// 	
// 	for(vector<Family *>::iterator oneFam = families.begin(); oneFam != families.end(); oneFam++) {
// 		bool auMoinsUnTaxonPourCetteFamille = false;
// 		for(set<string>::iterator unTaxon = tMembers.begin(); !auMoinsUnTaxonPourCetteFamille && unTaxon != tMembers.end(); unTaxon++)
// 			auMoinsUnTaxonPourCetteFamille = auMoinsUnTaxonPourCetteFamille || (*oneFam)->getSpecies()->find(*unTaxon) != (*oneFam)->getSpecies()->end();
// 		if(auMoinsUnTaxonPourCetteFamille) result++;
// 	}
// 	nbFC.insert(pair<string,int>(pTax,result));
// 	return(result);
// }


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

tpms::Taxon* DataBase::nameToTaxon(string taxonName)
{
    boost::to_upper(taxonName);
    map<string,Taxon*>::iterator foundTaxon = taxa.find(taxonName);
    if(foundTaxon != taxa.end()) return(foundTaxon->second);
    cout << "Unable to find the taxon named " << taxonName << " in the database" << endl;
    return(00);
}



}
