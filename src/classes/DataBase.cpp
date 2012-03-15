#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>


#define BOOST_FILESYSTEM_VERSION 2
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

DataBase::DataBase(string path, unsigned int nbThreads): mappingDone_NodesToTaxa(false), mappingDone_LeavesToSpecies(false), mappingDone_NodesToUnicityScores(false), nbThreads(nbThreads) {
	
	// checking this file exists
	fs::path dbFile(path);
	if(!fs::exists(dbFile)) {
	    cout << "The file "<< path << " does not exists, exiting."<< endl;
	    exit(1);
	}
	
	// using filesystem 3:
	// filename = dbFile.filename().string();
	
	// but for compatibility, we use filesystem 2:
	filename = dbFile.filename();
	
    
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

std::string DataBase::getStatus(std::string preText){
    ostringstream out;
    out << preText << "Database Status" << endl;
    out << preText << "  name\t" << filename<<endl;
    out << preText << "  reconciled\t" << (reconciled?"yes":"no")<< endl;
    out << preText << "  unicity scores computed\t" << (mappingDone_NodesToUnicityScores?"yes":"no")<< endl;
    return(out.str());
}

void DataBase::doFamiliesMapping_LeavesToSpecies() {
    // multithreading version
    if(!mappingDone_LeavesToSpecies){
	cout << "Mapping leaves to species:" << endl;
	Family::threadedWork_launchJobs(families,&Family::doMapping_LeavesToSpecies,nbThreads);
    }
    mappingDone_LeavesToSpecies=true;
}


void DataBase::doFamiliesMapping_NodesToTaxa() {
    // to use the mapping “node on taxon”, we must ensure the mapping “species on leave” has been performed
    doFamiliesMapping_LeavesToSpecies();
    if(!mappingDone_NodesToTaxa){
	cout << "Mapping nodes to Taxa:" << endl;
	Family::threadedWork_launchJobs(families,&Family::doMapping_NodesToTaxa,nbThreads);
	mappingDone_NodesToTaxa = true;
    }

}

void DataBase::doFamiliesRerooting_Taxonomy() {
    // to use the mapping “node on taxon”, we must ensure the mapping “species on leave” has been performed
    doFamiliesMapping_LeavesToSpecies();
    if(!mappingDone_NodesToTaxa){
	cout << "Re-rooting family trees using the taxonomic criteria:" << endl;
	Family::threadedWork_launchJobs(families,&Family::doRerooting_Taxonomy,nbThreads);
	mappingDone_NodesToTaxa = true;
    }

}

void DataBase::doFamiliesRerooting_UnicityTaxonomy() {
    // to use the mapping “node on taxon”, we must ensure the mapping “species on leave” has been performed
    doFamiliesMapping_LeavesToSpecies();
    if(!mappingDone_NodesToTaxa){
        cout << "Re-rooting families trees using a combo method Unicity+Taxonomy:" << endl;
        Family::threadedWork_launchJobs(families,&Family::doRerooting_UnicityTaxonomy,nbThreads);
        mappingDone_NodesToTaxa = true;
    }

}

void DataBase::doFamiliesRerooting_Unicity() {
    doFamiliesMapping_LeavesToSpecies();
    cout << "Re-rooting families trees with the unicity criteria:" << endl;
    Family::threadedWork_launchJobs(families,&Family::doRerooting_Unicity,nbThreads);
    mappingDone_NodesToUnicityScores = true;
}


void DataBase::doFamiliesMapping_NodesToUnicityScores() {
    doFamiliesMapping_LeavesToSpecies();
    cout << "UnicityScores computing:" << endl;
    Family::threadedWork_launchJobs(families,&Family::doMapping_NodesToUnicityScores,nbThreads);
    mappingDone_NodesToUnicityScores = true;
}


void DataBase::doFamiliesMapping_NodesToTaxonomicShift(){
    doFamiliesMapping_LeavesToSpecies();
    cout << "Nodes to taxonomic shift (grandfather affectation change):" << endl;
    Family::threadedWork_launchJobs(families,&Family::doMapping_NodesToTaxonomicShift,nbThreads);
}

void DataBase::doFamiliesComputation_detectTransfers(ostream * output){
    doFamiliesMapping_LeavesToSpecies();
    cout << "Transfers Detection:" << endl;
    Family::threadedWork_launchJobs(families,&Family::compute_detectTransfers,nbThreads,output);
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
	cout << "\n --  Reading the families from the collection file:" << endl;
	
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
	
	// initializing families
	
        cout << "\n --  Building families from the read data:" << endl ;
	Family::threadedWork_launchJobs(families,&Family::initialize, nbThreads);

	
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
	    string taxonToCreateName;
	    if ((*it)->hasName()) taxonToCreateName = (*it)->getName();
	    Taxon* currTaxon = new tpms::Taxon(taxonToCreateName,(*it),*this);
	    nodeToTaxon.insert(  pair<unsigned int, Taxon*> ( (*it)->getId(),currTaxon )  );
	    if(!taxonToCreateName.empty()) {
		taxa.insert(pair<string,Taxon*>(taxonToCreateName,currTaxon));
	    }
	}
	
	for(map<unsigned int,Taxon*>::iterator ct = nodeToTaxon.begin(); ct != nodeToTaxon.end(); ct++){
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

Taxon* DataBase::nameToTaxon(string taxonName)
{
    boost::to_upper(taxonName);
    map<string,Taxon*>::iterator foundTaxon = taxa.find(taxonName);
    if(foundTaxon != taxa.end()) return(foundTaxon->second);
    cout << "Unable to find the taxon named " << taxonName << " in the database" << endl;
    return(00);
}

Taxon* DataBase::nodeIdToTaxon(unsigned int nodeId){
    map<unsigned int,Taxon*>::iterator foundTaxon = nodeToTaxon.find(nodeId);
    if(foundTaxon != nodeToTaxon.end()) return(foundTaxon->second);
    cout << "Unable to find node "  << nodeId << endl;
    return(00);
}


}
