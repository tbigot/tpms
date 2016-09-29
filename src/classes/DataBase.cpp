//
// File: DataBase.cpp
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


#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>


//#define BOOST_FILESYSTEM_VERSION 2
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string.hpp>


#include "TreeTools.hpp"
#include "DataBase.hpp"
#include "Waiter.hpp"
#include "Family.hpp"
#include "Taxon.hpp"
#include "Pattern.hpp"


#include <Bpp/Phyl/Tree/Tree.h>


using namespace std;
using namespace bpp;
using namespace tpms;

namespace fs = boost::filesystem;


namespace tpms{

DataBase::DataBase(string path, bool expectSynonyms, unsigned int nbThreads):
mappingDone_LeavesToSpecies_(false),
mappingDone_NodesToTaxa_(false),
mappingDone_NodesToUnicityScores_(false),
mappingDone_NodesToMaxDepth_(false),
expectSynonyms_(expectSynonyms),
nbThreads_(nbThreads)
{
	
	// checking this file exists
	fs::path dbFile(path);
	if(!fs::exists(dbFile)) {
	    cout << "The file "<< path << " does not exist, exiting."<< endl;
	    exit(1);
	}

	filename_ = dbFile.filename().string();

    // ouverture du fichier
	ifstream RAPfile(path.c_str(), ifstream::in);
	
	// chargement à partir du fichier
	loadFromFile_(RAPfile);
	
	// on ferme le fichier proprement :
	RAPfile.close();
	
	//DEBUG OUT
	cout << "\nLoad sucessful!\n" << endl;
	
	std::cout << "Species tree loaded: " << speciesTree_->getNumberOfLeaves() << " species." << std::endl;
	
	cout << "Loaded " << families_.size() << " families / " << nbFamilies_ << " expected." << endl;
	if(families_.size() != nbFamilies_) {
		cout << "WARNING : The file " << path << " is supposed to contain more trees." << endl;
	}
	if(reconciled_) cout << "Trees are reconciled."; else cout << "Trees are not reconciled.";
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

bpp::TreeTemplate<bpp::Node> * DataBase::getSpeciesTree() { return(speciesTree_); }

std::string DataBase::getStatus(std::string preText){
    ostringstream out;
    out << preText << "Database Status" << endl;
    out << preText << "  name\t" << filename_<<endl;
    out << preText << "  reconciled\t" << (reconciled_?"yes":"no")<< endl;
    out << preText << "  unicity scores computed\t" << (mappingDone_NodesToUnicityScores_?"yes":"no")<< endl;
    return(out.str());
}

void DataBase::doFamiliesMapping_NodesToMaxDepth() {
    // multithreading version
    if(!mappingDone_NodesToMaxDepth_){
        cout << "Mapping nodes to max depth:" << endl;
        Family::threadedWork_launchJobs(families_,&Family::doMapping_NodesToMaxDepth,nbThreads_);
    }
    mappingDone_NodesToMaxDepth_=true;
}


void DataBase::doFamiliesMapping_LeavesToSpecies() {
    // multithreading version
    if(!mappingDone_LeavesToSpecies_){
	cout << "\n\nMapping leaves to species:" << endl;
	Family::threadedWork_launchJobs(families_,&Family::doMapping_LeavesToSpecies,nbThreads_);
    }
    mappingDone_LeavesToSpecies_=true;
}


void DataBase::doFamiliesMapping_NodesToTaxa() {
    // to use the mapping “node on taxon”, we must ensure the mapping “species on leave” has been performed
    doFamiliesMapping_LeavesToSpecies();
    if(!mappingDone_NodesToTaxa_){
	cout << "\n\nMapping nodes to Taxa:" << endl;
	Family::threadedWork_launchJobs(families_,&Family::doMapping_NodesToTaxa,nbThreads_);
	mappingDone_NodesToTaxa_ = true;
    }

}

void DataBase::doFamiliesRerooting_Taxonomy() {
    // to use the mapping “node on taxon”, we must ensure the mapping “species on leave” has been performed
    doFamiliesMapping_LeavesToSpecies();
    if(!mappingDone_NodesToTaxa_){
	cout << "Re-rooting family trees using the taxonomic criteria:" << endl;
	Family::threadedWork_launchJobs(families_,&Family::doRerooting_Taxonomy,nbThreads_);
	mappingDone_NodesToTaxa_ = true;
    }

}

void DataBase::doFamiliesRerooting_UnicityTaxonomy() {
    // to use the mapping “node on taxon”, we must ensure the mapping “species on leave” has been performed
    doFamiliesMapping_LeavesToSpecies();
    if(!mappingDone_NodesToTaxa_){
        cout << "Re-rooting families trees using a combo method Unicity+Taxonomy:" << endl;
        Family::threadedWork_launchJobs(families_,&Family::doRerooting_UnicityTaxonomy,nbThreads_);
        mappingDone_NodesToTaxa_ = true;
    }

}

void DataBase::doFamiliesRerooting_Daubin() {
    // to use the mapping “node on taxon”, we must ensure the mapping “species on leave” has been performed
    doFamiliesMapping_LeavesToSpecies();
    cout << "Re-rooting families trees using Daubin’s criteria:" << endl;
    Family::threadedWork_launchJobs(families_,&Family::doRerooting_Daubin,nbThreads_);

}

void DataBase::doFamiliesRerooting_Unicity() {
    doFamiliesMapping_LeavesToSpecies();
    cout << "Re-rooting families trees with the unicity criteria:" << endl;
    Family::threadedWork_launchJobs(families_,&Family::doRerooting_Unicity,nbThreads_);
    mappingDone_NodesToUnicityScores_ = true;
}


void DataBase::doFamiliesRerooting_LessTransfers(ostream * output) {
    doFamiliesMapping_LeavesToSpecies();
    cout << "Re-rooting families trying to minimize the number of transfers:" << endl;
    Family::threadedWork_launchJobs(families_,&Family::doRerooting_LessTransfers,nbThreads_,output);
}


void DataBase::doFamiliesMapping_NodesToUnicityScores() {
    doFamiliesMapping_LeavesToSpecies();
    cout << "UnicityScores computing:" << endl;
    Family::threadedWork_launchJobs(families_,&Family::doMapping_NodesToUnicityScores,nbThreads_);
    mappingDone_NodesToUnicityScores_ = true;
}


void DataBase::doFamiliesMapping_NodesToTaxonomicShift(){
    doFamiliesMapping_LeavesToSpecies();
    cout << "Nodes to taxonomic shift (grandfather affectation change):" << endl;
    Family::threadedWork_launchJobs(families_,&Family::doMapping_NodesToTaxonomicShift,nbThreads_);
}

void DataBase::doFamiliesComputation_detectTransfers(ostream * output){
    doFamiliesMapping_LeavesToSpecies();
    cout << "Transfers Detection:" << endl;
    Family::threadedWork_launchJobs(families_,&Family::compute_detectTransfers,nbThreads_,output);
}

void DataBase::loadFromFile_(ifstream & RAPfile) {
	// on extrait les lignes
	
	cout << "Loading preamble..." << flush;
	
	string currLigne;
	// ===PREMIERE LIGNE=== nous donne le nombre d'arbres contenus et leur réconciliation
	getline(RAPfile, currLigne);
	
	// Nombre d'arbres
	int currEspace = currLigne.find(" ");
	stringstream ssNbFam;
	ssNbFam << currLigne.substr(0,currEspace);
	ssNbFam >> nbFamilies_;
	
	// Réconciliés ?
	int espacePrecedent = currEspace;
	currEspace = currLigne.find(" ",espacePrecedent+1);
	string sReconciled = currLigne.substr(espacePrecedent+1,currEspace-espacePrecedent-1);
	
	if(sReconciled=="unreconciled") reconciled_ = false; else reconciled_ = true;
	
	cout << "    [DONE]" << endl;
	// la seconde ligne contient l'arbre exhaustif des espèces
	cout << "Loading species tree..." << flush;
	
	getline(RAPfile, currLigne);
	loadSpeciesTree_(currLigne);
	
	cout << "    [DONE]" << endl;

	// maintenant, prise en charge des familles : blocs suivants :
	cout << "\n --  Reading the families from the collection file:" << endl;
	
	unsigned int famillesChargees = 0;
	stringstream * curPreambule;
	string *curNewick;
	bool familleSuivante;
	bool crochetFermant;
	Family * laFamille;
	Waiter patienteur(&cout, nbFamilies_, '#');
	while(famillesChargees < nbFamilies_) {
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
		families_.push_back(laFamille);
		famillesChargees++;
		patienteur.doStep();
		
	}
	patienteur.drawFinal();
	
	// initializing families
	
        cout << "\n --  Building families from the read data:" << endl ;
	Family::threadedWork_launchJobs(families_,&Family::initialize, nbThreads_);

	
}

void DataBase::loadSpeciesTree_(string newickLine)
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
	speciesTree_ = tpms::TreeTools::newickToTree(snp,expectSynonyms_);
	
    cout << "NODES = " << speciesTree_->getNodes().size() << " ; LEAVES = " << speciesTree_->getLeaves().size() << endl;
	
	//Debug : toutes les espèces :
	//cout << "\n\n\n### Voici la liste de toutes les espèces de l'abre des espèces :" << endl;	
	
	vector<Node *> listeNoeuds = speciesTree_->getNodes();
		
	for(vector<Node *>::iterator it = listeNoeuds.begin(); it != listeNoeuds.end(); it++) {
	    string taxonToCreateName;
	    if ((*it)->hasName()) taxonToCreateName = (*it)->getName();
        boost::to_upper(taxonToCreateName);
	    Taxon* currTaxon = new tpms::Taxon(taxonToCreateName,(*it),*this);
	    nodeToTaxon_.insert(  pair<unsigned int, Taxon*> ( (*it)->getId(),currTaxon )  );
	    if(!taxonToCreateName.empty()) {
		taxa_.insert(pair<string,Taxon*>(taxonToCreateName,currTaxon));
	    }
	}
	
	for(map<unsigned int,Taxon*>::iterator ct = nodeToTaxon_.begin(); ct != nodeToTaxon_.end(); ct++){
	    ct->second->genRelations();
	}
}

vector<Family *> & DataBase::getFamilies() {return(families_);}

unsigned int DataBase::getNbFamilies() {
	return(families_.size());
}

string DataBase::getParentTaxon(string pTaxon, unsigned int level) {
	Node * pNode = speciesTree_->getNode(pTaxon);
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
    return(&species_);
}

Taxon* DataBase::nameToTaxon(string taxonName)
{
    boost::to_upper(taxonName);
    map<string,Taxon*>::iterator foundTaxon = taxa_.find(taxonName);
    if(foundTaxon != taxa_.end()) return(foundTaxon->second);
    cout << "Unable to find the taxon named " << taxonName << " in the database" << endl;
    return(00);
}

Taxon* DataBase::nodeIdToTaxon(unsigned int nodeId){
    map<unsigned int,Taxon*>::iterator foundTaxon = nodeToTaxon_.find(nodeId);
    if(foundTaxon != nodeToTaxon_.end()) return(foundTaxon->second);
    cout << "Unable to find node "  << nodeId << endl;
    return(00);
}


}
