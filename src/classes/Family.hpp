#ifndef FAMILY_HPP_TPMS
#define FAMILY_HPP_TPMS

#include <fstream>
#include <string>
#include <sstream>
#include <set>

//inclusions bio++
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>

//inclusions personnelles
#include "DataBase.hpp"

class Family {
    
    private:
	std::string _name;
	//! Base de données d'appartenance
	DataBase * db;
	
	bpp::TreeTemplate<bpp::Node> * tree;
	bpp::TreeTemplate<bpp::Node> * spTree;
	std::map<std::string,std::string> mne2spec;
	std::set<std::string> species;
	bpp::TreeTemplate<bpp::Node> * refTree;
	std::vector<unsigned int> unicityScores;
	
	void renameNodes(bpp::TreeTemplate<bpp::Node> *);
	
	
	void deleteFromLeavesToBif(bpp::Node * pnode);
	
	void computeUnicity(std::map<std::string, unsigned int> &thisNodeCount, bpp::Node * node);
	
	//! Suppression des fils uniques
	/*!
	Raccourci les embranchements en supprimant les nœuds n'ayant qu'un fils, inutiles du point de vue topologique.
	*/
	bpp::Node * removeUniqueSons(bpp::Node* localRoot);
	
	std::string mapNodeOnTaxon(bpp::Node& nodep);
	
	void writeRefTreeToFile(std::string path);
	void writeSpTreeToFile(std::string path);
	

    public:
	//constructeur à partir d'un fichier
	Family(std::stringstream* sIntro, std::string sNewick, DataBase* dbp);
	
	bpp::TreeTemplate<bpp::Node> * getTree();
	bpp::TreeTemplate<bpp::Node> * getSpTree();
	void genUnicityScores();
	
	bool containsSpecie(std::string specie);
	bool containsSpecies(std::set<std::string> speciesList);
	std::set<std::string> * getSpecies();
	static void getLeavesFromNode(bpp::Node * pnode, std::vector< bpp::Node* >& leaves);
	static void getLeavesFromNode(bpp::Node * pnode, std::set< bpp::Node* >& leaves, int &leavesNumber);
	static void getLeavesFromNode(bpp::Node * pnode, std::set< bpp::Node* >& leaves);
	
	std::set<std::string> getLeavesNamesFromNode(bpp::Node * pnode);
	static std::set<std::string> nodesToNames(std::set<bpp::Node *> &nodes);
	std::string getName();
	void writeTreeToStream(bpp::Node* root, std::ostream& sortie, unsigned int deep);
	
	//! Génération du Vrai arbre
	/*!
	Ce générateur va partir d'un arbre vrai de référence et n'en garder que les espèces présentes dans la famille
	*/
	
	void genRefTree(bool save=true, std::string path="");
	void genSpTree(bool save=true, std::string path="");
	std::vector<unsigned int> &getUnicityScores();
	
	int numberOfNodesBetween(bpp::Node * ancestor, bpp::Node * pnode);
	
	void loadRefTreeFromFile(std::string path);
	void loadSpTreeFromFile(std::string path);
	
	void addSequencesNames(bpp::Node* currNode);
	
	/**
	    * @brief Replace names of the pieces of pattern trees nodes by the name of the sequences
	    * 
	    * @param currNode node to start with
	    */
	void labelWithSequencesNames(bpp::Node* currNode);
	
	
};

#else

class Family;

#endif
