#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <vector>

#include "NodeConstraints.hpp"
#include "Family.hpp"
#include "DataBase.hpp"
#include "CandidateNode.hpp"

//inclusions bio++
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>

class Pattern {
	
	private:
		int id;
		bpp::TreeTemplate<bpp::Node> & tree;
		std::set<std::string> allSpecies;
		std::map<int, std::set<std::string> > species;
		DataBase & refDB;
		// associe à chaque taxon porté par un nœud de Pattern, une liste ordonnée d'espèces autorisée
		void developSpecies();
		// transforme un nom de nœud qui contient une contrainte en un nom de taxon + un objet contrainte
		void extractConstraints();
		// fonction récursive de recherche de motifs (Pattern Matching)
		bool patternMatch(bpp::Node* target, bpp::Node* pattern, tpms::CandidateNode* fatherCandidate);
		
		static bool isLeaf(bpp::Node * pNode);
		bool isInTaxa(bpp::Node * pNode,std::string pSpecie);
		
		std::vector<int> getIdWithTaxaList(bpp::TreeTemplate<bpp::Node> * sTree, std::set<std::string> * taxa);
		int enumerateTaxon(bpp::Node * localRoot, std::set<std::string> * tsMembers, bool tsComplementary, bpp::Node * impasse, unsigned int * totsize);
		
		std::vector<NodeConstraints *> constraints;
		void print(bpp::Node * noeud, int spc);
		
		bool isTreeBinary();
		
	public:
		//constructeur à partir d'un fichier
		Pattern(bpp::TreeTemplate< bpp::Node >& tree, DataBase& pRefDB);
		~Pattern();
		static bool nodeOnlyContainsTheseTaxa(bpp::Node * localRoot, std::set<std::string> & taxonMembers, bool invert=false);
		static bool isUnder(bpp::Node * descendant, bpp::Node * ascendant);
		unsigned int search(std::vector< Family* >& families, std::vector< std::pair< Family*, tpms::CandidateNode* > >& result);
		std::vector<int> xferDetected(std::map<int,Family *> * families, std::vector<int> famList, std::string sourceTaxon, std::string targetTaxon, std::string monophylyTaxon, unsigned int verifDeep, std::vector<unsigned int> bootstraps,  std::vector<unsigned int> sourceRates);
		std::vector<int> xferGapDetect(std::map<int,Family *> families, std::set<std::string>);
};
