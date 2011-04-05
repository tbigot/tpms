#ifndef DATABASE_HPP_TPMS
#define DATABASE_HPP_TPMS


#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <map>

#include "Family.hpp"
#include "Taxon.hpp"

//inclusions bio++
#include <Bpp/Phyl/Tree.h>

class DataBase {
	
	private:
		std::map<std::string,tpms::Taxon*> taxa;
		unsigned int nbFamilies;
		bool reconciled;
		bpp::TreeTemplate<bpp::Node> * speciesTree;
		std::vector<Family *> families;
		std::map<std::string,int> nbFC; // cache de recherche de la fonction nbFamiliesContaining
		std::set<std::string> species;
		
		void loadFamily(std::stringstream * intro, std::string newick);
		void loadSpeciesTree(std::string newickLine);
		void loadFromFile(std::ifstream& RAPfile);
		void renameNodes(bpp::TreeTemplate<bpp::Node> *);
		
		bool speciesTreesBuilded;
		bool refTreesBuilded;
		bool unicityScoresComputed;

				
		
				
	public:
		//constructeur à partir d'un fichier
		DataBase(std::string path);
		bpp::TreeTemplate<bpp::Node> * getSpeciesTree();
		std::vector<Family *> & getFamilies();
		unsigned int getNbFamilies();
		std::string getParentTaxon(std::string pTaxon, unsigned int level);
		bool taxonExists(std::string ptax);
		
		//nombre de familles contenant un taxon donné
		int nbFamiliesContaining(std::string pTax);
		
		std::set<std::string> getAllNodes(bpp::Node * localRoot,bool nodesWanted = true);
		std::set<std::string> getDescendants(std::string taxon, bool nodesWanted = true);
		std::set<std::string> getDescendants(std::vector<std::string> taxaList, bool nodesWanted = true);
		
		
		void iNeedSpeciesTrees(bool verbose, std::string path,bool generate=false);
		void iNeedMapping(bool verbose, std::string path, bool generate=false);
		void genUnicityScores();
		void genBestUnicityScores();
		
		std::set<std::string> * getSpecies();
		tpms::Taxon* nameToTaxon(std::string taxonName);

};
#else

class DataBase;

#endif
