
#ifndef DATABASE_HPP_TPMS
#define DATABASE_HPP_TPMS


#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <map>

#include <boost/thread.hpp>

#include "Family.hpp"
#include "Taxon.hpp"
#include "Waiter.hpp"

//inclusions bio++
#include <Bpp/Phyl/Tree.h>



namespace tpms{
class DataBase {
	
	private:
		std::map<std::string,tpms::Taxon*> taxa;
		unsigned int nbFamilies;
		bool reconciled;
		bpp::TreeTemplate<bpp::Node> * speciesTree;
		std::vector<Family *> families;
		std::map<std::string,int> nbFC; // cache de recherche de la fonction nbFamiliesContaining
		std::set<std::string> species;
		
		// number of threads that can be used to compute data
		unsigned int nbThreads;
		boost::mutex m;
		boost::mutex waiterMutex;
		
		void loadFamily(std::stringstream * intro, std::string newick);
		void loadSpeciesTree(std::string newickLine);
		void loadFromFile(std::ifstream& RAPfile);
		void renameNodes(bpp::TreeTemplate<bpp::Node> *);
		
		bool mappingDone_LeavesToSpecies;
		bool mappingDone_NodesToTaxa;
		bool mappingDone_NodesToUnicityScores;
		

				
		
				
	public:
		//constructeur à partir d'un fichier
		DataBase(std::string path, unsigned int nbThreads = 1);
		bpp::TreeTemplate<bpp::Node> * getSpeciesTree();
		std::vector<Family *> & getFamilies();
		unsigned int getNbFamilies();
		std::string getParentTaxon(std::string pTaxon, unsigned int level);
		
		

		//nombre de familles contenant un taxon donné
		// int nbFamiliesContaining(std::string pTax);
		
		std::set<std::string> getAllNodes(bpp::Node * localRoot,bool nodesWanted = true);
		
		static void doFamiliesMapping_LeavesToSpecies_oneThread(Waiter * waiter, boost::mutex *waiterMutex, std::vector<tpms::Family*>::iterator &familiesBegin, std::vector<tpms::Family*>::iterator &familiesEnd);

		
		void doFamiliesMapping_LeavesToSpecies();
		void doFamiliesMapping_NodesToTaxa();
		void doFamiliesMapping_NodesToUnicityScores();
		void doFamiliesMapping_NodesToBestUnicityScores();
		void doFamiliesMapping_NodesToLowestTaxa();
		void doFamiliesMapping_NodesToTaxonomicShift();
		
		void doFamiliesComputation_detectTransfers(std::ofstream *output);
		
		std::set<std::string> * getSpecies();
		tpms::Taxon* nameToTaxon(std::string taxonName);
		
};}
#else

namespace tpms{
    class DataBase;
}

#endif