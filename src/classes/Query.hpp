#ifndef QUERY_HPP_TPMS
#define QUERY_HPP_TPMS

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <set>

#include "DataBase.hpp"
namespace tpms{

class Query {
	
	private:
		std::string sourceTaxon;
		std::string targetTaxon;
		unsigned int verifDepth;
		unsigned int monophilyLevel;
		
		std::vector<unsigned int> bootstraps;
		std::vector<unsigned int> sourceRates;
	
	public:
		Query(std::string queryLine,std::string defaultsParam);
		
		std::string getPatternString();
		std::string getSource();
		std::string getTarget();
		unsigned int getDepth();
		unsigned int getMonophilyLevel();
		
		std::vector<unsigned int> getBootstraps();
		std::vector<unsigned int> getSourceRates();
		
    
};
}
#else
namespace tpms{

class Family;

}
#endif
