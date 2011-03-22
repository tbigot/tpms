#ifndef NODECONSTRAINTS
#define NODECONSTRAINTS

#include "DataBase.hpp"

#include <vector>
#include <set>
#include <string>

class NodeConstraints{
    
    public:
	enum NodeNature {ANY, DUPLICATION, SPECIATION};
	
    
    private:
	DataBase &refDB;
	NodeNature nature;
	std::set<std::string> authorisedSpecies;
	void buildAuthorisedSpecies(std::string speciesList);
	void addTaxon(std::string);
	void deleteTaxon(std::string);
	std::string initString;
	void cachParams();
	bool speciesRestrictions;
	
    public:
	NodeConstraints(DataBase & pRefDB);
	NodeConstraints(DataBase & pRefDB,std::string constraintString);
	void setConstraints(DataBase &pRefDB, std::string constraintString);
	void setNature(NodeNature);
	void setAuthorisedSpecies(std::set<std::string> &authorisedSpecies);
	void setAuthorisedSpecies(std::string authSpString);
	std::set<std::string> & getAuthorisedSpecies();
	bool hasSpeciesRestrictions();
	std::string getString();
	NodeNature getNature();
	bool isAuthorizedNature(NodeNature nature);
	
	
    
};

#else

class NodeConstraints;

#endif
