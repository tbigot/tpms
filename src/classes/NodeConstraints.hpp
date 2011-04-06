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
	std::set<tpms::Taxon*> fromFatherAllowedSpecies;
	std::set<tpms::Taxon*> allowedSpeciesOnNode;
	
	void buildAllowedSpecies(std::set<tpms::Taxon*>& spset,std::string spstr);
	
	void addTaxon(set<Taxon*>& spset,tpms::Taxon* taxon);
	void deleteTaxon(set<Taxon*>& spset,tpms::Taxon* taxon);
	
	std::string initString;
	void catchParams();
	bool speciesRestrictionsFromFather;
	
    public:
	NodeConstraints(DataBase & pRefDB);
	NodeConstraints(DataBase & pRefDB,std::string constraintString);
	void setConstraints(DataBase &pRefDB, std::string constraintString);
	void setNature(NodeNature);
	void setAuthorisedSpecies(std::set<std::string> &authorisedSpecies);
	void setAuthorisedSpecies(std::string authSpString);
	std::set<std::string> & getAuthorisedSpecies();
	bool hasSpeciesRestrictionsFromFather();
	std::string getString();
	NodeNature getNature();
	bool isAuthorizedNature(NodeNature nature);
	
	
    
};

#else

class NodeConstraints;

#endif
