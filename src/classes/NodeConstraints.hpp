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
	
	void addTaxon(std::set<Taxon*>& spset,tpms::Taxon* taxon);
	void deleteTaxon(std::set<Taxon*>& spset,tpms::Taxon* taxon);
	
	std::string initString;
	void catchParams();
	
	bool speciesRestrictionsFromFather;
	bool speciesRestructions;
	
    public:
	NodeConstraints(DataBase & pRefDB);
	NodeConstraints(DataBase & pRefDB,std::string constraintString);
	void setConstraints(DataBase &pRefDB, std::string constraintString);
	std::set<std::string> & getAuthorisedSpecies();
	bool hasSpeciesRestrictionsFromFather();
	std::string getString();
	NodeNature getNature();
	bool isAuthorizedNature(NodeNature nature);
	bool allowsAsSon(Family& family, bpp::Node* node);
	bool allows(Family& family, bpp::Node * node);
	
	
    
};

#else

class NodeConstraints;

#endif
