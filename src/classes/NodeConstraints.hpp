#ifndef TPMS_NODECONSTRAINTS
#define TPMS_NODECONSTRAINTS

#include "DataBase.hpp"
#include "Family.hpp"
#include "Taxon.hpp"


#include <vector>
#include <set>
#include <string>

class NodeConstraints{
    
    public:
	enum NodeNature {ANY, DUPLICATION, SPECIATION};
	enum Type {NODE,LEAF};
    
    private:
	DataBase &refDB;
	NodeNature nature;
	std::set<tpms::Taxon*> allowedSpeciesOnSubtree;
	std::set<tpms::Taxon*> allowedSpeciesOnNode;
	
	Type type;
	
	std::string constraintsOnNodeString;
	std::string subtreeConstraintsString;
	
	/**
	 * @brief True means the target node must be direct son of the node matching to his pattern father node.
	 **/
	bool direct;
	
	void buildAllowedSpecies(std::set<tpms::Taxon*>& spset,std::string spstr);
	
	void addTaxon(std::set<tpms::Taxon*>& spset,std::string);
	void deleteTaxon(std::set<tpms::Taxon*>& spset,std::string);
	
	std::string initString;
	void readConstraintsString();
	
	bool speciesRestrictionsAsSon;
	bool speciesRestrictions;
	
    public:
	NodeConstraints(DataBase & pRefDB);
	void setConstraints(DataBase &pRefDB, std::string constraintString, Type type);
	bool allowsAsSon(Family& family, bpp::Node* node);
	bool allows(Family& family, bpp::Node * node);
	std::set<tpms::Taxon*>& getAllowedSpecies();
	bool isLeaf();
	std::string getStr();
	
	
    
};

#else

class NodeConstraints;

#endif
