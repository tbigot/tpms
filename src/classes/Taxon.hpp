
#ifndef TPMS_TAXON
#define TPMS_TAXON

#include <string>
#include <Bpp/Phyl/TreeTemplate.h>


#include "DataBase.hpp"
namespace tpms{
    class Taxon{
	private:
	    std::string name;
	    std::set<tpms::Taxon*> ancestors;
	    std::set<tpms::Taxon*> descendants;
	    DataBase &db;
	    bpp::Node *nodeInSpTree;
	    Taxon* directAncestor;
	    
	    void genDescendantsList(bpp::Node* localNode);
	    void genAncestorsList(bpp::Node* localNode);
	    
	    
	public:
	    Taxon(std::string name,bpp::Node *nodeInSpTree, DataBase &database);
	    std::string getName();
	    
	    bool contains(Taxon *);
	    bool belongsTo(Taxon *);
	    
	    std::set<tpms::Taxon*>& getDescendants();
	    std::set<tpms::Taxon*>& getAncestors();
	    tpms::Taxon* getDirectAncestor();
	    bool hasAncestor();
	    bool containsAllTheseSpecies(std::set<tpms::Taxon*> species);
	    
	    void genRelations();
	    static tpms::Taxon* findSmallestCommonTaxon(std::set<tpms::Taxon*> taxa);
    };
}
#else

namespace tpms{

    class Taxon;
}
#endif
