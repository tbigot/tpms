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
	    
	    void genDescendantsList(bpp::Node* localNode, DataBase& database);
	    void genAncestorsList(bpp::Node* localNode, DataBase& database);
	    
	    
	public:
	    Taxon(std::string name,bpp::Node *nodeInSpTree, DataBase &database);
	    std::string getName();
	    
	    bool contains(Taxon *);
	    bool belongsTo(Taxon *);
    };
}

#else

namespace tpms{
    class Taxon;
}
#endif
