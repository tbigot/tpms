#ifndef TPMS_TAXON
#define TPMS_TAXON
#include <string>
#include <Bpp/Phyl/TreeTemplate.h>
#include "DataBase.hpp"


namespace tpms{
    class Taxon{
	private:
	    std::string _name;
	    std::set<tpms::Taxon> _ancestors;
	    std::set<tpms::Taxon> _descendants;
	    
	    void genDescendantsList(bpp::TreeTemplate<bpp::Node*>* speciesTree);
	    void genAncestorsList(bpp::TreeTemplate<bpp::Node*>* speciesTree);
	    
	    
	public:
	    Taxon(std::string name,bpp::Node *nodeInSpTree, DataBase &database);
	    std::string getName();
	    static const Taxon* notFound = 00;
    };
}

#else

namespace tpms{
    class Taxon;
}
#endif
