#ifndef TPMS_TAXON
#define TPMS_TAXON
#include <string>


namespace tpms{
    class Taxon{
	private:
	    std::string _name;
	    
	    
	public:
	    Taxon(std::string name);
	    std::string getName();
	    
    };
}

#else

namespace tpms{
    class Taxon;
}
#endif
