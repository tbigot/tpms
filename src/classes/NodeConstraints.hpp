//
// File: NodeConstraints.hpp
// Created by: Thomas Bigot
//

/*
   Copyright or © or Copr. Thomas Bigot 2012

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef TPMS_NODECONSTRAINTS
#define TPMS_NODECONSTRAINTS

#include "DataBase.hpp"
#include "Family.hpp"
#include "Taxon.hpp"


#include <vector>
#include <set>
#include <string>
namespace tpms{
    
class NodeConstraints{
    
    private:
	DataBase &refDB;
	NodeNature nature;
	std::set<tpms::Taxon*> allowedSpeciesOnSubtree;
	std::set<tpms::Taxon*> allowedSpeciesOnNode;
	
	NodeType type;
	unsigned int minBootstrap;
	
	std::string constraintsOnNodeString;
	std::string subtreeConstraintsString;
        
        bool ok;
	
	/**
	 * @brief True means the target node must be direct son of the node matching to his pattern father node.
	 **/
	bool direct;
	
	void buildAllowedSpecies(std::set<tpms::Taxon*>& spset,std::string spstr, DataBase &pRefDB);
	
	void addTaxon(std::set<tpms::Taxon*>& spset,std::string);
	void deleteTaxon(std::set<tpms::Taxon*>& spset,std::string);
		
	bool speciesRestrictionsAsSon;
	bool speciesRestrictions;
	
    public:
        bool isOk();
	NodeConstraints(DataBase & pRefDB);
	void setConstraints(DataBase &pRefDB, std::string constraintString, NodeType type);
	bool allowsAsSon(Family& family, bpp::Node* node);
	bool allows(Family& family, bpp::Node * node);
	std::set<tpms::Taxon*>& getAllowedSpecies();
	bool isLeaf();
	std::string getStr();
        bool isDirect();
	
	
    
};
    
}

#else
namespace tpms{

class NodeConstraints;
}
#endif
