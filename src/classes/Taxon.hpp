//
// File: Taxon.hpp
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


#ifndef TPMS_TAXON
#define TPMS_TAXON

#include <string>
#include <Bpp/Phyl/TreeTemplate.h>


#include "DataBase.hpp"
namespace tpms{
    class Taxon{
	private:
	    std::string name;
	    unsigned int depth;
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
            std::string getTaxonomy();
	    
	    bool contains(Taxon *);
	    bool belongsTo(Taxon *);
	    
	    std::set<tpms::Taxon*>& getDescendants();
	    std::set<tpms::Taxon*>& getAncestors();
	    tpms::Taxon* getDirectAncestor();
	    bool hasDirectAncestor();
	    bool hasAncestor();
	    bool containsAllTheseSpecies(std::set<tpms::Taxon*> species);
	    
	    void genRelations();
	    static tpms::Taxon* findLCA(std::set<tpms::Taxon*> taxa);
	    static unsigned int computeRelativeDepthDifference(tpms::Taxon* ancestor, tpms::Taxon* descendant, std::set<tpms::Taxon*>* taxaList);
        unsigned int getRelativeDepth(std::set<Taxon*> &taxaList);
	    
	    // computes distance between this taxon and "taxon"
	    // positive distance = "taxon" is an ancestor of this taxon
	    unsigned int getDepth();
	    
    };
}
#else

namespace tpms{

    class Taxon;
}
#endif
