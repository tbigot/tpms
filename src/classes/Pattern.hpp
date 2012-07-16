//
// File: Pattern.hpp
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

#ifndef TPMS_PATTERN
#define TPMS_PATTERN

#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <vector>

#include "Family.hpp"
#include "NodeConstraints.hpp"
#include "DataBase.hpp"
#include "CandidateNode.hpp"
#include "Taxon.hpp"

//inclusions bio++
#include <Bpp/Phyl/Tree.h>
namespace tpms{

class Pattern {
    
private:
    
    /**
     * @brief Are all the node constraints Ok?
     **/
    bool ok;
    
    /**
     * @brief The database associated to this pattern.
     **/
    tpms::DataBase &db;
    
    
     /**
     * @brief The pattern tree itself.
     **/
    bpp::TreeTemplate<bpp::Node> &tree;
    
    /**
     * @brief All species that can be in this pattern.
     **/
    std::set<tpms::Taxon*> allSpecies;
    

    /**
     * @brief Transform the first part of a node name into a constraint node
     *
     **/
    void extractConstraints();
    
    /** @} */
    

    
 /**
     * @brief Pattern Maching recursive function
     *
     **/
    // fonction récursive de recherche de motifs (Pattern Matching)
    bool patternMatch(Family& family,bpp::Node* target, bpp::Node* pattern, tpms::CandidateNode* fatherCandidate);
    
    
        
    /**
     * @name Small utilities
     * @{
     *
     **/
    static bool isLeaf(bpp::Node * pNode);
    
    
    /** @} */
    
    std::vector<int> getIdWithTaxaList(bpp::TreeTemplate<bpp::Node> * sTree, std::set<std::string> * taxa);
    int enumerateTaxon(bpp::Node * localRoot, std::set<std::string> * tsMembers, bool tsComplementary, bpp::Node * impasse, unsigned int * totsize);
    
    std::vector<NodeConstraints *> constraints;
    void toString(bpp::Node * noeud, int spc, std::ostream &outputStream);
    
    NodeConstraints* constraintsOf(bpp::Node* node);
    bool isTreeBinary();
    
    std::vector<unsigned int> mapping_NodesToMaxDepths;
    
public:
    
    unsigned int mapNodeToMaxDepth(bpp::Node*);
    
    bool patternMatchInit(Family& family, tpms::CandidateNode* initCnode);

    
    bool isOk();
        
    Pattern(bpp::TreeTemplate< bpp::Node >& tree, DataBase& pRefDB);
    ~Pattern();
    
    void toString(std::ostream&);
    
    static bool nodeOnlyContainsTheseTaxa(bpp::Node * localRoot, std::set<tpms::Taxon*> & taxonMembers, bool invert=false);
    
    static bool isUnder(bpp::Node * descendant, bpp::Node * ascendant);
    unsigned int search(std::vector< Family* >& families, std::vector< std::pair< Family*, tpms::CandidateNode* > >& result);
    std::vector<int> xferDetected(std::map<int,Family *> * families, std::vector<int> famList, std::string sourceTaxon, std::string targetTaxon, std::string monophylyTaxon, unsigned int verifDeep, std::vector<unsigned int> bootstraps,  std::vector<unsigned int> sourceRates);
    std::vector<int> xferGapDetect(std::map<int,Family *> families, std::set<std::string>);
};
}

#else
namespace tpms{

class Pattern;
}
#endif