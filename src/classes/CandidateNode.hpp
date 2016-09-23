//
// File: CandidateNode.hpp
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

#ifndef TPMS_CANDIDATENODE
#define TPMS_CANDIDATENODE

#include <set>
#include <map>

#include <Bpp/Phyl/Tree/Node.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>

namespace tpms{
    class CandidateNode{
	private:
	    bpp::Node * treeNode_;
	    bpp::Node * patternNode_;
	    tpms::CandidateNode * father_;
	    double distanceToFather_;
	    std::map<bpp::Node *, std::vector<tpms::CandidateNode *> > sons_;
	    bool isRoot_;
	    unsigned int addsons_;
	    bool isLeaf_();
	    
	    /**
	    * @brief Recursively builds result trees from CandidateNode trees
	    * 
	    * @return the nodes of trees that can be built from this CandidateNode
	    */
	    std::vector<bpp::Node *> recGenTrees_();
	    
	public:
	    CandidateNode(tpms::CandidateNode * father, bpp::Node * treeNode, bpp::Node * patternNode);
	    CandidateNode();
	    ~CandidateNode();
	    
	    bpp::Node * getTreeNode();
	    
	    void addSon(tpms::CandidateNode * son, bpp::Node * patternNode);
	    void confirm();
	    // cette fonction ne fonctionne que sur des Candidates qui sont racines
	    unsigned int genTrees(std::vector<bpp::TreeTemplate<bpp::Node> *> &trees);
	    
	    
	    /**
	    * @brief Get whole matching subtrees.
	    * 
	    * @param trees a reference of a tree vector that will be filled
	    * @return the number of generated trees
	    */
	    unsigned int getWholeMatchingSubtrees(std::vector<bpp::TreeTemplate<bpp::Node> *> &trees);

	    void print(unsigned int lvl);
    };
}

#else

namespace tpms{
    class CandidateNode;
}
#endif
