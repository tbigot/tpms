//
// File: TreeTools.hpp
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

#ifndef TPMS_TREETOOLS
#define TPMS_TREETOOLS

#include <string>

#include <Bpp/Phyl/Tree/TreeTemplate.h>


namespace tpms{
    class TreeTools{
    public:
	static bpp::Node * newickToNode(std::string& in, unsigned int* jogger, bpp::Node* father, unsigned int& nodeId, bool distanceDelimitation);
	static bpp::TreeTemplate<bpp::Node> * newickToTree(std::string& in, bool distanceDelimitation=false);
	// extract the newick line from a file, stripping possible comments.
	static std::string extractNewickLineFromFile(std::ifstream& in);
	static void treePrint(bpp::Node * noeud, unsigned int cpt);
	static bool isBinaryTree(bpp::Node * node);
	static bool isAtLeastBinaryTree(bpp::Node * node);
	static std::string nodeToNewick(bpp::Node * node);
	static void multifurcated2binary(bpp::TreeTemplate< bpp::Node >* currTree, unsigned int parentID, std::vector< bpp::TreeTemplate< bpp::Node > * >& trees);
	static double getDistanceBetweenTwoNodes(bpp::Node * ancestor, bpp::Node * descendant);
	/**
	 * @brief Make the list of all the nodes contained in a subtree
	 * 
	 * @param node the root node of the subtree
	 * @param nodeListToFill the vector to be filled
	 * 
	 * @return a vector of all the nodes contained in this subtree
	 */
	static void getNodesOfTheSubtree(std::vector<bpp::Node*> &nodeListToFill, bpp::Node* node, bool leavesOnly = false, bpp::Node* ignoreSubtree = 00);
	
    static std::vector<unsigned int> toNodesIDs(std::vector<bpp::Node*> &nodes);
    static std::vector<bpp::Node*> toNodes(bpp::TreeTemplate<bpp::Node> &tree, std::vector<unsigned int> &IDs);
    
	static void destroySubtree(bpp::Node* node);
	static unsigned int depthOfANode(bpp::Node* node);
    };
}

#else

namespace tpms{
    class TreeTools;
}
#endif
