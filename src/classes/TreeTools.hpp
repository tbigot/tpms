#ifndef TPMS_TREETOOLS
#define TPMS_TREETOOLS

#include <string>

#include <Bpp/Phyl/TreeTemplate.h>


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
	static void getNodesOfTheSubtree(std::vector<bpp::Node*> &nodeListToFill,bpp::Node* node);
	
	static void destroySubtree(bpp::Node* node);
	static unsigned int depthOfANode(bpp::Node* node);
    };
}

#else

namespace tpms{
    class TreeTools;
}
#endif
