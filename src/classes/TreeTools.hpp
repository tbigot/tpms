#ifndef TPMS_TREETOOLS
#define TPMS_TREETOOLS

#include <string>

#include <Phyl/Tree.h>
#include <Phyl/Newick.h>
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
};
}

#else

namespace tpms{
class TreeTools;
}
#endif
