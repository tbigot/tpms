#ifndef TPMS_CANDIDATENODE
#define TPMS_CANDIDATENODE

#include <set>
#include <map>

#include <Bpp/Phyl/TreeTemplate.h>

namespace tpms{
    class CandidateNode{
	private:
	    bpp::Node * treeNode;
	    bpp::Node * patternNode;
	    tpms::CandidateNode * father;
	    std::map<bpp::Node *, std::vector<tpms::CandidateNode *> > sons;
	    bool isRoot;
	    unsigned int addsons;
	    bool isLeaf();
	    
	    /**
	    * @brief Recursively builds result trees from CandidateNode trees
	    * 
	    * @return the nodes of trees that can be built from this CandidateNode
	    */
	    std::vector<bpp::Node *> recGenTrees();
	    
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
