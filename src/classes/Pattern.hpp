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
    
public:
    
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