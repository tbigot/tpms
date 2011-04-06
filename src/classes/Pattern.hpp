#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <vector>

#include "NodeConstraints.hpp"
#include "Family.hpp"
#include "DataBase.hpp"
#include "CandidateNode.hpp"
#include "Taxon.hpp"

//inclusions bio++
#include <Bpp/Phyl/Tree.h>

class Pattern {
    
private:
    /**
     * @brief The database associated to this pattern.
     **/
    DataBase &db;
    
    
     /**
     * @brief The pattern tree itself.
     **/
    bpp::TreeTemplate<bpp::Node> &tree;
    
    /**
     * @brief All species that can be in this pattern.
     **/
    std::set<tpms::Taxon*> allSpecies;
    
    
    /**
     * @brief Associates a set of authorized species to each pattern node.
     **/
    std::vector<tpms::Taxon*> species;
    

    /**
     * @name Construction functions
     * @{
     * @brief Associates a set of authorized species to each pattern node
     * OLD DOC: cette fonction associe, pour chaque feuille du pattern, une liste d'espèces
     * quand on va entamer une recherche, et qu'il y a un taxon dans le pattern, il faut savoir exactement quelles espèces de l'arbre il recouvre. Si un nœud de l'abre est « BACTERIA », on doit avoir la liste de toutes les espèces correspondantes, c'est-à-dire toutes les feulles qui sont sous BACTERIA dans l'abre des espèces.
     * species est donc une map qui associe à chaque ID de nœud un set d'espèces (liste ordonnée) sous forme de chaînes.
     *
     **/
    void fillSpeciesFromLeavesNames();
    
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
    bool patternMatch(Family& family ,bpp::Node* target, bpp::Node* pattern, tpms::CandidateNode* fatherCandidate);
    
    
    /**
     * @name Small utilities
     * @{
     *
     **/
    static bool isLeaf(bpp::Node * pNode);
    
    /**
     * @brief Does the node pNode accept the species parameter? Does species belong to the authorised species for this node?
     *
     **/
    bool pNodeAuthorisesThisSpecies(bpp::Node * pNode, tpms::Taxon* species);
    
    
    /** @} */
    
    std::vector<int> getIdWithTaxaList(bpp::TreeTemplate<bpp::Node> * sTree, std::set<std::string> * taxa);
    int enumerateTaxon(bpp::Node * localRoot, std::set<std::string> * tsMembers, bool tsComplementary, bpp::Node * impasse, unsigned int * totsize);
    
    std::vector<NodeConstraints *> constraints;
    void print(bpp::Node * noeud, int spc);
    
    bool isTreeBinary();
    
public:
    
    Pattern(bpp::TreeTemplate< bpp::Node >& tree, DataBase& pRefDB);
    ~Pattern();
    
    
    static bool nodeOnlyContainsTheseTaxa(bpp::Node * localRoot, std::set<tpms::Taxon*> & taxonMembers, bool invert=false);
    
    static bool isUnder(bpp::Node * descendant, bpp::Node * ascendant);
    unsigned int search(std::vector< Family* >& families, std::vector< std::pair< Family*, tpms::CandidateNode* > >& result);
    std::vector<int> xferDetected(std::map<int,Family *> * families, std::vector<int> famList, std::string sourceTaxon, std::string targetTaxon, std::string monophylyTaxon, unsigned int verifDeep, std::vector<unsigned int> bootstraps,  std::vector<unsigned int> sourceRates);
    std::vector<int> xferGapDetect(std::map<int,Family *> families, std::set<std::string>);
};
