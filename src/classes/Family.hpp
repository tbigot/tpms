#ifndef FAMILY_HPP_TPMS
#define FAMILY_HPP_TPMS

#include <fstream>
#include <string>
#include <sstream>
#include <ostream>
#include <set>

//inclusions bio++
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTemplate.h>


//inclusions personnelles
#include "DataBase.hpp"
#include "Taxon.hpp"
#include "Waiter.hpp"


namespace tpms{
    
    enum NodeNature {ANY, DUPLICATION, SPECIATION};
    enum NodeType {NODE,LEAF};

class Family {
    
struct transfer {
	tpms::Taxon* donnor;
	tpms::Taxon* receiver;
	unsigned int perturbationIndex;
    };
    
private:
    
    bool doneMapping_NodeToTaxa;
    
    std::string name;
    DataBase * db;
    bpp::TreeTemplate<bpp::Node> * tree;
    std::map<std::string,tpms::Taxon*> mne2tax;
    std::set<tpms::Taxon*> taxa;
    /**
     * @brief contains the list of original leaves of the tree. Even after a subtree deletion, we can know weather a node matches to a sequence (which is associated to a species).
     */
    std::set<bpp::Node*> leaves;
    
    /**
     * @brief during a threaded work, all text results will be stored here
     */
    std::ostringstream results;
    
    
    bool containsUndefinedSequences;
    
    // before initialization: these objects are deleted after initialization
    std::stringstream *preamble;
    std::string *newick;
   
    
    std::vector<tpms::NodeNature> mapping_NodesToNatures;
    std::vector<float> mapping_NodesToUnicityScores;
    std::vector<Taxon *> mapping_NodesToTaxa;
    std::vector<unsigned int> mapping_NodesToTaxonomicShift;
    
     /**
     * @brief contains, after each step of doMapping_NodesToTaxonomicShift, the taxonomic affectation of a grandfather node of a node
     */
    std::vector<Taxon*> mapping_grandFatherWithoutThisNode;
    
    /**
     * @brief contains, after each step of doMapping_NodesToTaxonomicShift, the list of nodes with a taxonomic shift > 1
     */
    std::set<bpp::Node*> computed_nodesInducingPerturbation;
    
    /**
     * @brief returns true if at least one node of the tree is mapped to a taxonomic perturbation > 1
     */
    bool transfersRemaining();
    
    /**
     * @brief this mapping list is filled by the function doMapping_detectTransfers; It contains the list of each detected transfer
     */
    std::vector<transfer> computed_detectedTransfers;
        
    
    std::map<tpms::Taxon*, unsigned int> compute_UnicityScoreOnNode(std::vector<float> &scores, bpp::Node * node, bpp::Node * orignNode, bool virtualRootOnTheBranch = false);
    tpms::Taxon* mapNodeOnTaxon(bool virtualRootOnTheBranch, bool recordResults, bpp::Node& node, bpp::Node* origin=00, bool recursive=true, bpp::Node* ignoredNode=00);
    
    /**
     * @brief gives the difference of taxonomic assignation depth (in the species tree), of a grandfather of a certain node, removing this node
     * 
     * @param node the node we want to test
     * @return (init depth - final depth) of the grandfather node without the node given as parameter
     */
    unsigned int computeMappingShiftWithoutTheNode(bpp::Node* node);
    
    
    std::vector<bpp::Node*> getUnicityBestRoots(std::vector<bpp::Node *> nodes);
    std::vector<bpp::Node*> getTaxonomyBestRoots(std::vector<bpp::Node *> nodes);

    
    
    void writeRefTreeToFile(std::string path);
    void writeSpTreeToFile(std::string path);
    
    
    // utilities, unrelated to family, have to be moved to TreeTools
    void deleteFromLeavesToBif(bpp::Node * pnode);
    bpp::Node * removeUniqueSons(bpp::Node* localRoot);
    
    
    
public:
    void initialize();
    
    
    //constructeur à partir d'un fichier
    Family(std::stringstream* sIntro, std::string* sNewick, DataBase* dbp);
    
    std::string getName();
    bpp::TreeTemplate<bpp::Node> * getTree();
    
    bool getContainsUndefinedSequences();
    
    void doMapping_NodesToUnicityScores();
    void doMapping_LeavesToSpecies();
    void doMapping_NodesToNatures();
    void doMapping_NodesToTaxa();
    void doMapping_NodesToTaxonomicShift();
    
    void doRerooting_Unicity();
    void doRerooting_Taxonomy();
    void doRerooting_UnicityTaxonomy();

    void compute_detectTransfers();
    
    
    tpms::Taxon* getSpeciesOfNode(bpp::Node * node);
    NodeNature getNatureOfNode(bpp::Node * node);
    
    bool containsSpecie(tpms::Taxon* taxon);

    static void getLeavesFromNode(bpp::Node * pnode, std::vector< bpp::Node* >& leaves);
    static void getLeavesFromNode(bpp::Node * pnode, std::set< bpp::Node* >& leaves, int &leavesNumber);
    static void getLeavesFromNode(bpp::Node * pnode, std::set< bpp::Node* >& leaves);
    
    std::set<std::string> getLeavesNamesFromNode(bpp::Node * pnode);
    void writeTreeToStream(bpp::Node* root, std::ostream& sortie, unsigned int deep);
    
    
    std::vector<float> &getUnicityScores();
    
    int numberOfNodesBetween(bpp::Node * ancestor, bpp::Node * pnode);
    
    void loadRefTreeFromFile(std::string path);
    void loadSpTreeFromFile(std::string path);
    
    void addSequencesNames(bpp::Node* currNode);
    
    void getTaxaOnThisSubtree(bpp::Node * node, std::vector<tpms::Taxon*>& speciesList);
    
    std::set<tpms::Taxon *> &getSpecies();
    
    /**
     * @brief Replace names of the pieces of pattern trees nodes by the name of the sequences
     * 
     * @param currNode node to start with
     */
    void labelWithSequencesNames(bpp::Node* currNode);
    
    NodeNature getNatureOf(bpp::Node* node);
    
    
    
    static void threadWork_initialize(Waiter *waiter, boost::mutex *waiterMutex, std::vector<tpms::Family*>::iterator &familiesBegin, std::vector<tpms::Family*>::iterator &familiesEnd);
    
    static void threadedWork_launchJobs(std::vector<Family *> families, void (Family::*function)(), unsigned int nbThreads, std::ostream *output = 00);
    static void threadedWork_oneThread(void(Family::*function)(),Waiter *progressbar, boost::mutex *progressbarMutex, std::ostream *output, boost::mutex *outputMutex , std::vector<tpms::Family*>::iterator &currPartBegin, std::vector<tpms::Family*>::const_iterator &currPartEnd);
    
    std::ostringstream &threadedWork_getResults();
    
};}

#else
namespace tpms{

class Family;
}
#endif
