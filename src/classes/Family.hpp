//
// File: Family.hpp
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

#ifndef FAMILY_HPP_TPMS
#define FAMILY_HPP_TPMS

#include <fstream>
#include <string>
#include <sstream>
#include <ostream>
#include <set>

//inclusions bio++
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>

namespace tpms{
    
    enum NodeNature {ANY, DUPLICATION, SPECIATION};
    enum NodeType {NODE,LEAF};
}

//inclusions personnelles
#include "DataBase.hpp"
#include "Taxon.hpp"
#include "Waiter.hpp"
#include "CandidateNode.hpp"
#include "Pattern.hpp"


namespace tpms{
class Family {
    
struct transfer {
	tpms::Taxon* donnor;
        std::vector<std::string> donnorLeaves;
	std::vector<tpms::Taxon*> acceptors;
        std::vector<std::vector<std::string> > acceptorsLeaves;
	unsigned int perturbationIndex;
        unsigned int bootstrap;
        tpms::Taxon* formerAncestorTaxon;
    };
    
private:
    
    unsigned int highestID_;
    
    bool doneMapping_NodeToTaxa_;
    
    std::string name_;
    DataBase * db_;
    bpp::TreeTemplate<bpp::Node> * tree_;
    std::map<std::string,tpms::Taxon*> mne2tax_;
    std::set<tpms::Taxon*> taxa_;
    /**
     * @brief contains the list of original leaves of the tree. Even after a subtree deletion, we can know weather a node matches to a sequence (which is associated to a species).
     */
    std::set<bpp::Node*> leaves_;
    
    /**
     * @brief during a threaded work, all text results will be stored here
     */
    std::ostringstream results_;
    
    unsigned int getTaxonomicSum_(std::vector<Taxon*> &taxa);
    
    
    bool containsUndefinedSequences_;
    
    // before initialization: these objects are deleted after initialization
    std::stringstream *preamble_;
    std::string *newick_;
   
    
    std::vector<tpms::NodeNature> mapping_NodesToNatures_;
    std::vector<float> mapping_NodesToUnicityScores_;
    std::vector<Taxon *> mapping_NodesToTaxa_;
    std::vector<unsigned int> mapping_NodesToTaxonomicShift_;
    std::vector<unsigned int> mapping_NodesToMaxDepths_;
    
    
    std::vector<std::set<tpms::Taxon*> > cacheMapping_SubtreesToTaxalist_;
    
     /**
     * @brief contains, after each step of doMapping_NodesToTaxonomicShift, the taxonomic affectation of a grandfather node of a node
     */
    std::vector<Taxon*> mapping_grandFatherWithoutThisNode_;
    
    /**
     * @brief contains, after each step of doMapping_NodesToTaxonomicShift, the list of nodes with a taxonomic shift > 1
     */
    std::set<bpp::Node*> computed_nodesInducingPerturbation_;
    
    /**
     * @brief returns true if at least one node of the tree is mapped to a taxonomic perturbation > 1
     */
    bool transfersRemaining_();
    
    /**
     * @brief this mapping list is filled by the function doMapping_detectTransfers; It contains the list of each detected transfer
     */
    std::vector<transfer> computed_detectedTransfers_;
        
    
    std::map<tpms::Taxon*, unsigned int> compute_UnicityScoreOnNode_(std::vector<float> &scores, bpp::Node * node, bpp::Node * orignNode, bool virtualRootOnTheBranch = false);
    
    void computeTaxonomicShift_(bpp::Node* node, bpp::Node* father, bpp::Node* grandFather, bpp::Node* greatGrandFather, std::vector<tpms::Taxon*>* local_nodeToTaxon, std::vector<tpms::Taxon*>* local_GFmappingWithoutTheNode, std::vector<unsigned int>* local_perturbationsInduced, std::set<bpp::Node*>* local_nodesInducingPerturbations);
    
    tpms::Taxon* mapNodeOnTaxon_(std::vector< tpms::Taxon* >* referenceMapping, std::vector< tpms::Taxon* >* recordResult, bpp::Node* node, bpp::Node* origin = 00, std::set< bpp::Node* >* ignoredNodes = 00, bool recursive = true, bpp::Node* ignoredNode = 00);
    unsigned int getTaxonomicSum_(bpp::Node* node, bpp::Node* father, std::vector<tpms::Taxon*>* local_nodesToTaxa, std::set<bpp::Node*>* nodesToIgnore);
    
    std::vector<bpp::Node*> getUnicityBestRoots_(std::vector<bpp::Node *> nodes);
    std::vector<bpp::Node*> getTaxonomyBestRoots_(std::vector<bpp::Node *> nodes);
    std::vector<unsigned int> getLessTransfersBestRoots_(std::vector<unsigned int> nodes);
    std::vector<unsigned int> getDaubinCriteriaBestRoots_(std::vector<unsigned int> nodes);
  
    
    void writeRefTreeToFile_(std::string path);
    void writeSpTreeToFile_(std::string path);
    
    
    // utilities, unrelated to family, have to be moved to TreeTools
    void deleteFromLeavesToBif_(bpp::Node * pnode);
    bpp::Node * removeUniqueSons_(bpp::Node* localRoot);
    
    void reRootAt_(std::vector<bpp::Node*> bestRoots);
    void reRootAt_(std::vector<unsigned int> bestRoots);

    void atomizeTaxon_(std::vector< tpms::Taxon* >& acceptors, std::vector< std::vector<std::string> >& acceptorsLeaves, tpms::Taxon* acceptor, bpp::Node* subtree);
    
    void print_tax_tree_(bpp::Node& node, unsigned int depth, bpp::Node* origin, std::vector<tpms::Taxon*>* mapping, std::set<bpp::Node*>* ignoredNodes, bool subtreeIgnored);
    
    void updateTaxa_();
    
    unsigned int mapNodeToMaxDepth_(bpp::Node*);
    
    unsigned int compute_detectTransfers_(bool writeResults,unsigned int &gainSum);

    
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
    void doMapping_NodesToMaxDepth();
    
    void doRerooting_Unicity();
    void doRerooting_LessTransfers();
    void doRerooting_Taxonomy();
    void doRerooting_UnicityTaxonomy();
    void doRerooting_Daubin();
    

    void compute_detectTransfers();
    
    
    tpms::Taxon* getTaxonOfNode(bpp::Node * node);
    unsigned int getMaxDepthOfSubtree(bpp::Node * node);

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
    
    std::set<tpms::Taxon*> & getTaxaOnThisSubtree(bpp::Node * node);
    
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
    static void threadedWork_patternMatching(std::vector<Family *> & families, tpms::Pattern * pattern, std::vector<std::pair<tpms::Family*, tpms::CandidateNode*> > * matchingFamilies, unsigned int nbThreads);
    static void threadedWork_oneThread(void(Family::*function)(),Waiter *progressbar, boost::mutex *progressbarMutex, std::ostream *output, boost::mutex *outputMutex , std::vector<tpms::Family*>::iterator &currPartBegin, std::vector<tpms::Family*>::const_iterator &currPartEnd);
    static void threadedWork_onePMThread(tpms::Pattern * pattern,std::vector<std::pair<tpms::Family *, tpms::CandidateNode*> > * resultVector,Waiter *progressbar, boost::mutex *progressbarMutex, boost::mutex *outputMutex , std::vector<tpms::Family*>::iterator &currPartBegin, std::vector<tpms::Family*>::const_iterator &currPartEnd);
    
    std::ostringstream &threadedWork_getResults();
    
    void clearCache();
    void initCache();
};}

#else
namespace tpms{

class Family;
}
#endif
