#ifndef FAMILY_HPP_TPMS
#define FAMILY_HPP_TPMS

#include <fstream>
#include <string>
#include <sstream>
#include <set>

//inclusions bio++
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTemplate.h>


//inclusions personnelles
#include "NodeConstraints.hpp"
#include "DataBase.hpp"
#include "Taxon.hpp"
namespace tpms{

class Family {
    
private:
    
    enum NodeNature {ANY, DUPLICATION, SPECIATION};
    enum Type {NODE,LEAF};
    
    std::string name;
    DataBase * db;
    bpp::TreeTemplate<bpp::Node> * tree;
    std::map<std::string,tpms::Taxon*> mne2tax;
    std::set<tpms::Taxon*> species;
   
    
    std::vector<tpms::Taxon*> mapping_LeavesToSpecies;
    std::vector<NodeNature> mapping_NodesToNatures;
    std::vector<unsigned int> mapping_NodesToUnicityScores;
    std::vector<Taxon *> mapping_NodesToTaxa;
        
    
    std::map<tpms::Taxon*, unsigned int> compute_UnicityScoreOnNode(std::vector<unsigned int> &scores, bpp::Node * node, bpp::Node * orignNode);

    
    void mapNodeOnTaxon(bpp::Node& node);
    
    void writeRefTreeToFile(std::string path);
    void writeSpTreeToFile(std::string path);
    
    
    // utilities, unrelated to family, has to be moved to TreeTools
    void deleteFromLeavesToBif(bpp::Node * pnode);
    bpp::Node * removeUniqueSons(bpp::Node* localRoot);   
    
    
public:
    
    
    //constructeur à partir d'un fichier
    Family(std::stringstream* sIntro, std::string sNewick, DataBase* dbp);
    
    std::string getName();
    bpp::TreeTemplate<bpp::Node> * getTree();
    
    
    void doMapping_NodesToUnicityScores();
    void doMapping_NodesToBestUnicityScores();
    void doMapping_LeavesToSpecies();
    void doMapping_NodesToNatures();
    void doMapping_NodesToTaxa();
    
    
    tpms::Taxon* getSpeciesOfNode(bpp::Node * node);
    NodeNature getNatureOfNode(bpp::Node * node);
    
    bool containsSpecie(tpms::Taxon* taxon);

    static void getLeavesFromNode(bpp::Node * pnode, std::vector< bpp::Node* >& leaves);
    static void getLeavesFromNode(bpp::Node * pnode, std::set< bpp::Node* >& leaves, int &leavesNumber);
    static void getLeavesFromNode(bpp::Node * pnode, std::set< bpp::Node* >& leaves);
    
    std::set<std::string> getLeavesNamesFromNode(bpp::Node * pnode);
    void writeTreeToStream(bpp::Node* root, std::ostream& sortie, unsigned int deep);
    
    
    std::vector<unsigned int> &getUnicityScores();
    
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
    
    
};}

#else
namespace tpms{

class Family;
}
#endif
