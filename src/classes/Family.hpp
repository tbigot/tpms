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
#include "DataBase.hpp"
#include "Taxon.hpp"
#include "NodeConstraints.hpp"

class Family {
    
private:
    std::string name;

    
    DataBase * db;
    std::vector<tpms::Taxon*> leave2spe;
    bpp::TreeTemplate<bpp::Node> * tree;
    std::map<std::string,tpms::Taxon*> mne2tax;
    std::set<tpms::Taxon*> species;
    std::vector<NodeConstraints::NodeNature> node2nat;
    bpp::TreeTemplate<bpp::Node> * refTree;
    std::vector<unsigned int> unicityScores;
    
    void renameNodes(bpp::TreeTemplate<bpp::Node> *);
    
    
    void deleteFromLeavesToBif(bpp::Node * pnode);
    
    std::map<std::string, unsigned int> computeUnicity(std::vector<unsigned int> &scores, bpp::Node * node, bpp::Node * originNode);
    
    bpp::Node * removeUniqueSons(bpp::Node* localRoot);
    
    std::string mapNodeOnTaxon(bpp::Node& nodep);
    
    void writeRefTreeToFile(std::string path);
    void writeSpTreeToFile(std::string path);
    
    
   
    
    
public:
    //constructeur à partir d'un fichier
    Family(std::stringstream* sIntro, std::string sNewick, DataBase* dbp);
    
    bpp::TreeTemplate<bpp::Node> * getTree();
    void genUnicityScores();
    void genBestUnicityScores();
    void genLeaveToSpecies();
    
    tpms::Taxon* getSpeciesOfNode(bpp::Node * node);
    NodeConstraints::NodeNature getNatureOfNode(bpp::Node * node);
    
    bool containsSpecie(tpms::Taxon* taxon);
    // bool containsSpecies(std::set<std::string> speciesList);
    // std::set<std::string> * getSpecies();
    static void getLeavesFromNode(bpp::Node * pnode, std::vector< bpp::Node* >& leaves);
    static void getLeavesFromNode(bpp::Node * pnode, std::set< bpp::Node* >& leaves, int &leavesNumber);
    static void getLeavesFromNode(bpp::Node * pnode, std::set< bpp::Node* >& leaves);
    
    std::set<std::string> getLeavesNamesFromNode(bpp::Node * pnode);
    std::string getName();
    void writeTreeToStream(bpp::Node* root, std::ostream& sortie, unsigned int deep);
    
    //! Génération du Vrai arbre
    /*!
     *	Ce générateur va partir d'un arbre vrai de référence et n'en garder que les espèces présentes dans la famille
     */
    void genRefTree(bool save=true, std::string path="");
    
    std::vector<unsigned int> &getUnicityScores();
    
    int numberOfNodesBetween(bpp::Node * ancestor, bpp::Node * pnode);
    
    void loadRefTreeFromFile(std::string path);
    void loadSpTreeFromFile(std::string path);
    
    void addSequencesNames(bpp::Node* currNode);
    
    /**
     * @brief Replace names of the pieces of pattern trees nodes by the name of the sequences
     * 
     * @param currNode node to start with
     */
    void labelWithSequencesNames(bpp::Node* currNode);
    
    
};

#else

class Family;

#endif
