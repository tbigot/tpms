//
// File: DataBase.hpp
// Created by: Thomas Bigot
//

/*
 *   Copyright or © or Copr. Thomas Bigot 2012
 * 
 *   This software is a computer program whose purpose is to provide classes
 *   for phylogenetic data analysis.
 * 
 *   This software is governed by the CeCILL  license under French law and
 *   abiding by the rules of distribution of free software.  You can  use,
 *   modify and/ or redistribute the software under the terms of the CeCILL
 *   license as circulated by CEA, CNRS and INRIA at the following URL
 *   "http://www.cecill.info".
 * 
 *   As a counterpart to the access to the source code and  rights to copy,
 *   modify and redistribute granted by the license, users are provided only
 *   with a limited warranty  and the software's author,  the holder of the
 *   economic rights,  and the successive licensors  have only  limited
 *   liability.
 * 
 *   In this respect, the user's attention is drawn to the risks associated
 *   with loading,  using,  modifying and/or developing or reproducing the
 *   software by the user in light of its specific status of free software,
 *   that may mean  that it is complicated to manipulate,  and  that  also
 *   therefore means  that it is reserved for developers  and  experienced
 *   professionals having in-depth computer knowledge. Users are therefore
 *   encouraged to load and test the software's suitability as regards their
 *   requirements in conditions enabling the security of their systems and/or
 *   data to be ensured and,  more generally, to use and operate it in the
 *   same conditions as regards security.
 * 
 *   The fact that you are presently reading this means that you have had
 *   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef DATABASE_HPP_TPMS
#define DATABASE_HPP_TPMS


#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <map>

#include <boost/thread.hpp>

#include "Family.hpp"
#include "Taxon.hpp"
#include "Waiter.hpp"

//inclusions bio++
#include <Bpp/Phyl/Tree/Tree.h>



namespace tpms{
  class DataBase {

  private:

    bool mappingDone_LeavesToSpecies_;
    bool mappingDone_NodesToTaxa_;
    bool mappingDone_NodesToUnicityScores_;
    bool mappingDone_NodesToMaxDepth_;

    bool expectSynonyms_;
    
    unsigned int nbThreads_;

    
    std::map<std::string,tpms::Taxon*> taxa_;
    std::map<unsigned int,tpms::Taxon*> nodeToTaxon_;
    unsigned int nbFamilies_;
    bool reconciled_;
    std::string filename_;
    bpp::TreeTemplate<bpp::Node> * speciesTree_;
    std::vector<Family *> families_;
    std::map<std::string,int> nbFC_; // cache de recherche de la fonction nbFamiliesContaining
    std::set<std::string> species_;

    boost::mutex m_;
    boost::mutex waiterMutex_;

    void loadFamily_(std::stringstream * intro, std::string newick);
    void loadSpeciesTree_(std::string newickLine);
    void loadFromFile_(std::ifstream& RAPfile);
    void renameNodes_(bpp::TreeTemplate<bpp::Node> *);





  public:


    /**
     * @brief Function to display status of the collection: mapping that have been done, others that are not yet.
     * @return a string containing the text to display.
     **/

    std::string getStatus(std::string);

    //constructeur à partir d'un fichier
    DataBase(std::string path, bool expectSynonyms, unsigned int nbThreads = 1);
    bpp::TreeTemplate<bpp::Node> * getSpeciesTree();
    std::vector<tpms::Family *> & getFamilies();
    unsigned int getNbFamilies();
    std::string getParentTaxon(std::string pTaxon, unsigned int level);



    //nombre de familles contenant un taxon donné
    // int nbFamiliesContaining(std::string pTax);

    std::set<std::string> getAllNodes(bpp::Node * localRoot,bool nodesWanted = true);

    static void doFamiliesMapping_LeavesToSpecies_oneThread(Waiter * waiter, boost::mutex *waiterMutex, std::vector<tpms::Family*>::iterator &familiesBegin, std::vector<tpms::Family*>::iterator &familiesEnd);


    void doFamiliesMapping_LeavesToSpecies();
    void doFamiliesMapping_NodesToTaxa();
    void doFamiliesMapping_NodesToMaxDepth();
    void doFamiliesMapping_NodesToUnicityScores();
    void doFamiliesMapping_NodesToTaxonomicShift();

    void doFamiliesRerooting_Unicity();
    void doFamiliesRerooting_Taxonomy();
    void doFamiliesRerooting_UnicityTaxonomy();
    void doFamiliesRerooting_LessTransfers(std::ostream * output);
    void doFamiliesRerooting_Daubin();


    void doFamiliesComputation_detectTransfers(std::ostream * output);

    std::set<std::string> * getSpecies();
    tpms::Taxon* nameToTaxon(std::string taxonName);
    tpms::Taxon* nodeIdToTaxon(unsigned int nodeId);


  };}
  #else
  
  namespace tpms{
    class DataBase;
  }
  
  #endif
  