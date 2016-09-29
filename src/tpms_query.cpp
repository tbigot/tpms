//
// File: tpms_query.cpp
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

#include <iostream>
#include <fstream>
#include <string>

#include "classes/DataBase.hpp"
#include "classes/Pattern.hpp"
#include "classes/CandidateNode.hpp"
#include "classes/TreeTools.hpp"
#include "classes/Waiter.hpp"


#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/App/BppApplication.h>

using namespace std;
using namespace tpms;
using namespace bpp;

int main(int argc, char *argv[]) {
  BppApplication thisApp(argc, argv, "tpms Query");
  try{
    string msg25 = "     version 2 beta";
    cout << "      __                                 \n     /\\ \\__   " << msg25 << "\n     \\ \\ ,_\\  _____    ___ ___     ____\n      \\ \\ \\/ /\\ '__`\\/' __` __`\\  /',__\\\n       \\ \\ \\_\\ \\ \\L\\ \\\\ \\/\\ \\/\\ \\/\\__, `\\\n        \\ \\__\\\\ \\ ,__/ \\_\\ \\_\\ \\_\\/\\____/\n         \\/__/ \\ \\ \\/ \\/_/\\/_/\\/_/\\/___/\n                \\ \\_\\\n                 \\/_/    W E L C O M E\n" << endl;
    string nbThreadsS = thisApp.getParam("threads");
    if (nbThreadsS.empty()) nbThreadsS = "1";
    istringstream nbThreadsSS(nbThreadsS);
    unsigned int nbThreads;
    nbThreadsSS >> nbThreads;
    if(nbThreads < 1) nbThreads = 1;
    cout << nbThreads << " threads will be used." << endl;
    // Do we have to use synonyms?
    bool expectSynonyms = false;
    if(thisApp.getParam("synonyms") == "yes")
      expectSynonyms = true;
    DataBase dbTest(thisApp.getParam("collection"),expectSynonyms,nbThreads);
    dbTest.doFamiliesMapping_LeavesToSpecies();
    
    
    string fileName;
    cout << "\nPlease enter the name in which the results will be recorded\nFILENAME> " << flush;
    getline(cin,fileName);
    
    
    while(fileName != ""){
      ofstream ooo(string(thisApp.getParam("output-dir")+"/"+fileName).c_str(), ofstream::out);
      string queryLine;
      cout << "Now enter the pattern in pseudo-newick formalism (please read the doc)\nPATTERN> " << flush;
      getline(cin,queryLine);
      
      ooo << "; PATTERN ENTERED :" << queryLine << endl;
      
      if(queryLine.find(';') == string::npos){
        cout << "The pattern you entered does not contain a semi-colon character (;) are you sure it is complete? Please do it again." << endl;
      }
      
      vector<pair<Family *,CandidateNode *> > familles;
      try{
        // REQUÊTE
        TreeTemplate<Node> * initTree = tpms::TreeTools::newickToTree(queryLine);
        
        if(!tpms::TreeTools::isBinaryTree(initTree->getRootNode())){
          cout << "\n**** Non binary query tree detected!" << endl;
          if(tpms::TreeTools::isAtLeastBinaryTree(initTree->getRootNode()))
            cout << "****   Binary trees will be generated instead." << endl;
          else
            cout << "****   Unique sons detected: imminent fail, sorry." << endl;
        }
        
        vector<TreeTemplate<Node> *> trees;
        trees.push_back(initTree);
        tpms::TreeTools::multifurcated2binary(initTree,initTree->getRootId(),trees);
        if(trees.size() > 1)
          cout << "Finaly, " << trees.size() << " patterns have been generated." << endl;
        
        bool isOk = true;
        
        for(vector<TreeTemplate<Node> * >::iterator currpt = trees.begin(); currpt != trees.end(); currpt++){
          Pattern tp(**currpt,dbTest);
          if (!tp.isOk()){
            cout << "\n\n [!!!] Your pattern contains invalid species / taxa: unable to perform the search. Please retry."<< endl;
            isOk = false;
            break;
          }
          cout << "\nPerforming tree pattern matching search in the collection, please wait" << endl;
          
          dbTest.doFamiliesMapping_NodesToTaxa();
          
          dbTest.doFamiliesMapping_NodesToMaxDepth();
          
          Family::threadedWork_patternMatching(dbTest.getFamilies(),&tp,&familles,nbThreads);
          
          //tp.search(dbTest.getFamilies(),familles);
        }
        
        if(isOk){
          unsigned int trouve = 0;
          // DÉPOUILLAGE DES RÉSULTATS
          cout << "\n   Writing results to the file " << thisApp.getParam("output-dir")+"/"+fileName << endl;
          
          cout << "   - Building all matching subtrees, this could take a while..." << endl;
          
          Waiter waitB(&cout, familles.size(),'-');
          
          for(vector<pair<Family *,CandidateNode * > >::iterator curFam = familles.begin(); curFam != familles.end(); curFam++) {
            waitB.doStep();
            
            ooo << "; family_name " << curFam->first->getName() << endl;
            ooo << "; number_of_species " << curFam->first->getSpecies().size() << endl;
            
            // retrieving all trees corresponding to CandidateNode :
            vector<TreeTemplate<Node> *> resultTrees;
            
            if(thisApp.getParam("onlyFamNames").empty())
            {
              if(thisApp.getParam("extractSubtree").empty()){
                trouve += curFam->second->genTrees(resultTrees);
                for(vector<TreeTemplate<Node> *>::iterator currTree = resultTrees.begin(); currTree != resultTrees.end(); currTree++){
                  curFam->first->addSequencesNames((*currTree)->getRootNode());
                  ooo << curFam->first->getName() << '\n';
                  ooo << tpms::TreeTools::nodeToNewick((*currTree)->getRootNode()) << ";\n";
                }
                
              } else {
                trouve += curFam->second->getWholeMatchingSubtrees(resultTrees);
                for(vector<TreeTemplate<Node> *>::iterator currTree = resultTrees.begin(); currTree != resultTrees.end(); currTree++){
                  curFam->first->labelWithSequencesNames((*currTree)->getRootNode());
                  ooo << curFam->first->getName() << '\n';
                  ooo << tpms::TreeTools::nodeToNewick((*currTree)->getRootNode()) << ";\n";
                }
              }
            }
            
            
          }
          ooo << endl;
          waitB.drawFinal();
          cout << "   -[Done]\n" << endl;
          ooo.close();
          cout << "SELECTED: " << familles.size() << " / " << dbTest.getNbFamilies() << " families in the collection." << endl;
          cout << trouve << " patterns found." << endl;
        } 
        
      }catch (bpp::NodeNotFoundException e) {
        cout << "Node not found: \n" << e.what() << endl;
      }
      
      cout << "\n\n\nFILENAME> " << flush;
      getline(cin,fileName);
    }
  }catch(Exception missingArguments){
    cout << "Missing arguments: " << missingArguments.what() << '.' << endl;
    return(1);
  }
}
