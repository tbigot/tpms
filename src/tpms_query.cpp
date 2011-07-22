#include <iostream>
#include <fstream>
#include <string>

#include "classes/DataBase.hpp"
#include "classes/Pattern.hpp"
#include "classes/CmdLineArgs.hpp"
#include "classes/CandidateNode.hpp"
#include "classes/TreeTools.hpp"
#include "classes/Waiter.hpp"


#include <Bpp/Phyl/TreeTemplate.h>

using namespace std;
using namespace tpms;
using namespace bpp;

int main(int argc, char *argv[]) {
   try{
    CmdLineArgs args(argc, argv, "database,output",cerr);
    args.print(cout);
	istringstream nbThreadsSS(args.getArg("threads"));
	unsigned int nbThreads;
	nbThreadsSS >> nbThreads;
	if(nbThreads < 1 || nbThreads > 128) nbThreads = 1;
	cout << nbThreads << " threads will be used." << endl;
	DataBase dbTest(args.getArg("database"),nbThreads);
	dbTest.doFamiliesMapping_LeavesToSpecies();
	dbTest.doFamiliesMapping_NodesToTaxa();
	dbTest.doFamiliesMapping_NodesToUnicityScores();

	
	
	
	
	string fileName;
	cout << "\nPlease enter the name in which the results will be recorded\nFILENAME> " << flush;
	getline(cin,fileName);

	
	while(fileName != ""){
		ofstream ooo(string(args.getArg("output")+"/"+fileName).c_str(), ofstream::out);
		string queryLine;
		cout << "Now enter the pattern in pseudo-newick formalism (please read the doc)\nPATTERN> " << flush;
		getline(cin,queryLine);
		
		ooo << "; PATTERN ENTERED :" << queryLine << endl;
		
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
			
			
			
			for(vector<TreeTemplate<Node> * >::iterator currpt = trees.begin(); currpt != trees.end(); currpt++){
			    Pattern tp(**currpt,dbTest);
			    cout << "\nPerforming tree pattern matching search in the collection, please wait" << endl;
			    tp.search(dbTest.getFamilies(),familles);
			}
			
			
			unsigned int trouve = 0;
			// DÉPOUILLAGE DES RÉSULTATS
			cout << "\n   Writing results to the file " << args.getArg("output")+"/"+fileName << endl;
			
			cout << "   - Building all matching subtrees, this could take a while..." << endl;
			
			Waiter waitB(&cout, familles.size(),'-');
			
			for(vector<pair<Family *,CandidateNode * > >::iterator curFam = familles.begin(); curFam != familles.end(); curFam++) {
				waitB.step();
				
				ooo << "; family_name " << curFam->first->getName() << endl;
				ooo << "; number_of_species " << curFam->first->getSpecies().size() << endl;
				
				// retrieving all trees corresponding to CandidateNode :
				vector<TreeTemplate<Node> *> resultTrees;
				
				if(args.getArg("onlyFamNames").empty()) if(args.getArg("pleaseReturnWholeMatchingSubtrees").empty()){
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
			ooo << endl;
			waitB.drawFinal();
			cout << "   -[Done]\n" << endl;
			ooo.close();
			cout << "SELECTED: " << familles.size() << " / " << dbTest.getNbFamilies() << " families in the database." << endl;
			cout << trouve << " patterns found." << endl;
			
		}catch (bpp::NodeNotFoundException e) {
			cout << "Node not found: \n" << e.what() << endl;
		}
		
		cout << "\n\n\nquery name> " << flush;
		getline(cin,fileName);
	}
    }catch(string missingArguments){
	cout << "Missing arguments: " << missingArguments.substr(0,missingArguments.size()-2) << '.' << endl;
	return(1);
    }
}
