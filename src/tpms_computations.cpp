//
// File: tpms_computations.cpp
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


/*
 * The aim of this program is to provide a command line style interface to 
 * perform some scoring and re-rooting operations on collections.
 * 
 * Usage is:
 * tpms_computations -collection=RAP-File -threads=N
 * 
 * threads
 * 	The number of threads used by the program during specific multithreaded operations.
 * 
 * collection
 * 	A file in RAP format containing the collection of trees.
 * 
 * */


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


enum annotationsTypes {NONE,UNICITY,BIPARTITION};

void print_menu();
std::string get_choice();
void print_status(tpms::DataBase &db);
void save_to_file(tpms::DataBase& collection, std::string path, annotationsTypes type);

void print_menu(){
    cout << "\n  +----------------------------\n  |        MENU" << endl;
    cout << "  |\n  |SCORING\n  |=======\n  |   CU) compute unicity scores\n  |   CB) compute bipartition scores"<< endl;
    cout << "  |\n  |ROOTING\n  |==========\n  |   RU) root with unicity criteria\n  |   RT) root with taxonomic criteria\n  |   RC) root with the unicity+taxonomy criteria\n  |   RXD) root with the minimum of transfers criteria\n  |   RVD) root Daubin's criteria" << endl;
    cout << "  |\n  |MISC\n  |====\n  |    T) collection status\n  |   XD) LGT detection\n  |    M) help menu\n  |    S) Save collection\n  |   SU) Save collection with Unicity scores"<< endl;
    cout << "  |    Q) QUIT"<< endl;
    
}

void print_status(DataBase& db)
{
    cout << "\n   +-------------------------------------+" << endl;
    cout << db.getStatus("   |");
    cout << "   +-------------------------------------+" << endl;
}


std::string get_choice(){
    cout << "\ncommand? (m for help) > " << flush;
    string answer;
    getline(cin,answer);
    std::transform(answer.begin(), answer.end(), answer.begin(), (int(*)(int)) toupper);
    return(answer);
}

void save_to_file(DataBase& collection, std::string path, annotationsTypes type)
{
    ofstream out(path.c_str(),ofstream::out);
    // liste des familles
    vector<Family*> & dbFamilies = collection.getFamilies();
    Waiter patienteur(&std::cout, dbFamilies.size(), '#');
    for(vector<Family*>::iterator currFamily = dbFamilies.begin(); currFamily != dbFamilies.end(); currFamily++){
	patienteur.step();
	
	vector<float> scores;
	if(type == UNICITY){
	    (*currFamily)->doMapping_NodesToUnicityScores();
	    scores = (*currFamily)->getUnicityScores();}
	
	// re labelling the node names
	if(type != NONE){
	    vector<Node *> nodesList = (*currFamily)->getTree()->getNodes();
	    for(vector<Node *>::iterator nit = nodesList.begin(); nit != nodesList.end(); nit++){
		ostringstream iscore;
		iscore << scores.at((*nit)->getId());
		(*nit)->setName("[" + iscore.str() + "]" + ((*nit)->hasName() ? (*nit)->getName() : ""));
	    }
	}
	out << (*currFamily)->getName()<< '\n' <<tpms::TreeTools::nodeToNewick((*currFamily)->getTree()->getRootNode())<< ';' <<endl;
    }
    patienteur.drawFinal();
}

int main(int argc, char *argv[]) {
   try{
       string msg25 = "     version 0.95";
       cout << "      __                                 \n     /\\ \\__   " << msg25 << "\n     \\ \\ ,_\\  _____    ___ ___     ____\n      \\ \\ \\/ /\\ '__`\\/' __` __`\\  /',__\\\n       \\ \\ \\_\\ \\ \\L\\ \\\\ \\/\\ \\/\\ \\/\\__, `\\\n        \\ \\__\\\\ \\ ,__/ \\_\\ \\_\\ \\_\\/\\____/\n         \\/__/ \\ \\ \\/ \\/_/\\/_/\\/_/\\/___/\n                \\ \\_\\\n                 \\/_/    W E L C O M E\n" << endl;
    CmdLineArgs args(argc, argv, "collection,output-dir",cerr);
    args.print(cout);
    
    // Getting Threads Number
    istringstream nbThreadsSS(args.getArg("threads"));
    unsigned int nbThreads;
    nbThreadsSS >> nbThreads;
    if(nbThreads < 1   || nbThreads > 128) nbThreads = 1;
    cout << nbThreads << " threads will be used." << endl;
    
    
    // Opening Collection   
    DataBase collection(args.getArg("collection"),nbThreads);
    
    cout << "To perform these operations, I need to associate each sequence with a species." << endl;
    collection.doFamiliesMapping_LeavesToSpecies();
    
    
    string command("?");
    print_menu();
    while(command != "Q" && !command.empty()){
	
	command = get_choice();
    
	if(command == "CU"){
	    cout << "Scoring: annotate with unicity scores" << endl;
	    collection.doFamiliesMapping_NodesToUnicityScores();
	    
	}else if(command =="CB"){
	    cout << "Scoring: annotate with bipartition scores" << endl;
	    cout << "not available yet" << endl;
	}else if(command == "RU"){
	    cout << "Rooting with the unicity criteria" << endl;
	    collection.doFamiliesRerooting_Unicity();
	}else if(command == "RC"){
            cout << "Rooting with the combo unicity+taxonomy criteria" << endl;
            collection.doFamiliesRerooting_UnicityTaxonomy();
	}else if(command == "RT"){
        cout << "Rooting with the taxonomy criteria" << endl;
        collection.doFamiliesRerooting_Taxonomy();
    }else if(command == "RVD"){
        cout << "Rooting with the daubin’s criteria" << endl;
        collection.doFamiliesRerooting_Daubin();
    }else if(command == "XD"){
	    cout << "Processing LGT detection" << endl;
	    string path = args.getArg("output-dir")+"/detectedTransfers.txt";
	    ofstream resultFile(path.c_str());
	    collection.doFamiliesComputation_detectTransfers(&resultFile);
    }else if(command == "RXD"){
        cout << "Rooting to get less transfers" << endl;
        string path = args.getArg("output-dir")+"/detectedTransfers.txt";
        ofstream resultFile(path.c_str());
        collection.doFamiliesRerooting_LessTransfers(&resultFile);
	}else if(command == "Q" || command.empty()){
	    cout << "Bye" << endl;
	}else if(command == "M"){
	    print_menu();
	}else if(command == "T"){
	    print_status(collection);
	}else if(command == "SU"){
	    string path = args.getArg("output-dir")+"/ANNOTATEDunicityscoresCollection";
	    cout << "Saving to file " << path << " with unicity scores." << endl;
	    save_to_file(collection,path,UNICITY);
	} else if(command == "S"){
	    string path = args.getArg("output-dir")+"/Collection";
	    cout << "Saving to file " << path << endl;
	    save_to_file(collection,path,NONE);
	} else {
	    cout << "The choice " << command << " is not possible." <<endl;
	}
    
    }
    
    }catch(string missingArguments){
	cout << "Missing arguments: " << missingArguments.substr(0,missingArguments.size()-2) << '.' << endl;
	return(1);
    }
}
