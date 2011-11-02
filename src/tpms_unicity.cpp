#include <iostream>
#include <fstream>
#include <string>

#include <Bpp/Phyl/Tree.h>

#include "classes/Family.hpp"
#include "classes/DataBase.hpp"
#include "classes/Pattern.hpp"
#include "classes/Query.hpp"
#include "classes/Waiter.hpp"
#include "classes/CmdLineArgs.hpp"
#include "classes/TreeTools.hpp"


using namespace std;
using namespace bpp;
using namespace tpms;

int main(int argc, char *argv[]) {
    try{
	CmdLineArgs args(argc, argv, "database,output",cerr);
	
	args.print(cout);
	
	DataBase currDB(args.getArg("database"));
	ofstream out(args.getArg("output").c_str(),ofstream::out);
	
	cout << "Building loading species trees." << endl;
	currDB.doFamiliesMapping_LeavesToSpecies();
	
	cout << "Generating UnicityScores (may take a while, please have a coffe - ask Jos)." << endl;
	currDB.doFamiliesMapping_NodesToBestUnicityScores();
	
	// liste des familles
	vector<Family*> & dbFamilies = currDB.getFamilies();
	
	cout << "Yew! I've finished the work, now I can write resulting trees into newick files" << endl;
	
	Waiter patienteur(&std::cout, dbFamilies.size(), '#');
	for(vector<Family*>::iterator currFamily = dbFamilies.begin(); currFamily != dbFamilies.end(); currFamily++){
	    patienteur.step();
	    
	    vector<unsigned int> scores = (*currFamily)->getUnicityScores();
	    
	    // labelling the node names
	    
	    vector<Node *> nodesList = (*currFamily)->getTree()->getNodes();
	    for(vector<Node *>::iterator nit = nodesList.begin(); nit != nodesList.end(); nit++){
		ostringstream iscore;
		iscore << scores.at((*nit)->getId());
		(*nit)->setName("[" + iscore.str() + "]" + ((*nit)->hasName() ? (*nit)->getName() : ""));
	    }
	    
	    out << (*currFamily)->getName()<< '\n' <<tpms::TreeTools::nodeToNewick((*currFamily)->getTree()->getRootNode())<< ';' <<endl;
	}
	patienteur.drawFinal();
	
	
	cout << "Terminay, bye." << endl;
	
    } catch(string  missingArguments){
	cout << "Missing arguments: " << missingArguments.substr(0,missingArguments.size()-2) << '.' << endl;
	return(1);
    }
    return(0);
    
    
}	

