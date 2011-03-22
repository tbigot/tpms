#include <iostream>
#include <fstream>
#include <string>

#include "classes/DataBase.hpp"
#include "classes/Pattern.hpp"
#include "classes/Query.hpp"


using namespace std;

void stGet(string * source, string * target, int * ctrlLevel, int * bootstrap,DataBase * database);
std::string queryTaxon(std::string message,DataBase * database);
std::string taxonCheck(std::string qTaxon,DataBase * database);

int main(int argc, char *argv[]) {
	
	if(argc != 4){
		cout << "Usage : " << argv[0] << " dbFile testList output" << endl;
		exit(1);
	}
	
	// DataBase dbTest(argv[1]);

	
	string currLine;
	
	ifstream in(argv[2],ifstream::in);
	ofstream out(argv[3],ofstream::out);
	ofstream lastResults("lastResults",ofstream::out);
	
	
	// on lit la première ligne qui contient les valeurs par défaut
	string defaultsParams;
	getline(in,defaultsParams);
	
	vector<Query*> queries;
	
	while(getline(in,currLine)) {
		Query * cQ = new Query(currLine,defaultsParams);
		queries.push_back(cQ);
	}
	
	vector<int> xFamilies;
	
	DataBase currDB(argv[1]);
	
	
	//ligne d'en-tête pour le fichier de sortie csv
	out << "source,target,monophily,rate1,rate2,detected,total" <<endl;
	
	double trates = 0;
	
	for(vector<Query*>::iterator oneQuery = queries.begin(); oneQuery < queries.end(); oneQuery++) {
		
		
		
		string monophilyTaxon = currDB.getParentTaxon((*oneQuery)->getTarget(),(*oneQuery)->getMonophilyLevel());
		Pattern currPattern((*oneQuery)->getPatternString(),currDB.getSpeciesTree());
		cout << "Recherche " << (*oneQuery)->getTarget() << " chez les NON-" << (*oneQuery)->getSource() << "." << endl;
		xFamilies = currPattern.search(currDB.getFamilies());
		// xFamilies.empty();
		// xFamilies.push_back(421474);
		cout << "Terminé. " << xFamilies.size() << " familles trouvées.\n  ----" << endl;
		cout << "Deuxième passe : contrôle des nœuds parents selon paramétrage. Monophylie demandée : " << monophilyTaxon << endl;
		xFamilies = currPattern.xferDetected(currDB.getFamilies(),
									   xFamilies,
									   "!"+(*oneQuery)->getSource(), //source
									   (*oneQuery)->getTarget(), // target
									   monophilyTaxon,
									   // currDB.getParentTaxon((*oneQuery)->getTarget(),0), // monophily
									   (*oneQuery)->getDepth(),
									   (*oneQuery)->getBootstraps(),
									   (*oneQuery)->getSourceRates());
		
		int total = currDB.nbFamiliesContaining((*oneQuery)->getTarget());							   
									   
		cout << "Terminé. " << xFamilies.size() << " familles trouvées sur "<< total << ".\n  ----" << endl;
									   
		// Résultats pour cette requete
		out << (*oneQuery)->getSource() << "," << (*oneQuery)->getTarget() << "," << monophilyTaxon << "," << (*oneQuery)->getSourceRates().at(0) << "," << (*oneQuery)->getSourceRates().at(1)<< "," << xFamilies.size() << ',' << total << endl;
		
		
		
		trates += (double)xFamilies.size()/ (double)total;
		
		stringstream tempFamNum;
		string famNum;
		for(vector<int>::iterator it = xFamilies.begin(); it != xFamilies.end(); it++) {
			tempFamNum.str("");
			tempFamNum << *it;
			famNum = tempFamNum.str();
			while(famNum.size() < 6) {famNum = "0" + famNum;}
			lastResults << "HBG" << famNum << endl;
		}
	}
	
	
	
	cout << "#### #### " << trates / (double)queries.size() << endl;
	cout << defaultsParams << endl;
	
}	

