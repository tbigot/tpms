#include <iostream>
#include <fstream>
#include <string>

#include "classes/DataBase.hpp"
#include "classes/Pattern.hpp"
#include "classes/Query.hpp"

using namespace std;

Query * stGet(DataBase * db);
std::string queryTaxon(std::string message,DataBase * database);
std::string taxonCheck(std::string qTaxon,DataBase * database);



int main(int argc, char *argv[]) {
	DataBase currDB(argv[1]);
	ofstream lastResults("tpms_LastResults",ofstream::out);
	vector<int> xFamilies;
	
	Query * currQuery = stGet(&currDB);
	
	// boucle d'interaction
	while(currQuery != 00){
		vector<int> familles;
		try{
		string monophilyTaxon = currDB.getParentTaxon(currQuery->getTarget(),currQuery->getMonophilyLevel());
		Pattern currPattern(currQuery->getPatternString(),&currDB);
		cout << "Recherche " << currQuery->getTarget() << " chez les NON-" << currQuery->getSource() << "." << endl;
		xFamilies = currPattern.search(currDB.getFamilies());
		// xFamilies.empty();
		// xFamilies.push_back(421474);
		cout << "Terminé. " << xFamilies.size() << " familles trouvées.\n  ----" << endl;
		cout << "Deuxième passe : contrôle des nœuds parents selon paramétrage. Monophylie demandée : " << monophilyTaxon << endl;
		xFamilies = currPattern.xferDetected(currDB.getFamilies(),
									   xFamilies,
									   "!"+currQuery->getSource(), //source
									   currQuery->getTarget(), // target
									   monophilyTaxon,
									   // currDB.getParentTaxon((*oneQuery)->getTarget(),0), // monophily
									   currQuery->getDepth(),
									   currQuery->getBootstraps(),
									   currQuery->getSourceRates());
		
		int total = currDB.nbFamiliesContaining(currQuery->getTarget());							   
									   
		cout << "Terminé. " << xFamilies.size() << " familles trouvées sur "<< total << ".\n  ----" << endl;
									   

		stringstream tempFamNum;
		string famNum;
		
		cout << '\n';
		
		for(vector<int>::iterator it = xFamilies.begin(); it != xFamilies.end(); it++) {
			tempFamNum.str("");
			tempFamNum << *it;
			famNum = tempFamNum.str();
			while(famNum.size() < 6) {famNum = "0" + famNum;}
			lastResults << "HBG" << famNum << endl;
			cout << "HBG" << famNum << " :: ";
		}
		cout << endl;
			
		} catch (bpp::NodeNotFoundException e) {
			cout << "Noeud non trouvé : \n" << e.what() << endl;
		}
		currQuery = stGet(&currDB);
		
		
	}
}

Query * stGet(DataBase * db) {
	string source = queryTaxon("Entrez le taxon source (vide pour quitter) : ",db);
	if(source.empty()) return(00);
	
	string target = queryTaxon("Entrez le taxon cible : ",db);
	
	cout << "Entrez le niveau de tolérance à la monophylie : " << flush;
	string sMonophyly;
	getline(cin,sMonophyly);
	
	cout << "Entrez le nombre de noeuds a remonter : " << flush;
	string sCtrlLevel;
	getline(cin,sCtrlLevel);
	
	cout << "Entrez la liste des paramètres en %  séparés par des espaces (bs1,taux1 bs2,taux2) : " << flush;
	string sParams;
	getline(cin,sParams);
	
	return(new Query(source + ":" + target + ":" + sMonophyly + ":" + sCtrlLevel + ":" + sParams,""));
	
}

string queryTaxon(string message,DataBase * db) {
	string result, candidates;
	bool taxonOK = false;
	while(!taxonOK){
		cout << message << flush;
		getline(cin,result);
		if(result.empty()) return("");
		std::transform(result.begin(), result.end(), result.begin(), (int(*)(int)) toupper);
		candidates = taxonCheck(result,db);
		if(candidates.empty()) taxonOK = true; else cout << candidates << endl;
	}
	return(result);
}

string taxonCheck(string qTaxon,DataBase * db) {
	string result;
	try{
		db->getSpeciesTree()->getNode(qTaxon);
	} catch (bpp::NodeNotFoundException e) {
		result += " ";
		std::vector<bpp::Node *> nodes = db->getSpeciesTree()->getNodes();
		for(std::vector<bpp::Node *>::iterator oneNode = nodes.begin(); oneNode < nodes.end(); oneNode++) {
			if((*oneNode)->hasName() && (*oneNode)->getName().substr(0,qTaxon.size()) == qTaxon)
				result = result + (*oneNode)->getName() + " :: ";
		}
	}
	return(result);
}
