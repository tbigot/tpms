#include <iostream>
#include <fstream>
#include <string>

#include "classes/DataBase.hpp"
#include "classes/Pattern.hpp"

using namespace std;

void stGet(string * source, string * target, int * ctrlLevel, int * bootstrap,DataBase * database);
std::string queryTaxon(std::string message,DataBase * database);
std::string taxonCheck(std::string qTaxon,DataBase * database);

int main(int argc, char *argv[]) {
	DataBase dbTest(argv[1]);
		
	string taxonSource;
	string taxonTarget;
	int ctrlLevel;
	int bootstrap;
	
	string sPattern;
	
	stGet(&taxonSource,&taxonTarget,&ctrlLevel,&bootstrap,&dbTest);
	
	// boucle d'interaction
	while(taxonSource != ""){
		vector<int> familles;
		try{
			// generation du pattern
			sPattern = "(" + taxonSource + ",(!" + taxonSource + ",(!" + taxonSource + "," + taxonTarget + ")));";
			cout << "\nPattern généré : " << sPattern << endl;
			Pattern tp(sPattern,dbTest.getSpeciesTree());
			cout << "Pattern créé. On recherche maintenant le motif " << sPattern << " dans la base." << endl;
			
			familles = tp.search(dbTest.getFamilies());
			
			
			//debug: on les prend tous
			/*map<int,Family *> * toutesFamilles = dbTest.getFamilies();
			for(map<int,Family *>::iterator laFamNr = toutesFamilles->begin(); laFamNr != toutesFamilles->end(); laFamNr++) {
				familles.push_back(laFamNr->first);
			}*/
			
			cout << "Deuxième étape : on cherche tous les « " << taxonTarget << " » au milieu des « non " << taxonSource << " » avec " << ctrlLevel << " niveaux de vérification, et un minimum de bootstrap de " << bootstrap << "%." << endl;
			// on liste les réponses
			cout << '[';
			for(vector<int>::iterator it = familles.begin(); it != familles.end(); it++) {
				cout << "HBG" << *it << ',';
			}
			cout << ']' << endl;
			
			// on met les réponses dans un fichier lastResults
			ofstream out("./lastResults", ofstream::out);
			stringstream tempFamNum;
			string famNum;
			for(vector<int>::iterator it = familles.begin(); it != familles.end(); it++) {
				tempFamNum.str("");
				tempFamNum << *it;
				famNum = tempFamNum.str();
				while(famNum.size() < 6) {famNum = "0" + famNum;}
				out << "HBG" << famNum << endl;
			}
			
			
			cout << "TOTAL: " << familles.size() << " / " << dbTest.getNbFamilies() << endl;
			
		} catch (bpp::NodeNotFoundException e) {
			cout << "Noeud non trouvé : \n" << e.what() << endl;
		}
		
		stGet(&taxonSource,&taxonTarget,&ctrlLevel,&bootstrap,&dbTest);
		
	}
}

void stGet(string * source, string * target, int * ctrlLevel, int * bootstrap,DataBase * db) {
	*source = queryTaxon("Entrez le taxon source : ",db);
	*target = queryTaxon("Entrez le taxon cible : ",db);
	
	cout << "Entrez le niveau de vérification : " << flush;
	string sCtrlLevel;
	getline(cin,sCtrlLevel);
	stringstream ssCtrlLevel(sCtrlLevel);
	ssCtrlLevel >> *ctrlLevel;
	cout << "Entrez le bootstrap minimum : " << flush;
	string sBootstrap;
	getline(cin,sBootstrap);
	stringstream ssBootstrap(sBootstrap);
	ssBootstrap >> *bootstrap;
	
}

string queryTaxon(string message,DataBase * db) {
	string result, candidates;
	bool taxonOK = false;
	while(!taxonOK){
		cout << message << flush;
		getline(cin,result);
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
