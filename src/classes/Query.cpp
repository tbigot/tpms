#include "Query.hpp"

#include <sstream>
#include <iostream>


using namespace std;

namespace tpms{
Query::Query(string queryLine, string defaultsParams) {
	// à partir d'ici on décompose la ligne queryLine
	// rappel : elle est du format :
	// gammaproteobacteria:ESCHERICHIA COLI K12 K12:2:2:90,70 90,50
	// source:cible:monophilyLevel:nombre de noeuds à remonter:param noeud1,param noeud2
	
	
	stringstream ssql;
	
	ssql << queryLine;
	
	if(queryLine.at(queryLine.size()-1) == ':') ssql << defaultsParams;
	
	// cout << ">> "  << ssql.str() << endl;
	
	getline(ssql,sourceTaxon,':');
	std::transform(sourceTaxon.begin(), sourceTaxon.end(), sourceTaxon.begin(), (int(*)(int)) toupper); // on passe en majuscules
	getline(ssql,targetTaxon,':');
	std::transform(targetTaxon.begin(), targetTaxon.end(), targetTaxon.begin(), (int(*)(int)) toupper); // on passe en majuscules
	
	string monophilyLevelS;
	getline(ssql,monophilyLevelS,':');
	istringstream monophilyLevelSS(monophilyLevelS);
	monophilyLevelSS >> monophilyLevel;
		
	string verifDepthS;
	getline(ssql,verifDepthS,':');
	istringstream verifDepthSS(verifDepthS);
	verifDepthSS >> verifDepth;

	string params;
	getline(ssql,params,':');
	
	
	// remplissage des vecteurs de paramètres
	istringstream paramsSS(params);
	string oneNodeParam;
	string oneBS,oneSR;
	int oneBSi,oneSRi;
	while(getline(paramsSS,oneNodeParam,' ')) {
		istringstream oneNodeParamSS(oneNodeParam);
		getline(oneNodeParamSS,oneBS,',');
		getline(oneNodeParamSS,oneSR,',');
		istringstream oneBSss(oneBS);
		istringstream oneSRss(oneSR);
		oneBSss >> oneBSi;
		oneSRss >> oneSRi;
		bootstraps.push_back(oneBSi);
		sourceRates.push_back(oneSRi);
	}
	
	// affichage debug
	// cout << "Query chargée avec " << sourceTaxon << " ::: " << targetTaxon << " ::: " << verifDepth << " ::: " << params << endl;
	// cout << "Bootstraps : ";
	//for(vector<unsigned int>::iterator bsp = bootstraps.begin(); bsp < bootstraps.end(); bsp++) cout << *bsp << " ::: ";
	// cout << "\nSourceRates : ";
	//for(vector<unsigned int>::iterator bsp = sourceRates.begin(); bsp < sourceRates.end(); bsp++) cout << *bsp << " ::: ";
	//cout << endl;
	
	
	
}

string Query::getPatternString() {
	// return("(" + sourceTaxon + ",(!" + sourceTaxon + ",(!" + sourceTaxon + "," + targetTaxon + ")));");
	return("(!" + sourceTaxon + ",(!" + sourceTaxon + "," + targetTaxon + "));");
}

string Query::getTarget() { return(targetTaxon); }

string Query::getSource() { return(sourceTaxon); }

vector<unsigned int> Query::getBootstraps() { return(bootstraps); }

vector<unsigned int> Query::getSourceRates() { return(sourceRates); }

unsigned int Query::getMonophilyLevel() { return(monophilyLevel); }

unsigned int Query::getDepth() { return(verifDepth); }

bool Query::check(DataBase * db){
	if(!db->taxonExists(sourceTaxon)) {
		cerr << "###ERR : Taxon source " << sourceTaxon << " inconnu" << endl;
		return(false);
	}
	if(!db->taxonExists(targetTaxon)) {
		cerr << "###ERR : Taxon cible " << targetTaxon << " inconnu" << endl;
		return(false);
	}
	
	return(true);
}
}