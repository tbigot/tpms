//
// File: Query.cpp
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
	
	getline(ssql,sourceTaxon_,':');
	std::transform(sourceTaxon_.begin(), sourceTaxon_.end(), sourceTaxon_.begin(), (int(*)(int)) toupper); // on passe en majuscules
	getline(ssql,targetTaxon_,':');
	std::transform(targetTaxon_.begin(), targetTaxon_.end(), targetTaxon_.begin(), (int(*)(int)) toupper); // on passe en majuscules
	
	string monophilyLevelS;
	getline(ssql,monophilyLevelS,':');
	istringstream monophilyLevelSS(monophilyLevelS);
	monophilyLevelSS >> monophilyLevel_;
		
	string verifDepthS;
	getline(ssql,verifDepthS,':');
	istringstream verifDepthSS(verifDepthS);
	verifDepthSS >> verifDepth_;

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
		bootstraps_.push_back(oneBSi);
		sourceRates_.push_back(oneSRi);
	}
	
	// affichage debug
	// cout << "Query chargée avec " << sourceTaxon << " ::: " << targetTaxon << " ::: " << verifDepth << " ::: " << params << endl;
	// cout << "Bootstraps : ";
	//for(vector<unsigned int>::iterator bsp = bootstraps.begin(); bsp < bootstraps.end(); bsp++) cout << *bsp << " ::: ";
	// cout << "\nSourceRates : ";
	//for(vector<unsigned int>::iterator bsp = sourceRates.begin(); bsp < sourceRates.end(); bsp++) cout << *bsp << " ::: ";
	//cout << endl;
	
	
	
}



string Query::getTarget() { return(targetTaxon_); }

string Query::getSource() { return(sourceTaxon_); }

vector<unsigned int> Query::getBootstraps() { return(bootstraps_); }

vector<unsigned int> Query::getSourceRates() { return(sourceRates_); }

unsigned int Query::getMonophilyLevel() { return(monophilyLevel_); }

unsigned int Query::getDepth() { return(verifDepth_); }


}