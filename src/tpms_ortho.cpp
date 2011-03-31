#include <iostream>
#include <fstream>
#include <string>

#include "classes/DataBase.hpp"
#include "classes/Pattern.hpp"
#include "classes/Query.hpp"
#include "classes/Waiter.hpp"
#include "classes/CmdLineArgs.hpp"

#include <Phyl/Tree.h>

using namespace std;
using namespace bpp;

int main(int argc, char *argv[]) {
  try{
    CmdLineArgs args(argc, argv, "database,testlist,outputdir",cerr);
    args.print(cout);
  
  // DataBase dbTest(argv[1]);
  
  
  string currLine;
  
  ifstream in(args.getArg("testlist").c_str(),ifstream::in);
  
  string path(argv[3]); // contient l'adresse du dossier dans lequel enregistrer les résultats
  //ofstream out(argv[3],ofstream::out);
  
  
  // on lit la première ligne qui contient les valeurs par défaut
  vector<vector<string> > queries;
  vector<string> queryNames;
  
  
  // on commence par extraire chaque ligne
  while(getline(in,currLine)) {
    vector<string> thisLine;
    // currline est la ligne courante, on la découpe selon ESPACE
    string currQueryName;
    
    stringstream currLineSS(currLine);
    
    getline(currLineSS,currQueryName,':');
    queryNames.push_back(currQueryName);
    
    // deuxième partie de la ligne après l'espace : liste des taxons
    string listeTaxons;
    getline(currLineSS,listeTaxons,':');
    istringstream listeTaxonsSS(listeTaxons);
    string currTaxon;
    while(getline(listeTaxonsSS,currTaxon,','))
      thisLine.push_back(currTaxon);
    queries.push_back(thisLine);
  }
  
  DataBase currDB(args.getArg("database"));
  
  cout << "Building loading species trees." << endl;
  currDB.iNeedSpeciesTrees(true,string("cache"),true);
  
  //on suit les queries avec leurs noms : itérateur
  
  vector<string>::iterator currQueryName = queryNames.begin();
  
  for(vector<vector<string> >::iterator oneQuery = queries.begin(); oneQuery < queries.end(); oneQuery++) {
    
      int nombreSetsTrouves = 0;
      
    ofstream out(string(args.getArg("outputdir")+"/"+(*currQueryName)).c_str(),ofstream::out);
    
    cout << "================\nREQUETE  " << *currQueryName << endl;
    
    // pour chaque requete, on doit développer la liste proposée càd passer d'une liste mixe taxons/espèces à une liste espèces
    cout << "Taille de la requete = " << oneQuery->size();
    set<string> species = currDB.getDescendants(*oneQuery,false);
   
    
    //DEBUG: on affiche les espèces de query
    for(set<string>::iterator debugIT = species.begin(); debugIT != species.end(); debugIT++)
	cout << *debugIT << endl;
    
    
    // on va itérer chaque famille disponible dans la banque pour voir si elle contient des orthologues
    vector<Family*> & dbFamilies = currDB.getFamilies();
    Waiter patienteur(&std::cout, dbFamilies.size(), '#');
    for(vector<Family*>::iterator currFamily = dbFamilies.begin(); currFamily != dbFamilies.end(); currFamily++){
      patienteur.step();
      // première étape, on vérifie que la famille contient bien tous les taxons demandés
      if((*currFamily)->containsSpecies(species)){
	// si c'est bon, on explore les nœuds un à un
	// on charge la liste des Feuilles desquelles partir
	vector<Node *> leaves = (*currFamily)->getSpTree()->getLeaves();
	
	// on initialise la liste des feuilles interdites
	set<Node  *> forbiddenLeaves;
	
	// on compte le nombre de set d'orthologues trouvés dans chaque famille
	int compteurTrouves = 0;
	
	// on itère les feuilles
	for(vector<Node *>::iterator leave = leaves.begin(); leave < leaves.end(); leave++){
	  // on vérifie que la feuille n'est pas dans les interdites, et qu'elle appartient à la liste des espèces
	  if(forbiddenLeaves.find(*leave) == forbiddenLeaves.end() && species.find((*leave)->getName()) != species.end()){ // == si elle n'est pas interdite et qu'elle est dans species
	    // à partir de maintenant, on remonte jusqu'en haut de la monophylie
	    Node * currNode = (*leave);
	    while(currNode->hasFather() && Pattern::nodeOnlyContainsTheseTaxa(currNode->getFather(),species,false))
	      currNode = currNode->getFather();
	    // maintenant, on est en haut de la monophylie
	    // quelles sont les feuilles de ce sous-arbre ?
	    set<Node *> currLeaves;
	    int currLeavesNumber = 0;
	    Family::getLeavesFromNode(currNode,currLeaves,currLeavesNumber);
	    set<string> currLeavesNames = Family::nodesToNames(currLeaves);
	    //reste juste à vérifier si le sous-arbre contient tous les représentants de species
	    if(currLeavesNumber == species.size() &&  currLeavesNames == species) {
		compteurTrouves++;
		nombreSetsTrouves++;
		// TROUVE !!! si on est ici, c'est qu'on a une famille d'orthologues telle que recherchée, on doit stocker le résultat dans le fichier
		// il s'agit ici d'obtenir le nœud racine de notre orthologie, mais sur l'arbre des séquences.
		Node * racineSsArbreGenes = (*currFamily)->getTree()->getNode(currNode->getId());
		set<Node *> geneTreeLeaves;
		Family::getLeavesFromNode(racineSsArbreGenes,geneTreeLeaves);
		set<string> listeGenes = Family::nodesToNames(geneTreeLeaves);
		
		string nomFamilleComplet = (*currFamily)->getName();
		
		// à partir d'ici, on a tout ce qu'il faut pour écrire dans le fichier de résultats
		for(set<string>::iterator currMnemo = listeGenes.begin(); currMnemo != listeGenes.end(); currMnemo++){
		    out << nomFamilleComplet << " " << compteurTrouves << " " << *currMnemo << endl;
		}
	    }
	    // on interdit de regarder à nouveau dans les nœuds déjà explorés
	    forbiddenLeaves.insert(currLeaves.begin(),currLeaves.end());
	  }
	}
	
	
      
      }

    }
    
    patienteur.drawFinal();
    cout << "... Trouve : " << nombreSetsTrouves << " sets." << endl;
    out.close();
    currQueryName++;
  }
  
  } catch(string  missingArguments){
	cout << "Missing arguments: " << missingArguments.substr(0,missingArguments.size()-2) << '.' << endl;
	return(1);
    }
    return(0);

  
}	

