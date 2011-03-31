#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>

#include "Pattern.hpp"
#include "NodeConstraints.hpp"
#include "TreeTools.hpp"
#include "Waiter.hpp"
#include "CandidateNode.hpp"


using namespace std;
using namespace bpp;
using namespace tpms;

Pattern::Pattern(TreeTemplate<Node> &tree, DataBase &refDB): refDB(refDB), tree(tree){
    
    extractConstraints();
    developSpecies();
    cout << "\n -- Constraints visual representation:" << endl;
    print(tree.getRootNode(),0);
}

Pattern::~Pattern(){
    for(vector<NodeConstraints *>::iterator cc=constraints.begin(); cc != constraints.end(); cc++)
	delete(*cc);
}


void Pattern::extractConstraints(){
    // les noms de feuilles peuvent contenir des contraintes entre accolades
    // la contrainte est portée par le nœud fils alors qu'il s'agit d'une contrainte sur la branche menant au nœud fils
    // et ceci pour des raisons de commodité.
    
    // resize et remplissage : ici, on sait ce qu'on fait.
    constraints.resize(tree.getNumberOfNodes());
    for(unsigned int i = 0; i!= constraints.size();i++)
	constraints.at(i) = new NodeConstraints(refDB);
    
    vector<Node *> noeuds = tree.getNodes();
    for(vector<Node *>::iterator unNoeud = noeuds.begin(); unNoeud != noeuds.end(); unNoeud++){
	string initName;
	if((*unNoeud)->hasName()) initName = (*unNoeud)->getName();
	string nodeName = initName;
	// si le nom du nœud contient une accolade ouvrante, on définit la contrainte
	string::size_type found = initName.find('{');
	if (found != string::npos){ // = trouvé
	    if(!(*unNoeud)->hasFather())
		cout << "\n!!!!!! Constraints on root are NOT YET supported !" << endl ;
	    nodeName = initName.substr(0,found);
	    (constraints.at((*unNoeud)->getId()))->setConstraints(refDB,initName.substr(found+1, (initName.size()-found-2)));
	    if(!nodeName.empty()) (*unNoeud)->setName(nodeName); else (*unNoeud)->deleteName();
	}
	
    }
}


void Pattern::developSpecies() {
    // cette fonction associe, pour chaque feuille du pattern, une liste d'espèces
    
    // quand on va entamer une recherche, et qu'il y a un taxon dans le pattern, il faut savoir exactement quelles espèces de l'arbre il recouvre. Si un nœud de l'abre est « BACTERIA », on doit avoir la liste de toutes les espèces correspondantes, c'est-à-dire toutes les feulles qui sont sous BACTERIA dans l'abre des espèces.
    
    // species est donc une map qui associe à chaque ID de nœud un set d'espèces (liste ordonnée) sous forme de chaînes.
    // on récupère tous les nœuds, mais seules les feuilles nous intéressent.
    vector<Node *> feuilles = tree.getNodes();
    
    string currNodeName;
    
    set<string> allSons;
    for(vector<Node *>::iterator it = feuilles.begin(); it < feuilles.end(); it++) {
	if((*it)->hasName()) {
	    currNodeName = (*it)->getName();
	    if(currNodeName.at(0) == '!') currNodeName = currNodeName.substr(1,currNodeName.size()-1);
	    allSons = refDB.getDescendants(currNodeName);
	    species.insert(pair<int,set<string> >((*it)->getId(),allSons));
	    
	    /*TO KNOW THE SPECIES LIST THAT IS ACCEPTED BY A NODE
	     * 
	     * 
	     * for(set<string>::iterator lspc=allSons.begin(); lspc != allSons.end(); lspc++){
		cout << *lspc << ", ";
	    }
	    cout << endl;*/
	    
	}
    }
    
}


unsigned int Pattern::search(std::vector< Family* >& families, vector< pair< Family*, CandidateNode* > >& result){
    unsigned int found = 0;
    Waiter patienteur(&cout, families.size(), '#');	
    
    // cette fonction parcourt la liste des familles, et pour chaque famille, elle doit vérifier qu'au moins une espèce de chaque noeud du pattern est présent dans la liste d'espèce de la famille.
    
    for(vector<Family *>::iterator theFamily = families.begin(); theFamily != families.end(); theFamily++) {

	//premiere étape = vérification des taxons
	// -- on liste les noeuds du pattern : les especes sont contenues dans species
	bool tousNoeudsOntEspece = true; // tous les noeuds correspodent-ils à au moins une espèce de la famille ?
	// si un seul noeud n'a pas d'équivalence, la famille ne nous intéresse pas. d'où le « oneNodeList » dans les conditions du for
	for(map<int,set<string> >::iterator oneNodeList = species.begin(); oneNodeList != species.end() && tousNoeudsOntEspece ; oneNodeList++) {
	    // la valeur de chaque élément de la map est un vecteur de strings correspondant aux espèces
	    // -- on est dans un noeud du pattern avec une liste d'espèces candidates. On doit trouver si au moins une espèce est présente dans la famille.
	    bool uneEspeceTrouvee = false; // y a-t-il au moins une espèce en commun ?
	    // de deux choses l'une : soit on est dans le cas d'un nœud classique avec un nom de taxon, dans ce cas, on doit trouver au moins un représentant de ce taxon dans la famille
	    // soit c'est un nœud négatif, qui commence par un «!», auquel cas, il faudrait trouver au moins une espèce qui n'appartient pas au taxon donné. On n'a aucun vecteur qui contienne ce type de liste d'espèces. Souvent le complémentaire d'un taxon est vaste, donc on part sur le principe qu'il y a un représentant dans la liste. Au pire on aura fait du pattern matching pour rien, mais peu probable
	    
	    if(tree.getNode(oneNodeList->first)->getName().at(0) == '!') uneEspeceTrouvee = true;
	    
	    for(set<string>::iterator thePatternSpecie = oneNodeList->second.begin(); thePatternSpecie != oneNodeList->second.end() && !uneEspeceTrouvee; thePatternSpecie++) {
		
		
		if((*theFamily)->getSpecies()->find(*thePatternSpecie) != (*theFamily)->getSpecies()->end())
		    uneEspeceTrouvee = true;
	    }
	    tousNoeudsOntEspece = tousNoeudsOntEspece && uneEspeceTrouvee;
	    // si le nœud qui vient d'être fait n'a pas d'espèce, tousNoeudsOntEspece devient faux
	    // et dans ce cas, on arrête la boucle
	}
	// en arrivant ici, soit on a passé tous les test avec succès tousNoeudsOntEspece=vrai, soit on s'est arrêté en route parce qu'un nœud n'avait pas d'espèce dans la famille
	
	
	if(tousNoeudsOntEspece) {
	    CandidateNode * candRoot = new CandidateNode();
	    if(patternMatch((*theFamily)->getSpTree()->getRootNode(),tree.getRootNode(),candRoot)){
		result.push_back(pair<Family *, CandidateNode*>(*theFamily,candRoot));
		found++;
	    }
	    else
		delete(candRoot);
	    
	}
		
	
	
	patienteur.step();
    }
    return(found);
    
}

bool Pattern::patternMatch(Node * target, Node * pattern, CandidateNode * fatherCandidate) {
    // On reprend ici l'algoritme de JF Dufayard (Dufayard et al, 2005) duquel on a retiré la problématique de contraintes de branches
    
    
    if(isLeaf(target) && isLeaf(pattern)){
	if(
	    (pattern->getName().at(0)!='!' && isInTaxa(pattern,target->getName()))
	    || (pattern->getName().at(0)=='!' && !isInTaxa(pattern,target->getName()))
	)
	{
	    CandidateNode * currCandidate = new CandidateNode(fatherCandidate,target,pattern);
	    currCandidate->confirm();
	    return(true);
	} else return(false);
    } else if(isLeaf(target) && !isLeaf(pattern)) {
	return(false);
    } else if(!isLeaf(target) && isLeaf(pattern)) {
	bool result = patternMatch(target->getSon(0),pattern,fatherCandidate);
	result = patternMatch(target->getSon(1),pattern,fatherCandidate) || result;
	return(result);
    } else {
	    Node * tson1 = target->getSon(0);
	    Node * tson2 = target->getSon(1);
	    Node * pson1 = pattern->getSon(0);
	    Node * pson2 = pattern->getSon(1);
	    
	    NodeConstraints::NodeNature targetNature = NodeConstraints::SPECIATION;
	    if(target->hasName() && target->getName().at(0) == '#') targetNature = NodeConstraints::DUPLICATION;
	    
	    CandidateNode * candidate = new CandidateNode(fatherCandidate, target, pattern);
	    
	 if(
	    (constraints.at(pattern->getId())->isAuthorizedNature(targetNature) ) &&
	    (!constraints.at(pson1->getId())->hasSpeciesRestrictions()
	    || (nodeOnlyContainsTheseTaxa(tson1,constraints.at(pson1->getId())->getAuthorisedSpecies()) )
	    )
	    
	    &&
	    (!constraints.at(pson2->getId())->hasSpeciesRestrictions()
	    || (nodeOnlyContainsTheseTaxa(tson2,constraints.at(pson2->getId())->getAuthorisedSpecies())   )
	    )
	    
	    && ( patternMatch(tson1, pson1, candidate)
	    && patternMatch(tson2, pson2, candidate))
	    ) {
	     candidate->confirm(); 
	     return(true);}
	 else delete(candidate);
	    
	 candidate = new CandidateNode(fatherCandidate, target, pattern);
	 if	(
	    (!constraints.at(pson1->getId())->hasSpeciesRestrictions()
	    || (nodeOnlyContainsTheseTaxa(tson2,constraints.at(pson1->getId())->getAuthorisedSpecies()) )
	    )
	    
	    &&
	    (!constraints.at(pson2->getId())->hasSpeciesRestrictions()
	    || (nodeOnlyContainsTheseTaxa(tson1,constraints.at(pson2->getId())->getAuthorisedSpecies()) )
	    )
	    
	    && ( patternMatch(tson2, pson1, candidate)
	    && patternMatch(tson1, pson2, candidate))
	) {
	     candidate->confirm();
	     return(true);}
	else delete(candidate);
	
	bool result = patternMatch(tson1, pattern, fatherCandidate);
	result = patternMatch(tson2, pattern, fatherCandidate) || result;
	return(result);
	
	
	
    }
    
}

// extrait de tous les nœuds de l'arbre ceux qui appartiennent à un set de taxons
std::vector<int> Pattern::getIdWithTaxaList(bpp::TreeTemplate<bpp::Node> * sTree, std::set<string> * taxa) {
    // trouve tous les IDs des noeuds qui appartiennent à la liste taxons
    vector<int> retour;
    Node * currNode;
    
    vector<int> allIds = sTree->getNodesId();
    for(vector<int>::iterator currNodeId = allIds.begin(); currNodeId != allIds.end(); currNodeId++) {
	currNode = sTree->getNode(*currNodeId);
	if(currNode->hasName()) if(taxa->find(currNode->getName()) != taxa->end())
	    retour.push_back(*currNodeId);
    }
    return(retour);
}

bool Pattern::isLeaf(Node * pNode) { return(pNode->getNumberOfSons() == 0); }

// ce nœud autorise-t-il cette espèce ?
bool Pattern::isInTaxa(Node * pNode,string pSpecie){
    // species associe à chaque ID de noeud un set d'espèces
    map<int, set<string> >::iterator it = species.find(pNode->getId());
    // it pointe maintenant sur le vector qui nous intéresse. On cherche l'espèce dedans.
    set<string>::iterator it2 = it->second.find(pSpecie);
    return(it2 != it->second.end());
}


bool Pattern::nodeOnlyContainsTheseTaxa(Node * localRoot, set<string> & taxonMembers, bool invert){
    // monophyletique 
    //FONCTION RÉCURSIVE
    
    // - cas de base : feuille.
    // Soit elle n'appartient pas à tsMembers, soit c'est le nœud targetNode, sinon false.
    if(isLeaf(localRoot)){
	if(invert){
	    if(taxonMembers.find(localRoot->getName()) == taxonMembers.end()) return(true); else return(false);
	} else {
	    if(taxonMembers.find(localRoot->getName()) != taxonMembers.end()) return(true); else return(false);
	}
    } else { // -cas récursif = un noeud
	bool tousLesDescendantsOK = true; // par défaut si pas de fils, ils sont tous descendants du taxon. On opérera avec des ET
	for(unsigned int currFilsIndex = 0; tousLesDescendantsOK && currFilsIndex < localRoot->getNumberOfSons(); currFilsIndex++){
	    tousLesDescendantsOK = tousLesDescendantsOK && nodeOnlyContainsTheseTaxa(localRoot->getSon(currFilsIndex), taxonMembers, invert);
	}
	return(tousLesDescendantsOK);
}
}

int Pattern::enumerateTaxon(Node * localRoot, set<string> * tsMembers, bool tsComplementary, Node * impasse, unsigned int * totsize){
    //renvoie le taux d'espèces membres du taxon source (tsMembers)
    //FONCTION RÉCURSIVE
    
    // - cas de base : feuille.
    // Soit elle n'appartient pas à tsMembers, soit c'est le nœud targetNode, sinon false.
    if(isLeaf(localRoot)){ // cas de base numéro 1: une feuille
	*totsize = 1;
	if(tsComplementary){
	    if(tsMembers->find(localRoot->getName()) == tsMembers->end()) return(1); else return(0);
	} else {
	    if(tsMembers->find(localRoot->getName()) != tsMembers->end()) return(1); else return(0);
	}
} else if(localRoot == impasse) { // cas de base numéro 2: l'impasse
    *totsize = 0;
    return(0);
}else { // -cas récursif = un noeud
    unsigned int localTotSize = 0;
    unsigned int localNumber = 0;
    unsigned int currSonTotSize = 0;
    
    for(unsigned int currFilsIndex = 0; currFilsIndex < localRoot->getNumberOfSons(); currFilsIndex++){
	localNumber += enumerateTaxon(localRoot->getSon(currFilsIndex), tsMembers, tsComplementary, impasse, &currSonTotSize);
	localTotSize += currSonTotSize;
    }
    *totsize = localTotSize;
    return(localNumber);
}
}


vector<int> Pattern::xferDetected(map<int,Family *> * families, vector<int> selector, string sourceTaxon, string targetTaxon, string monophylyTaxon, unsigned int verifDeep, vector<unsigned int> bootstraps, vector<unsigned int> sourceRates) {
    // cette fonction doit trouver dans la liste famList toutes les familles qui groupent un membre de taxonTarget dans taxonSource
    bool tsComplementary = sourceTaxon.at(0) == '!'; if(tsComplementary) sourceTaxon = sourceTaxon.substr(1,sourceTaxon.size()-1);
    set<string> tsMembers = refDB.getDescendants(sourceTaxon);
    bool ttComplementary = targetTaxon.at(0) == '!'; if(ttComplementary) targetTaxon = targetTaxon.substr(1,targetTaxon.size()-1);
    set<string> ttMembers = refDB.getDescendants(targetTaxon);
    bool mtComplementary = monophylyTaxon.at(0) == '!'; if(mtComplementary) monophylyTaxon = monophylyTaxon.substr(1,monophylyTaxon.size()-1);
    set<string> mtMembers = refDB.getDescendants(monophylyTaxon);
    vector<int> listeTarget;
    Node * currNode,* currAncester;
    bool auMoinsUnTransfertDetecte = false;
    bool tousLesBootstrapSuffisants = true;
    bool tousLesTauxSuffisants = true;
    vector<int> resultat;
    
    Waiter patienteur(&cout, selector.size(), '#');
    
    // on itère les familles :
    for(vector<int>::iterator laFamilleId = selector.begin() ; laFamilleId != selector.end() ; laFamilleId++) {
	Family * currFamille = (*families)[*laFamilleId];
	auMoinsUnTransfertDetecte = false; // pas encore de transfert détecté pour cette famille
	
	// pour currFamille, on doit trouver tous les noeuds appartenant au taxon target
	listeTarget = getIdWithTaxaList(currFamille->getSpTree(), &ttMembers);
	
	for(vector<int>::iterator oneNodeID = listeTarget.begin(); !auMoinsUnTransfertDetecte && oneNodeID != listeTarget.end(); oneNodeID++) {
	    currNode = currFamille->getSpTree()->getNode(*oneNodeID);
	    currAncester = currNode;
	    tousLesBootstrapSuffisants = true; // par défaut tous les bootstrap sont ok, on va opérer avec des ET logiques
	    tousLesTauxSuffisants = true;
	    Node * impasse;
	    // le bloc suivant nous permet de sortir de la monophylie
	    // dans le while, on se protège des arbres de familles monophyletiques (première condition: empêche la remontée au delà de al racine), improbable, mais ne coûte rien
	    while(currAncester->getFather() != currFamille->getSpTree()->getRootNode() && nodeOnlyContainsTheseTaxa(currAncester->getFather(),mtMembers,mtComplementary)){
		currAncester = currAncester->getFather();
		//cout << "R-";
	    }
	    // cout << endl;
	    impasse = currAncester;
	    //ATTENTION: on ne peut pas remonter plus haut que la racine de l'arbre
	    for(unsigned int i=0; tousLesBootstrapSuffisants && tousLesTauxSuffisants && currAncester != currFamille->getSpTree()->getRootNode() && i < verifDeep; i++) {
		// on remonte après le bloc
		currAncester = currAncester->getFather();
		
		// on vérifie que chaque nœud est correctement bootstrapé
		tousLesBootstrapSuffisants = tousLesBootstrapSuffisants && currAncester->getBootstrapValue() >= bootstraps.at(i);
		unsigned int size = 0;
		//cout << "Rate=" << (double)enumerateTaxon(currAncester,&tsMembers,tsComplementary,impasse,&size)/(double)size << " / " << (double)sourceRates.at(i)/100.0 << endl;
		//cout << "BS=" << currAncester->getBootstrapValue() << endl;
		
		
		if(tousLesBootstrapSuffisants) tousLesTauxSuffisants = tousLesTauxSuffisants && ((double)enumerateTaxon(currAncester,&tsMembers,tsComplementary,impasse,&size)/(double)size >= (double)sourceRates.at(i)/100.0);
		// if(tousLesTauxSuffisants) cout << "Tous taux OK" << endl; else cout << "Taux NOK !!!!" << endl;
		// if(tousLesBootstrapSuffisants) cout << "Tous bs OK" << endl; else cout << "BS NOK !!!!" << endl;
		
	    }
	    
	    auMoinsUnTransfertDetecte = auMoinsUnTransfertDetecte || (tousLesBootstrapSuffisants && tousLesTauxSuffisants && currAncester != currFamille->getSpTree()->getRootNode());
	}
	// if(auMoinsUnTransfertDetecte) cout << "Au moins 1" << endl; else cout << "Aucun" << endl;
	if(auMoinsUnTransfertDetecte) resultat.push_back(*laFamilleId);
	patienteur.step();
	
    }
    return(resultat);
}

//une fonction récursive de plus !
bool Pattern::isUnder(Node * descendant, Node * ascendant){
    if(descendant == ascendant)
	return(true);
    else
	if(descendant->hasFather())
	    return(isUnder(descendant->getFather(),ascendant));
	else
	    return(false);
}


// Fonction de détection de transfert par détection d'incongruence mesurée par saut
// 2010 Daubin, Perriere, Bigot
vector<int> xferGapDetect(map<int,Family *> & families, set<string> startTaxa){
    vector<int> famillesSelectionnees;
    
    for(map<int,Family*>::iterator currFamille = families.begin(); currFamille != families.end(); currFamille++){
	
	
    }
    
    return(famillesSelectionnees);
}

void Pattern::print(Node * noeud, int cpt){
    for(unsigned int i=0; i<cpt; i++) cout << " ";
	cout << noeud->getId();
	if(noeud->hasName()) cout << '|' << noeud->getName(); else cout << "<no name>";
	if(constraints.at(noeud->getId())->hasSpeciesRestrictions()){
	    cout << " -> " << constraints.at(noeud->getId())->getString();
	}
	if(constraints.at(noeud->getId())->getNature() == NodeConstraints::DUPLICATION)
	    cout << " <DUPLICATION>";
	else if(constraints.at(noeud->getId())->getNature() == NodeConstraints::SPECIATION)
	    cout << " <SPECIATION>";
	else cout << " <ANY>";
	    
	
	cout << endl;
	for(unsigned int i=0; i< noeud->getNumberOfSons(); i++) print(noeud->getSon(i),cpt+1);
}

bool Pattern::isTreeBinary(){
    return(tpms::TreeTools::isBinaryTree(tree.getRootNode()));
}
