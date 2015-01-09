//
// File: Pattern.cpp
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


#include "Pattern.hpp"
#include "NodeConstraints.hpp"
#include "TreeTools.hpp"
#include "Waiter.hpp"
#include "CandidateNode.hpp"
#include "Taxon.hpp"
#include "Family.hpp"


#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>



using namespace std;
using namespace bpp;
using namespace tpms;


namespace tpms{
Pattern::Pattern(TreeTemplate<Node> &tree, DataBase &db): db(db), tree(tree), ok(true) {
    
    extractConstraints();
    cout << "\n -- Constraints visual representation:" << endl;
    toString(tree.getRootNode(),0,std::cout);
    mapping_NodesToMaxDepths.resize(tree.getNumberOfNodes());
    mapNodeToMaxDepth(tree.getRootNode());
}

Pattern::~Pattern(){
    for(vector<NodeConstraints *>::iterator cc=constraints.begin(); cc != constraints.end(); cc++)
	delete(*cc);
}


unsigned int Pattern::mapNodeToMaxDepth(Node* node)
{
    unsigned int maxDepth = 0;
    if(!node->isLeaf()){
        vector<Node*> sons = node->getSons();
        for(vector<Node*>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++){
            unsigned int currMaxDepth = mapNodeToMaxDepth(*currSon);
            if(currMaxDepth > maxDepth) maxDepth = currMaxDepth;
        }
    }
    mapping_NodesToMaxDepths.at(node->getId()) = maxDepth;
    return(maxDepth);
}


void Pattern::extractConstraints(){
    // les noms de feuilles peuvent contenir des contraintes entre accolades
    // la contrainte est portée par le nœud fils alors qu'il s'agit d'une contrainte sur la branche menant au nœud fils
    // et ceci pour des raisons de commodité.
    
    // resize et remplissage : ici, on sait ce qu'on fait.
    constraints.resize(tree.getNumberOfNodes());
    for(unsigned int i = 0; i!= constraints.size();i++)
	constraints.at(i) = new NodeConstraints(db);
    
    vector<Node *> noeuds = tree.getNodes();
    for(vector<Node *>::iterator unNoeud = noeuds.begin(); unNoeud != noeuds.end(); unNoeud++){
	string initName;
	if((*unNoeud)->hasName()) initName = (*unNoeud)->getName();
	string nodeName = initName;
	// if a pattern node has a name, it means it has constraints to extract
	tpms::NodeType nodeType;
	if((*unNoeud)->getNumberOfSons() == 0) nodeType = tpms::LEAF; else nodeType = tpms::NODE;
	(constraints.at((*unNoeud)->getId()))->setConstraints(db,initName,nodeType);
	
    }
    
    // now checking constraints
    for(vector<NodeConstraints*>::iterator cc = constraints.begin(); ok && cc != constraints.end(); cc++){
        if(*cc != 00)
            ok &= (*cc)->isOk();
    }
}


bool Pattern::isOk(){
    return(ok);
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
// 	for(vector<Taxon*>::iterator oneNodeList = species.begin(); oneNodeList != species.end() && tousNoeudsOntEspece ; oneNodeList++) {
// 	    // la valeur de chaque élément de la map est un vecteur de strings correspondant aux espèces
// 	    // -- on est dans un noeud du pattern avec une liste d'espèces candidates. On doit trouver si au moins une espèce est présente dans la famille.
// 	    bool uneEspeceTrouvee = false; // y a-t-il au moins une espèce en commun ?
// 	    // de deux choses l'une : soit on est dans le cas d'un nœud classique avec un nom de taxon, dans ce cas, on doit trouver au moins un représentant de ce taxon dans la famille
// 	    // soit c'est un nœud négatif, qui commence par un «!», auquel cas, il faudrait trouver au moins une espèce qui n'appartient pas au taxon donné. On n'a aucun vecteur qui contienne ce type de liste d'espèces. Souvent le complémentaire d'un taxon est vaste, donc on part sur le principe qu'il y a un représentant dans la liste. Au pire on aura fait du pattern matching pour rien, mais peu probable
// 	    
// 	    if(tree.getNode(oneNodeList->first)->getName().at(0) == '!') uneEspeceTrouvee = true;
// 	    
// 	    for(set<string>::iterator thePatternSpecie = oneNodeList->second.begin(); thePatternSpecie != oneNodeList->second.end() && !uneEspeceTrouvee; thePatternSpecie++) {
// 		
// 		
// 		if((*theFamily)->getSpecies()->find(*thePatternSpecie) != (*theFamily)->getSpecies()->end())
// 		    uneEspeceTrouvee = true;
// 	    }
// 	    tousNoeudsOntEspece = tousNoeudsOntEspece && uneEspeceTrouvee;
// 	    // si le nœud qui vient d'être fait n'a pas d'espèce, tousNoeudsOntEspece devient faux
// 	    // et dans ce cas, on arrête la boucle
// 	}
// 	// en arrivant ici, soit on a passé tous les test avec succès tousNoeudsOntEspece=vrai, soit on s'est arrêté en route parce qu'un nœud n'avait pas d'espèce dans la famille
// 	
	
	if(tousNoeudsOntEspece) {
	    CandidateNode * candRoot = new CandidateNode();
	    if(patternMatch(**theFamily,(*theFamily)->getTree()->getRootNode(),tree.getRootNode(),candRoot)){
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


// bool Pattern::nodeOnlyContainsTheseTaxa(Node * localRoot, set<string> & taxonMembers, bool invert){
//     // monophyletique 
//     //FONCTION RÉCURSIVE
//     
//     // - cas de base : feuille.
//     // Soit elle n'appartient pas à tsMembers, soit c'est le nœud targetNode, sinon false.
//     if(isLeaf(localRoot)){
// 	if(invert){
// 	    if(taxonMembers.find(localRoot->getName()) == taxonMembers.end()) return(true); else return(false);
// 	} else {
// 	    if(taxonMembers.find(localRoot->getName()) != taxonMembers.end()) return(true); else return(false);
// 	}
//     } else { // -cas récursif = un noeud
// 	bool tousLesDescendantsOK = true; // par défaut si pas de fils, ils sont tous descendants du taxon. On opérera avec des ET
// 	for(unsigned int currFilsIndex = 0; tousLesDescendantsOK && currFilsIndex < localRoot->getNumberOfSons(); currFilsIndex++){
// 	    tousLesDescendantsOK = tousLesDescendantsOK && nodeOnlyContainsTheseTaxa(localRoot->getSon(currFilsIndex), taxonMembers, invert);
// 	}
// 	return(tousLesDescendantsOK);
// }
// }

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


// vector<int> Pattern::xferDetected(map<int,Family *> * families, vector<int> selector, string sourceTaxon, string targetTaxon, string monophylyTaxon, unsigned int verifDeep, vector<unsigned int> bootstraps, vector<unsigned int> sourceRates) {
//     // cette fonction doit trouver dans la liste famList toutes les familles qui groupent un membre de taxonTarget dans taxonSource
//     bool tsComplementary = sourceTaxon.at(0) == '!'; if(tsComplementary) sourceTaxon = sourceTaxon.substr(1,sourceTaxon.size()-1);
//     set<string> tsMembers = db.getDescendants(sourceTaxon);
//     bool ttComplementary = targetTaxon.at(0) == '!'; if(ttComplementary) targetTaxon = targetTaxon.substr(1,targetTaxon.size()-1);
//     set<string> ttMembers = db.getDescendants(targetTaxon);
//     bool mtComplementary = monophylyTaxon.at(0) == '!'; if(mtComplementary) monophylyTaxon = monophylyTaxon.substr(1,monophylyTaxon.size()-1);
//     set<string> mtMembers = db.getDescendants(monophylyTaxon);
//     vector<int> listeTarget;
//     Node * currNode,* currAncester;
//     bool auMoinsUnTransfertDetecte = false;
//     bool tousLesBootstrapSuffisants = true;
//     bool tousLesTauxSuffisants = true;
//     vector<int> resultat;
//     
//     Waiter patienteur(&cout, selector.size(), '#');
//     
//     // on itère les familles :
//     for(vector<int>::iterator laFamilleId = selector.begin() ; laFamilleId != selector.end() ; laFamilleId++) {
// 	Family * currFamille = (*families)[*laFamilleId];
// 	auMoinsUnTransfertDetecte = false; // pas encore de transfert détecté pour cette famille
// 	
// 	// pour currFamille, on doit trouver tous les noeuds appartenant au taxon target
// 	listeTarget = getIdWithTaxaList(currFamille->getSpTree(), &ttMembers);
// 	
// 	for(vector<int>::iterator oneNodeID = listeTarget.begin(); !auMoinsUnTransfertDetecte && oneNodeID != listeTarget.end(); oneNodeID++) {
// 	    currNode = currFamille->getSpTree()->getNode(*oneNodeID);
// 	    currAncester = currNode;
// 	    tousLesBootstrapSuffisants = true; // par défaut tous les bootstrap sont ok, on va opérer avec des ET logiques
// 	    tousLesTauxSuffisants = true;
// 	    Node * impasse;
// 	    // le bloc suivant nous permet de sortir de la monophylie
// 	    // dans le while, on se protège des arbres de familles monophyletiques (première condition: empêche la remontée au delà de al racine), improbable, mais ne coûte rien
// 	    while(currAncester->getFather() != currFamille->getSpTree()->getRootNode() && nodeOnlyContainsTheseTaxa(currAncester->getFather(),mtMembers,mtComplementary)){
// 		currAncester = currAncester->getFather();
// 		//cout << "R-";
// 	    }
// 	    // cout << endl;
// 	    impasse = currAncester;
// 	    //ATTENTION: on ne peut pas remonter plus haut que la racine de l'arbre
// 	    for(unsigned int i=0; tousLesBootstrapSuffisants && tousLesTauxSuffisants && currAncester != currFamille->getSpTree()->getRootNode() && i < verifDeep; i++) {
// 		// on remonte après le bloc
// 		currAncester = currAncester->getFather();
// 		
// 		// on vérifie que chaque nœud est correctement bootstrapé
// 		tousLesBootstrapSuffisants = tousLesBootstrapSuffisants && currAncester->getBootstrapValue() >= bootstraps.at(i);
// 		unsigned int size = 0;
// 		//cout << "Rate=" << (double)enumerateTaxon(currAncester,&tsMembers,tsComplementary,impasse,&size)/(double)size << " / " << (double)sourceRates.at(i)/100.0 << endl;
// 		//cout << "BS=" << currAncester->getBootstrapValue() << endl;
// 		
// 		
// 		if(tousLesBootstrapSuffisants) tousLesTauxSuffisants = tousLesTauxSuffisants && ((double)enumerateTaxon(currAncester,&tsMembers,tsComplementary,impasse,&size)/(double)size >= (double)sourceRates.at(i)/100.0);
// 		// if(tousLesTauxSuffisants) cout << "Tous taux OK" << endl; else cout << "Taux NOK !!!!" << endl;
// 		// if(tousLesBootstrapSuffisants) cout << "Tous bs OK" << endl; else cout << "BS NOK !!!!" << endl;
// 		
// 	    }
// 	    
// 	    auMoinsUnTransfertDetecte = auMoinsUnTransfertDetecte || (tousLesBootstrapSuffisants && tousLesTauxSuffisants && currAncester != currFamille->getSpTree()->getRootNode());
// 	}
// 	// if(auMoinsUnTransfertDetecte) cout << "Au moins 1" << endl; else cout << "Aucun" << endl;
// 	if(auMoinsUnTransfertDetecte) resultat.push_back(*laFamilleId);
// 	patienteur.step();
// 	
//     }
//     return(resultat);
// }

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

void Pattern::toString(ostream &outputStream){
    toString(tree.getRootNode(),0,outputStream);
}

void Pattern::toString(Node * noeud, int cpt, ostream &outputStream){
    for(unsigned int i=0; i<cpt; i++) cout << " ";
	outputStream << noeud->getId();
	outputStream << constraintsOf(noeud)->getStr();
	outputStream << endl;
	for(unsigned int i=0; i< noeud->getNumberOfSons(); i++) toString(noeud->getSon(i),cpt+1,outputStream);
}

bool Pattern::isTreeBinary(){
    return(tpms::TreeTools::isBinaryTree(tree.getRootNode()));
}


NodeConstraints* Pattern::constraintsOf(Node* node){
    return(constraints.at(node->getId()));
}

bool Pattern::patternMatchInit(Family &family, CandidateNode * initCnode){
    family.initCache();
    bool result = patternMatch(family, family.getTree()->getRootNode(),tree.getRootNode(), initCnode);
    family.clearCache();
    return(result);
}


bool Pattern::patternMatch(Family& family,Node * target, Node * pattern, CandidateNode * fatherCandidate) {
    // Dufayard et al, 2005
    // Bigot et al,2013.
    
    
    if(isLeaf(target) && isLeaf(pattern)){
	if(constraintsOf(pattern)->allows(family,target))
	{
	    CandidateNode * currCandidate = new CandidateNode(fatherCandidate,target,pattern);
	    currCandidate->confirm();
	    return(true);
	} else return(false);
    } else if(isLeaf(target) && !isLeaf(pattern)) {
	return(false);
    } else if(!isLeaf(target) && isLeaf(pattern)) {
	return(patternMatch(family,target->getSon(0),pattern,fatherCandidate)
	|| patternMatch(family,target->getSon(1),pattern,fatherCandidate));
    } else {
        // first, checking if the rest of the tree is big enough to contain the pattern
        if(family.getMaxDepthOfSubtree(target) < mapping_NodesToMaxDepths.at(pattern->getId()))
            return(false);
        
        
	// here, the pattern node can accept the target node
        // so we shortcut the pattern in left&right; same for the target tree.
	Node * tson1 = target->getSon(0);
	Node * tson2 = target->getSon(1);
	Node * pson1 = pattern->getSon(0);
	Node * pson2 = pattern->getSon(1);
	
	// we need to check if the target node is accepted (nature) by the pattern
	
	if(constraintsOf(pattern)->allows(family,target)){
            // all the tests have to be done: we are not only wondering if the family matches, but
            // what are all the matching patterns.
            bool thisNodeMatches = false;
	
	    CandidateNode * candidate = new CandidateNode(fatherCandidate, target, pattern);
	    
	    if(constraintsOf(pson1)->allowsAsSon(family,tson1) && constraintsOf(pson2)->allowsAsSon(family,tson2)
		&& patternMatch(family,tson1, pson1, candidate) && patternMatch(family,tson2, pson2, candidate) ) {
                    candidate->confirm(); 
                    thisNodeMatches = true;
		} else delete(candidate);
		
		
		
	    candidate = new CandidateNode(fatherCandidate, target, pattern);
		
	    if(constraintsOf(pson1)->allowsAsSon(family,tson2) && constraintsOf(pson2)->allowsAsSon(family,tson1)
		&& patternMatch(family,tson2, pson1, candidate) && patternMatch(family,tson1, pson2, candidate) ){
		candidate->confirm();
                thisNodeMatches = true;
		} else delete(candidate);
		    
	    // if we haven't returned yet, it means the current target node is not matching to a pattern node
	    // we have to try the pattern node in the sons
                //Only if the direct link is not required
                // DIRECT LINK MANAGED HERE
	    thisNodeMatches |= ( !constraintsOf(pattern)->isDirect() && (patternMatch(family,tson1, pattern, fatherCandidate) || patternMatch(family,tson2, pattern, fatherCandidate)));
            return(thisNodeMatches);
	} else return(false);
		
    }
    
}

}