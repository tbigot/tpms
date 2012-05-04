//
// File: tpms_mkdb.cpp
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

#include <iostream>
#include <fstream>
#include <string>
#include <boost/filesystem.hpp>

#include "classes/DataBase.hpp"
#include "classes/Pattern.hpp"
#include "classes/CmdLineArgs.hpp"
#include "classes/CandidateNode.hpp"
#include "classes/TreeTools.hpp"

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Raa/RAA.h>


using namespace std;
using namespace tpms;
using namespace bpp;

namespace fs = boost::filesystem;

int main(int argc, char *argv[]) {
    try{
    CmdLineArgs args(argc, argv,"",cerr);
    // args.print(cout);
    
    if(!args.getArg("extract-seqnames").empty()){
	cout << "extract-seqnames option : I will just extract sequences names" << endl;
	args.requireArgs("families-trees,output");
	
	// i open files
	string path(args.getArg("families-trees"));
	if(!fs::is_directory(path)){
	    cout << path << " is not a directory. I give up" << endl ;
	    return(2);
	}
	
	cout << "" << endl;
	
	// opening output file
	ofstream output(args.getArg("output").c_str(),ofstream::out);
	
	vector<string> seqList;
	
	fs::directory_iterator end;
	for (fs::directory_iterator it(path); it != end; ++it)
	{
	    cout << "  - from file " << it->path().filename();
	    if (fs::is_regular_file(it->status())){
		ifstream curFamFile(it->path().string().c_str(),ifstream::in);
		string newickLine = tpms::TreeTools::extractNewickLineFromFile(curFamFile);
		TreeTemplate<Node> * curFamTree = tpms::TreeTools::newickToTree(newickLine);
		
		if(curFamTree!=00){
		    unsigned int seqNamesNumber = 0;
		    // extracting all leaves names
		    vector<Node *> currLeaves = curFamTree->getLeaves();
		    for(vector<Node *>::iterator currLeave = currLeaves.begin(); currLeave != currLeaves.end(); currLeave++)	{
			if((*currLeave)->hasName()){
			    seqNamesNumber++;
			    seqList.push_back((*currLeave)->getName());
			}
		    }
		    cout << ": " << seqNamesNumber << " sequences found." << endl;
		}
		else
		    cout << " failed." << endl ;
		
	    }
	}
	
	cout << "I have extracted " << seqList.size() << " from the species trees to the file " << args.getArg("output") << '.' << endl;
	vector<string> namesList(seqList.size());
	
	// guess from a db ?
	if(!args.getArg("guess-from-db").empty()){
	    cout << "As asked, I will try to guess species names using ACNUC from the database " << args.getArg("guess-from-db") << endl;
	    try{
		RAA raa(args.getArg("guess-from-db"));
		
		for(unsigned int i=0; i != seqList.size(); i++){
		    RaaSeqAttributes * currSeqAttributes = raa.getAttributes(seqList.at(i));
		    namesList[i] = currSeqAttributes->getSpeciesName();
		    delete(currSeqAttributes);
		}
	    } catch(int errorNumber){
		cout << "Error connecting to acnuc database. Error number " << errorNumber << "." << endl;
	    }
	}
	
	// results -> output file
	
	for(unsigned int i=0; i != seqList.size(); i++){
	    output << seqList.at(i) << ':' << namesList.at(i) << endl;
	}
	
	return(0);
    }
    
    
    // normal program (not extract seqnames)
    args.requireArgs("sp-tree,families-trees,output,seq-to-species");
    
    cout << "* First step: loading and checking species tree" << endl;
    ifstream spTreeFile(args.getArg("sp-tree").c_str(),ifstream::in);
    
    string spTreeNw;
    getline(spTreeFile,spTreeNw);
    
    
    // tests on species tree
    // SEEMS TO BE NEWICK ?
    if(spTreeNw.find(';') == string::npos){
	cerr << "   The species tree does not seem to be Newick formatted. Sorry, I give up." << endl;
	return(1);
    }
    
    // I try to load newick
    TreeTemplate<Node> * spTree = tpms::TreeTools::newickToTree(spTreeNw);
    vector<Node *> leaves = spTree->getLeaves();
    cout << "   - This tree comprises " << leaves.size() << " leaves." << endl;
    
    set<string> spList;
    for(vector<Node*>::iterator currLeave = leaves.begin(); currLeave != leaves.end(); currLeave++){
	if((*currLeave)->hasName()) {
	    string currSpecie = (*currLeave)->getName();
	    
	    if(currSpecie.at(0) != '"' || currSpecie.at(currSpecie.size()-1) != '"'){
		cout << "Species name must be beetwin quotes at leave "<< currSpecie << ". Sorry, I give up." << endl;
		return(1);
	    }
	    spList.insert(currSpecie.substr(1,currSpecie.size()-2));
	}
    }
    
    spTreeFile.close();
    
    cout << "   - Read " << spList.size() << " species in this tree." << endl;
    cout << "   DONE \n" << endl;
    
    
    
    // NOW, LOADING MNEMONIC <-> SPECIE ASSOCIATIONS TABLE
    map<string,string> mne2spec;
    
    cout << "* Second step: loading mnemonic-specie associations" << endl;
    ifstream seqToSpeciesFile(args.getArg("seq-to-species").c_str(),ifstream::in);
    string currLine;
    string mnemo, specie;
    while(getline(seqToSpeciesFile,currLine)){
	istringstream currLineSS(currLine);
	getline(currLineSS,mnemo,':');
	getline(currLineSS,specie,':');
	
	map<string,string>::iterator found = mne2spec.find(mnemo);
	if(found != mne2spec.end()) {
	    cout << "    The mnemonic " << mnemo << " is already associated to the specie " << found->second << ". Sorry, I give up." << endl;
	    return(1);
	}
	
	if(spList.find(specie) == spList.end()) {
	    cout << "    The specie " << specie << " does not seem to be in the species tree. Sorry, I give up." << endl;
	    return(1);
	}
	mne2spec.insert(pair<string,string>(mnemo,specie));
    }
    
    // filling output file ("DB" file)
    ofstream output(args.getArg("output").c_str(),ofstream::out);
    output << "000000000 unreconciled trees" << endl;
    output << spTreeNw << endl;
    
    // WRITING FAMILIES TREES TO CONCFile
    
    unsigned int numberOfTrees = 0;
    
    fs::directory_iterator end;
    for (fs::directory_iterator it(args.getArg("families-trees")); it != end; ++it)
    // = FOR EACH TREE
    {
	vector<string> seqList;
	
	cout << "  - including family tree from " << it->path().filename() << endl;
	if (fs::is_regular_file(it->status())){
	    ifstream curFamFile(it->path().string().c_str(),ifstream::in);
	    string newickLine = tpms::TreeTools::extractNewickLineFromFile(curFamFile);
	    TreeTemplate<Node> * curFamTree = tpms::TreeTools::newickToTree(newickLine);
	    
	    if(curFamTree!=00){
		// extracting all leaves names
		vector<Node *> currLeaves = curFamTree->getLeaves();
		
		// EXTRACTING LEAVES NAMES FROM THE TREE -> seqList
		for(vector<Node *>::iterator currLeave = currLeaves.begin(); currLeave != currLeaves.end(); currLeave++)	{
		    if((*currLeave)->hasName()){
			string seqName = (*currLeave)->getName();
			if(mne2spec.find(seqName) != mne2spec.end()){
			    seqList.push_back(seqName);
			} else {
			    cout << "   '" << seqName << "' not found in seqNames <-> species association. This family will not be included." << endl;
			    continue;
			}
		    }
		}
		
		numberOfTrees++;
		// LEAVES NAMES ARE EXTRACTED, WRITING “COMMENT” PART OF THE NEWICK CODE
		// printing file name
		output << it->path().filename() << endl;
		output << '[' << endl;
		for(vector<string>::iterator currSeqName = seqList.begin(); currSeqName != seqList.end(); currSeqName++)
		{
		    output << *currSeqName << '"' << mne2spec[*currSeqName] << '"' << endl;
		}
		output << ']' << endl;
		output << newickLine << endl;
		
		
	    }
	    else
	    {
		cout << "   failed to load tree. This tree will not be included." << endl ;
		continue;
	    }
	    
	}
    } // we treated all the trees filesystem
    // rewriting the number of trees in the preamble of the file
    output.seekp(0);
    string numberOfTreesStr;
    ostringstream numberOfTreesSs;
    numberOfTreesSs << numberOfTrees;
    numberOfTreesStr = numberOfTreesSs.str();
    for(; numberOfTreesStr.size()<9; numberOfTreesStr = "0" + numberOfTreesStr);
    
    output << numberOfTreesStr;
 
    
    output.close();
    
    
    
    
    
    
     }catch(string missingArguments){
	cout << "Missing arguments: " << missingArguments.substr(0,missingArguments.size()-2) << '.' << endl;
	return(1);
    }
    
}
    