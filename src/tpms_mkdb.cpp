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
    