//
// File: CmdLineArgs.cpp
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

#include "CmdLineArgs.hpp"

#include <sstream>
#include <map>
#include <vector>


using namespace std;


namespace tpms{
std::string CmdLineArgs::getArg(std::string argName)
{
    map<string,string>::iterator foundArg = argsList_.find(argName);
    if(foundArg == argsList_.end()) return "";
    else return foundArg->second;
}

CmdLineArgs::CmdLineArgs(int nbArgs, char* mainArgs[], string reqArgs, ostream & errorLog) throw(string)
{
    binPath_ = string(mainArgs[0]);
    // boucle de prise en charge des arguments en c-style
    string argName,argValue;
    istringstream currArg;
    for(unsigned int i = 1; i < nbArgs; i++){
	currArg.clear();
	currArg.str(mainArgs[i]);
	getline(currArg,argName,'=');
	getline(currArg,argValue,'=');
	if(argValue.empty()) argValue = "true";
	if(argName.at(0) == '-')
	    argsList_.insert(pair<string,string>(argName.substr(1,argName.size()-1), argValue));
	else
	    errorLog << "Arguments must start with a dash in:\n" << argName << endl;
    }
    requireArgs(reqArgs);
    
}

void CmdLineArgs::print(ostream & out){
    out << "### Argument list control ###" << endl;
    out << "## Binary = " << binPath_ << endl;
    
    for(map<string,string>::iterator currArgPair = argsList_.begin(); currArgPair != argsList_.end(); currArgPair++){
	out << "## ";
	if(reqArgs_.find(currArgPair->first) != reqArgs_.end()) out << '*';
	out << currArgPair->first << '=' << currArgPair->second << endl;
    }
}

vector<string> CmdLineArgs::getMissingArgs_(){
    //cette fonction vérifie simplement l'inclusion de reqArgs dans argsList.
    vector<string> missingArguments;
    for(set<string>::iterator currRA = reqArgs_.begin() ; currRA != reqArgs_.end(); currRA++){
	if(argsList_.find(*currRA) == argsList_.end()) {
	    missingArguments.push_back(*currRA);
	}
    }
    return(missingArguments);
}

void CmdLineArgs::addReqArgs_(std::string reqArgs){
    string currReqArg;
    istringstream reqArgsList(reqArgs);
    while(getline(reqArgsList,currReqArg,','))
	this->reqArgs_.insert(currReqArg);
}

void CmdLineArgs::requireArgs(string reqArgs)
{
    addReqArgs_(reqArgs);
    vector<string> missingArguments = getMissingArgs_();
    if(missingArguments.size() != 0){
	ostringstream missingArgumentsSs;
	for(vector<string>::iterator it = missingArguments.begin(); it != missingArguments.end(); it++)
	    missingArgumentsSs << *it << ", ";
	throw(missingArgumentsSs.str());
    }
}

}
