#include "CmdLineArgs.hpp"

#include <sstream>
#include <map>
#include <vector>


using namespace std;


namespace tpms{
std::string CmdLineArgs::getArg(std::string argName)
{
    map<string,string>::iterator foundArg = argsList.find(argName);
    if(foundArg == argsList.end()) return "";
    else return foundArg->second;
}

CmdLineArgs::CmdLineArgs(int nbArgs, char* mainArgs[], string reqArgs, ostream & errorLog) throw(string)
{
    binPath = string(mainArgs[0]);
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
	    argsList.insert(pair<string,string>(argName.substr(1,argName.size()-1), argValue));
	else
	    errorLog << "Arguments must start with a dash in:\n" << argName << endl;
    }
    requireArgs(reqArgs);
    
}

void CmdLineArgs::print(ostream & out){
    out << "### Argument list control ###" << endl;
    out << "## Binary = " << binPath << endl;
    
    for(map<string,string>::iterator currArgPair = argsList.begin(); currArgPair != argsList.end(); currArgPair++){
	out << "## ";
	if(reqArgs.find(currArgPair->first) != reqArgs.end()) out << '*';
	out << currArgPair->first << '=' << currArgPair->second << endl;
    }
}

vector<string> CmdLineArgs::getMissingArgs(){
    //cette fonction vérifie simplement l'inclusion de reqArgs dans argsList.
    vector<string> missingArguments;
    for(set<string>::iterator currRA = reqArgs.begin() ; currRA != reqArgs.end(); currRA++){
	if(argsList.find(*currRA) == argsList.end()) {
	    missingArguments.push_back(*currRA);
	}
    }
    return(missingArguments);
}

void CmdLineArgs::addReqArgs(std::string reqArgs){
    string currReqArg;
    istringstream reqArgsList(reqArgs);
    while(getline(reqArgsList,currReqArg,','))
	this->reqArgs.insert(currReqArg);
}

void CmdLineArgs::requireArgs(string reqArgs)
{
    addReqArgs(reqArgs);
    vector<string> missingArguments = getMissingArgs();
    if(missingArguments.size() != 0){
	ostringstream missingArgumentsSs;
	for(vector<string>::iterator it = missingArguments.begin(); it != missingArguments.end(); it++)
	    missingArgumentsSs << *it << ", ";
	throw(missingArgumentsSs.str());
    }
}

}
