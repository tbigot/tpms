#ifndef CMDLINEARGS
#define CMDLINEARGS

#include <set>
#include <map>
#include <string>
#include <iostream>
#include <vector>

namespace tpms{

class CmdLineArgs {
  
    private:
	// required arguments list
	std::set<std::string> reqArgs;
	// list containing arguments declared by user
	std::map<std::string, std::string> argsList;
	// path of the binary called (argv[0]);
	std::string binPath;
	// check that required arguments have all been declared
	std::vector<std::string> getMissingArgs();
	void addReqArgs(std::string reqArgs);
    
    
    public:
	CmdLineArgs(int nbArgs, char* mainArgs[], std::string reqArgs, std::ostream& errorLog) throw(std::string);
	std::string getArg(std::string argName);
	void print(std::ostream & output);
	void requireArgs(std::string reqArgs);
};
}

#else
namespace tpms{

class CmdLineArgs;
}
#endif
