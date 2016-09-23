//
// File: CmdLineArgs.hpp
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

#ifndef CMDLINEARGS
#define CMDLINEARGS

#include <set>
#include <map>
#include <string>
#include <iostream>
#include <vector>

// Usage example:
//  CmdLineArgs args(argc, argv, "collection,output-dir",cerr); à mettre au début du programme
// pour appeler les arguments : args.getArg("collection") -> une string est retournée


namespace tpms{

class CmdLineArgs {
  
    private:
	// required arguments list
	std::set<std::string> reqArgs_;
	// list containing arguments declared by user
	std::map<std::string, std::string> argsList_;
	// path of the binary called (argv[0]);
	std::string binPath_;
	// check that required arguments have all been declared
	std::vector<std::string> getMissingArgs_();
	void addReqArgs_(std::string reqArgs);
    
    
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
