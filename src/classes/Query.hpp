//
// File: Query.hpp
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


#ifndef QUERY_HPP_TPMS
#define QUERY_HPP_TPMS

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <set>

#include "DataBase.hpp"
namespace tpms{

class Query {
	
	private:
		std::string sourceTaxon_;
		std::string targetTaxon_;
		unsigned int verifDepth_;
		unsigned int monophilyLevel_;
		
		std::vector<unsigned int> bootstraps_;
		std::vector<unsigned int> sourceRates_;
	
	public:
		Query(std::string queryLine,std::string defaultsParam);
		
		std::string getPatternString();
		std::string getSource();
		std::string getTarget();
		unsigned int getDepth();
		unsigned int getMonophilyLevel();
		
		std::vector<unsigned int> getBootstraps();
		std::vector<unsigned int> getSourceRates();
		
    
};
}
#else
namespace tpms{

class Family;

}
#endif
