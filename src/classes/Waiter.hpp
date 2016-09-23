//
// File: Waiter.hpp
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


#ifndef WAITER
#define WAITER

#include <sstream>
#include <cmath>
#include <iostream>
class Waiter {
  
  private:
    enum FrameType { percent, undefined };
    FrameType type_;
    
    std::ostream * output_;
    int total_;
    int done_;
    char indic_;
    unsigned int displayStep_;
    void writeFrame_();
    int direction_;
    void drawProgressBar_();
    bool finished_;
    float realstep_;
    unsigned int cols_;
    
    
  public:
    //constructeur à partir d'un fichier
    Waiter(std::ostream * pOutput,int pTotal, char pIndic);
    Waiter(std::ostream * pOutput, char pIndic);
    ~Waiter();
    void doStep(unsigned int step = 1);
    void setDone(int newDone);
    int getDone();
    void drawFinal();
};



#else

class Waiter;

#endif
