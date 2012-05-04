//
// File: Waiter.cpp
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

#include <cstdlib>
#include <stdio.h>
#include <sys/ioctl.h>

#include "Waiter.hpp"

using namespace std;


Waiter::Waiter(ostream * pOutput,int pTotal, char pIndic):output(pOutput), total(pTotal), indic(pIndic),type(percent),finished(false) {
  displayStep = 1;
  direction = 1;
  done = 0;
  char * cols_c;
  string cols_s;
  
  struct winsize w;
  ioctl(0, TIOCGWINSZ, &w);

  cols = w.ws_col;
  
  if(cols > 200 || cols < 10) cols = 48;
  else if(cols >= 15) cols -= 10;
  writeFrame();
}

Waiter::Waiter(ostream * pOutput, char pIndic):output(pOutput), indic(pIndic),type(undefined),finished(false) {
  displayStep = 1;
  direction = 1;
  done = 0;
}

Waiter::~Waiter(){
  if(!finished) drawFinal();
  
}

void Waiter::step(unsigned int step){
  setDone(done+step);
}

int Waiter::getDone() {
  return(done);
}

void Waiter::setDone(int newDone){
  done = newDone;
  switch(type){
    case percent:
      if((float)done/(float)total >= ((float)displayStep)/(double)cols) {
	displayStep++;
	drawProgressBar();
      }
      break;
    case undefined:
      if(done%1000 ==0){
	if(displayStep == 20||displayStep ==0) direction = direction * -1;
	displayStep = displayStep + (direction * 1);
	drawProgressBar();
      }
      break;
  }
}	

void Waiter::writeFrame(){
  switch(type){
    case percent:
      *output << "|0%";
      for(unsigned int i = 0; i < (cols - 5); i++)
	*output << ' ';
      *output << "100%|" << endl;
      break;
    case undefined:
      *output << "|  ";
      for(unsigned int i = 0; i < (20 - 5); i++)
	*output << ' ';
      *output << "    |" << endl;
      break;
  }
}

void Waiter::drawProgressBar(){
  switch(type){
    case percent:
      *output << "[";
      for(unsigned int i = 0 ; i < displayStep; i++) {
	*output << indic;
      }
      *output << "] " << floor((float)done/(float)total*100.0) << " %\r" << flush;
      break;
    case undefined:
      *output << "[";
      for(unsigned int i = 0 ; i <= 20; i++) {
	if(i == displayStep) *output << indic;
	else *output << ' ';
      }
      *output << "] " << done << '\r' << flush;
      break;
  }
}

void Waiter::drawFinal(){
  finished = true;
  switch(type){
    case percent:
      *output << "\r[";
      for(unsigned int i = 0; i<=cols; i++)
	*output << indic;
      *output << "]      " << endl;
      break;
    case undefined:
      *output <<"[";
      for(unsigned int i = 0; i<=20; i++)
	*output << indic;
      *output << "] " << done << endl;
      break;
  } 
}
