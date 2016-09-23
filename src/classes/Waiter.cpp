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


Waiter::Waiter(ostream * pOutput,int pTotal, char pIndic):output_(pOutput), total_(pTotal), indic_(pIndic),type_(percent),finished_(false) {
  displayStep_ = 1;
  direction_ = 1;
  done_ = 0;
  char * cols_c;
  string cols_s;
  
  struct winsize w;
  ioctl(0, TIOCGWINSZ, &w);

  cols_ = w.ws_col;
  
  if(cols_ > 200 || cols_ < 10) cols_ = 48;
  else if(cols_ >= 15) cols_ -= 10;
  writeFrame_();
}

Waiter::Waiter(ostream * pOutput, char pIndic):output_(pOutput), indic_(pIndic),type_(undefined),finished_(false) {
  displayStep_ = 1;
  direction_ = 1;
  done_ = 0;
}

Waiter::~Waiter(){
  if(!finished_) drawFinal();
  
}

void Waiter::step(unsigned int step){
  setDone(done_+step);
}

int Waiter::getDone() {
  return(done_);
}

void Waiter::setDone(int newDone){
  done_ = newDone;
  switch(type_){
    case percent:
      if((float)done_/(float)total_ >= ((float)displayStep_)/(double)cols_) {
	displayStep_++;
	drawProgressBar_();
      }
      break;
    case undefined:
      if(done_%1000 ==0){
	if(displayStep_ == 20||displayStep_ ==0) direction_ = direction_ * -1;
	displayStep_ = displayStep_ + (direction_ * 1);
	drawProgressBar_();
      }
      break;
  }
}	

void Waiter::writeFrame_(){
  switch(type_){
    case percent:
      *output_ << "|0%";
      for(unsigned int i = 0; i < (cols_ - 5); i++)
	*output_ << ' ';
      *output_ << "100%|" << endl;
      break;
    case undefined:
      *output_ << "|  ";
      for(unsigned int i = 0; i < (20 - 5); i++)
	*output_ << ' ';
      *output_ << "    |" << endl;
      break;
  }
}

void Waiter::drawProgressBar_(){
  switch(type_){
    case percent:
      *output_ << "[";
      for(unsigned int i = 0 ; i < displayStep_; i++) {
	*output_ << indic_;
      }
      *output_ << "] " << floor((float)done_/(float)total_*100.0) << " %\r" << flush;
      break;
    case undefined:
      *output_ << "[";
      for(unsigned int i = 0 ; i <= 20; i++) {
	if(i == displayStep_) *output_ << indic_;
	else *output_ << ' ';
      }
      *output_ << "] " << done_ << '\r' << flush;
      break;
  }
}

void Waiter::drawFinal(){
  finished_ = true;
  switch(type_){
    case percent:
      *output_ << "\r[";
      for(unsigned int i = 0; i<=cols_; i++)
	*output_ << indic_;
      *output_ << "]      " << endl;
      break;
    case undefined:
      *output_ <<"[";
      for(unsigned int i = 0; i<=20; i++)
	*output_ << indic_;
      *output_ << "] " << done_ << endl;
      break;
  } 
}
