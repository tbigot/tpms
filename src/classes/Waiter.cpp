
#include "Waiter.hpp"

using namespace std;


Waiter::Waiter(ostream * pOutput,int pTotal, char pIndic):output(pOutput), total(pTotal), indic(pIndic),type(percent),finished(false) {
  displayStep = 1;
  direction = 1;
  done = 0;
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
      if((float)done/(float)total >= ((float)displayStep)/(double)48) {
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
      for(unsigned int i = 0; i < (48 - 5); i++)
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
      for(unsigned int i = 0; i<=48; i++)
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
