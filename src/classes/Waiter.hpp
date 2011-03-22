#ifndef WAITER
#define WAITER

#include <sstream>
#include <cmath>
#include <iostream>

class Waiter {
  
  private:
    enum FrameType { percent, undefined };
    FrameType type;
    
    std::ostream * output;
    int total;
    int done;
    char indic;
    unsigned int displayStep;
    void writeFrame();
    int direction;
    void drawProgressBar();
    bool finished;
    float realstep;
    
    
  public:
    //constructeur à partir d'un fichier
    Waiter(std::ostream * pOutput,int pTotal, char pIndic);
    Waiter(std::ostream * pOutput, char pIndic);
    ~Waiter();
    void step();
    void setDone(int newDone);
    int getDone();
    void drawFinal();
};



#else

class Waiter;

#endif
