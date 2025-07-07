

#ifndef Nodes_hpp
#define Nodes_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <vector>
#include <numeric>
#include <random>
#include "SignalCommon.hpp"

using namespace std ;

#define nnode 5                                                // need to be an odd number
#define nbacteria 200                                         // 2* n^2
#define points nbacteria*nnode
//#define domainx  1000.0
//#define domainy  1000.0
#define nPili   1


class node
{   public:
    //location
    double x ;
    double y ;
    //Linear spring force
    double fSpringx ;
    double fSpringy ;
    //Bending spring force
    double fBendingx ;
    double fBendingy ;
    //Force due to Lennard-Jones interaction
    double fljx ;
    double fljy ;
    //Random force due to thermal fluctuations
    double xdev ;
    double ydev ;
    //Motor force
    double fMotorx ;
    double fMotory ;
    
    double protein ;
    
    void test () ;
    
};

class pilus
{   public:
    double lFree ;                                // pili free length
    double retractionRate ;
    double vRet ;                                // retraction rate for connected pili
    double subAttachmentRate = 1 ;              //Attachment rate to the substrate
    double subDetachmentRate = 0.1 ;            //Detachment rate from the substrate
    double lContour ;
    double piliForce ;                  //magnitude of the force
    double pili_Fx ;
    double pili_Fy ;
    double xEnd ;                   //x component of the end of the pili
    double yEnd ;                   //y component of the end of the pili
    // true for pili-substrate connection and false for pili-pili connection
    bool piliSubs = true ;
    // we change this parameter for simulation
    bool attachment ;
    bool retraction ;
    
    //pili properties
    double piliMaxLength = 2.5 ;                       //pili maximum length
    double vProPili = 0.5 ;                                 // constant protrude velocity
    double vRetPili0 = 0.5 ;                                // constant retraction rate when the pili is not connected
    double kPullPili = 2000 ;
    double fStall = 180.0 ;
    
};

#endif /* Nodes_hpp */

