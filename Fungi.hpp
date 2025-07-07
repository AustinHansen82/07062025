//
//  Fungi.hpp
//  Myxobacteria
//
//  Created by Alireza Ramezani on 11/19/20.
//  Copyright © 2020 Alireza Ramezani. All rights reserved.
//

#ifndef Fungi_hpp
#define Fungi_hpp


#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

#include "growthFunctions.h"

class Fungi {
public:
    //---------------------------- Parameters and sub-classes ------------------------------
    vector<HyphaeSegment> hyphaeSegments ;
    vector<vector<double> > tips_Coord ;
    vector<int > tips_SegmentID ;
    vector<int> centerConnectedSegments ;
    vector<double> Hypahe_X1Values;
    vector<double> Hypahe_X2Values;
    vector<double> Hypahe_Y1Values;
    vector<double> Hypahe_Y2Values;
    //Angle between extended segment and primary segment
    //Used for manual network
    double extendingAngle = 13.0 * (constants::pi / 180.0) ;
    // Hyphae + Liquid around it
    double hyphaeWidth = 6.0 ;
    // What fraction of hyphaeWidth is hyphae? the rest would be liquid
    double hyphaeOverLiq = 0.5 ;
    // production rate at the tips.
    
    double production = pow(10, 2) ;     //This is going to update the number used in TissueGrid
    //Include the aging effect in production rates.
    double proDecayFactor = 0.5 ;       //Currently not in use
    //Center of the fungal network
    double initX = 500.0 ;
    double initY = 500 ;
    double Bacteria_inLiquid = 1.0 ; //Parameter that takes 1 if Simulation is In Liquid and 0 if with Fungi
    int init_Count = 4 ;
    bool loading_Network = true ;
    bool branchIsTip = false ;
    
    int machineID = 1 ;
    string folderName = "./animation/machine" + to_string(machineID) + "/" ;
    string statsFolder = "./dataStats/machine" + to_string(machineID) + "/" ;
    
    //---------------------------- Functions ---------------------------------------
    Fungi () ;
    void Find_Hyphae_Tips() ;               // Store the branch points and the tips in the "tips vector"
    void Find_Hyphae_Tips2() ;              //Only the tips of fungi network are stored in the "tips vector"
    void WriteSourceLoc(vector<vector<double> >) ;
    //Need calibration if you want to use it. (proDecayFactor, HyphaeSegment.p1 & p2, and etc. )
    //Scale production rates based on the network-distance from:
    void FindProductionNetwork ();          // the tips
    void FindProductionNetwork2();          //the Center of the network
    vector<int> FindNeighborList (int i) ;
    void Rec_CalNeighborProduction (int id, vector<int> ngbrList, int nbrElement, bool connection) ;
    void Rec_CalNeighborProduction2(int id, vector<int> ngbrList, int nbrElement ) ;
    
    //Find the neighbors of the segments at the center of the network
    vector<int> FindCenterPointConnection () ;
    
    void UpdateFungi_FromConfigFile () ;
    
};


#endif /* Fungi_hpp */

