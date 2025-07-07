

#include "TissueBacteria.hpp"
//#include "SignalCommon.hpp"
#include <chrono>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <limits.h>
#include <cmath>
using namespace std;


void initializeSlurmConfig(int argc, char* argv[]) ;

//-----------------------------------------------------------------------------------------------------
//domain and medium properties

GlobalConfigVars globalConfigVars ;
TissueBacteria tissueBacteria ;


//-----------------------------------------------------------------------------------------------------

int main (int argc, char* argv[])
{

   auto start = std::chrono::high_resolution_clock::now() ;
   
   
   initializeSlurmConfig(argc, argv);
   tissueBacteria.UpdateTissue_FromConfigFile() ;
   Fungi fungi;
  
   srand(time (0)) ;
   cout<<"Search area for slime is "<<tissueBacteria.searchAreaForSlime<<endl ;
    
    //-----------------------------Auto Saving Scripts -----------------------------------------------
    
    char currentDir[PATH_MAX]; // defines a predefined variable currentDir of size of maximal path length
    if (getcwd(currentDir, sizeof(currentDir)) != NULL) { //getcwd is a unix function - accepts a pointer to the buffer where the workind directory path will be stored and the maximum size
        std::cout << "Current directory: " << currentDir << std::endl;
    } else {
        std::cerr << "Error getting current directory." << std::endl;
        return 1;
    }
       
    std::string scriptPath = "\"" + std::string(currentDir) + "/Auto_Save_Model.sh\""; // the backslashes \ surround string in double quotes
    std::cout << "Executing script at: " << scriptPath << std::endl;

    // Also encapsulate the tissueBacteria.statsFolder in double quotes
    string fullCommand = scriptPath + " \"" + tissueBacteria.statsFolder + "\"";
    system(fullCommand.c_str()); // Execute the bash script
    

    //---------------------------- Defining output files and directories ------------------------------
    ofstream ProteinLevelFile ;
    ProteinLevelFile.open(tissueBacteria.statsFolder + "ProteinLevelFile.txt") ;
    ofstream FrequencyOfVisit ;
    FrequencyOfVisit.open(tissueBacteria.statsFolder +"FrequncyOfVisit.txt") ;
    ofstream strSwitchP (tissueBacteria.statsFolder +"SwitchProbablities.txt" ) ;
    ofstream trajectories ( tissueBacteria.statsFolder +"trajectories.txt") ;
    
    cout<<"program is running"<<endl  ;

    //Print Header to SwitchProbabilities Data File 
    for (int i=0; i<nbacteria; i++)
    {
    ofstream strSwitchP2;
    strSwitchP2.open(tissueBacteria.statsFolder + "WriteSwitchProbabilities_Bacteria" + to_string(i) + ".txt", ios::app);
    {
       // strSwitchP << bacteria[i].motilityMetabolism.switchProbability << '\t'  ;
        strSwitchP2 << setw(10) << "maxRunDurat" << '\t' 
                    << setw(10) << "receptAct" << '\t'
                    << setw(10) << "legand" << '\t'
                    << setw(10) << "methylation" << '\t'
                    << setw(10) << "changeRate" << '\t'
                    << setw(10) << "LegandEnergy" << '\t'
                    << setw(10) << "MethylEnergy" << '\t'
                    << setw(10) << "switchProb" << '\t'
                    << setw(10) << "switchMode" << '\t';
                    
    }
    strSwitchP2<< endl ;
    }

   //--------------------------- Time controlling parameters ------------------------------------------
    double nt= tissueBacteria.runTime/tissueBacteria.dt + tissueBacteria.initialTime/tissueBacteria.initialStep ;                   //total number of steps
    nt = std::round(nt);
    nt =static_cast<int>(nt) ;
    double initialNt = tissueBacteria.initialTime/tissueBacteria.initialStep ;                     // run time for initialization
    initialNt = std::round(initialNt);
    initialNt =static_cast<int>(initialNt) ;
    //    int inverseInitialStep = static_cast<int>(initialNt/initialTime) ;  // used for visualization(BacterialVisualization_ParaView)
    int inverseDt =static_cast<int>((nt-initialNt)/(tissueBacteria.runTime)) ;              // used for visualization(BacterialVisualization_ParaView)
    // each frame is 0.1 second
    inverseDt = inverseDt/ tissueBacteria.updatingFrequency ;
   
   //--------------------------- Initializations -----------------------------------------------------
    tissueBacteria.Bacteria_Initialization() ;
    tissueBacteria.Update_LJ_NodePositions() ;
    tissueBacteria.Initialize_BacteriaProteinLevel() ;
    tissueBacteria.Initialize_Distributions_RNG() ;
    tissueBacteria.Update_BacteriaMaxDuration() ;
    tissueBacteria.Initialize_ReversalTimes() ;
    tissueBacteria.Initialize_Pili () ;
    tissueBacteria.InitializeMatrix() ;
    tissueBacteria.Initialze_AllRandomForce() ;

    fungi = driver(fungi) ;
    vector<HyphaeSegment> hyphaeSegments_main = fungi.hyphaeSegments ;
    
    //Bacteria would try to follow hyphae as a highway
    tissueBacteria.sourceAlongHyphae = false ;
    tissueBacteria.Find_FungalNetworkTrace2(fungi) ;
   
   //--------------------------- Chemical diffusion setup -----------------------------------------------------
   
   
    vector<vector<double> > pointSources ;
    pointSources = tissueBacteria.Find_secretion_Coord_Rate(fungi) ;
    fungi.WriteSourceLoc( pointSources) ;
    
   //Solve diffusion equation using Euler method
   tissueBacteria.chemoProfile = tissueBacteria.TB_Cal_ChemoDiffusion2D(0.0, tissueBacteria.domainx, 0.0, tissueBacteria.domainy ,tissueBacteria.tGrids.numberGridsX , tissueBacteria.tGrids.numberGridsY ,pointSources, tissueBacteria.sourceProduction, tissueBacteria.tGrids.chemo_profile_type ) ;
   //Save source locations in bacteria class. Used to check if the bacteria is within a source region or not
   tissueBacteria.Pass_PointSources_To_Bacteria(pointSources) ;
    
    // Change damping coefficient based on the level of liquid/slime in the underlying grid
    tissueBacteria.Update_ViscousDampingCoeff ();
   
   
   //--------------------------- Main loop -----------------------------------------------------------
    for (int l=0; l< (nt+1); l++)
    {
        //Initialization phase. Let the bacteria relax
        if (l < initialNt)
        {
            
            tissueBacteria.Update_Bacteria_AllNodes() ;
            tissueBacteria.Cal_AllLinearSpring_Forces () ;
            tissueBacteria.Cal_AllBendingSpring_Forces() ;
            //tissueBacteria.u_lj() ;
            //tissueBacteria.TermalFluctiation_Forces() ;
            /*
             if (l%inverseInitialStep==0)
             {
             BacterialVisualization_ParaView() ;
             tissueBacteria.ParaView_Liquid() ;
             }
             */
            tissueBacteria.Update_MotilityMetabolism_Only2(.01) ;
            tissueBacteria.WriteSwitchProbabilitiesByBacteria();
            tissueBacteria.PositionUpdating(tissueBacteria.initialStep) ;
            
        }
        
        
        else
        {
            tissueBacteria.Check_Perform_AllReversing_andWrapping() ;
            tissueBacteria.Update_Bacteria_AllNodes() ;
            tissueBacteria.Cal_AllLinearSpring_Forces () ;
            tissueBacteria.Cal_AllBendingSpring_Forces() ;
            tissueBacteria.Cal_AllBacteriaLJ_Forces() ;
            //tissueBacteria.TermalFluctiation_Forces() ;
            tissueBacteria.Cal_MotorForce () ;
            tissueBacteria.Handle_BacteriaTurnOrientation() ;
            //  tissueBacteria.SlimeTrace() ;
            //  tissueBacteria.Update_BacterialConnection ;
            //  tissueBacteria.Bacterial_ProteinExchange() ;
            //  tissueBacteria.PiliForce() ;
            
            
            if (l%1000==0 && l%inverseDt!=0 && tissueBacteria.chemotacticMechanism == metabolism)
            {
               cout<<(l-initialNt)/inverseDt<<endl ;
               //tissueBacteria.UpdateReversalFrequency() ;       // test Effect of chemoattacrant on reversal motion
               tissueBacteria.Update_MotilityMetabolism_Only(.01) ;
            }
            //Write the outputs and call chemotaxis related functions
            else if  (l%inverseDt==0)
            {
               //tissueBacteria.Update_SurfaceCoverage(ProteinLevelFile, FrequencyOfVisit) ; //This would slow down the code bc of high number of grids
                tissueBacteria.ParaView_Liquid() ;
               tissueBacteria.BacterialVisualization_ParaView ()  ;
               tissueBacteria.WriteTrajectoryFile() ;
               cout<<(l-initialNt)/inverseDt<<endl ;
               tissueBacteria.WriteBacteria_AllStats() ;
               tissueBacteria.WriteSwitchProbabilitiesByBacteria();
               
                //    cout << coveragePercentage <<endl ;
               
               // reducing time step in the simulation will result in changing the reversal frequencies. Lower number as inputs
               tissueBacteria.UpdateReversalFrequency() ;       // test Effect of chemoattacrant on reversal motion
               tissueBacteria.Update_MotilityMetabolism(.01) ;
               
            }
            tissueBacteria.PositionUpdating(tissueBacteria.dt) ;
        }
        
    }
    
   tissueBacteria.WriteNumberReverse() ;
   std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
   std::chrono::seconds duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
   cout << "Time to simulate "<<nbacteria<<" bacteria is : " << duration.count() << " seconds" << endl;
    cout<<"No Bug"<<endl ;
    return 0 ;
}


//------------------------------------------------------------------------------------------------------------
void initializeSlurmConfig(int argc, char* argv[]) {
   
   // read configuration.
   ConfigParser parser;
   std::string configFileNameDefault = "./resources/bacteria_M.cfg";
   globalConfigVars = parser.parseConfigFile(configFileNameDefault);
   std::string configFileNameBaseL = "./resources/test";
   std::string configFileNameBaseR = ".cfg";

   // Unknown number of input arguments.
   if (argc != 1 && argc != 3) {
      std::cout << "ERROR: Incorrect input argument count.\n"
            << "Expect either no command line argument or three arguments"
            << std::endl;
      exit(0);
   }
   // one input argument. It has to be "-slurm".
   else if (argc == 3) {
      if (strcmp(argv[1], "-slurm") != 0) {
         std::cout
               << "ERROR: one argument received from commandline but it's not recognized.\n"
               << "Currently, the argument value must be -slurm"
               << std::endl;
         exit(0);
      } else {
         std::string configFileNameM(argv[2]);
         std::string configFileNameCombined = configFileNameBaseL
               + configFileNameM + configFileNameBaseR;
         parser.updateConfigFile(globalConfigVars, configFileNameCombined);
      }
   }
   // no input argument. Take default.
}
