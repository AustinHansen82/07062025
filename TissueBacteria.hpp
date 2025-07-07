
#ifndef TissueBacteria_hpp
#define TissueBacteria_hpp

#include <stdio.h>
#include "Bacteria.hpp"
#include "Diffusion2D.hpp"
#include "ranum2.h"

#endif /* TissueBacteria_hpp */

enum ChemotacticMechanism
{
    observational = 0 ,
    metabolism = 1
    
};
enum InitialCondition
{
    uniform = 0 ,
    circular = 1 ,
    swarm = 2 ,
    center = 3,
    alongNetwork = 4
};




class TissueBacteria
{
public:
    //---------------------------- Parameters and sub-classes ------------------------------
    bacterium bacteria[nbacteria];
    TissueGrid tGrids ;

    int machineID = 3 ;
    string folderName = "./animation/machine" + to_string(machineID) + "/" ;
    string statsFolder = "./dataStats/machine" + to_string(machineID) + "/" ;
    string animationName = "Bacteria_" ;
    
    bool inLiquid = true ;
    bool PBC = true ;
    ChemotacticMechanism chemotacticMechanism = observational ;
    InitialCondition initialCondition = circular ;
    //Used to calibrate durations with experimental data
    bool run_calibrated = 1 ;
    double lognormal_run_m = 0.7 ;
    double lognormal_run_s = 0.6 ;
    double lognormal_run_a = 1.0 ;
    double lognormal_wrap_m = -0.4 ;
    double lognormal_wrap_s = 0.9 ;
    double lognormal_wrap_a = 1.0 ;
    std::default_random_engine runDuration_seed ;
    std::default_random_engine wrapDuration_seed ;
    int seed = 123456789;
    std::lognormal_distribution<double> runDuration_distribution ;
    std::lognormal_distribution<double> wrapDuration_distribution ;
    
    std::mt19937 reversal_rng{std::random_device{}()};
    std::uniform_real_distribution<double> unif_distribution{0.0, 1.0};
    std::mt19937 multiplier_rng{std::random_device{}()};
    std::uniform_int_distribution<int> unif_int_distribution{0, 1};
    
    //Needed to calibrate angle distributions with experimental data.( Not calibrated)
    double normal_turnAngle_mean = 0.0 ;
    double normal_turnAngle_SDV = 3.1415/6.0 * (1.0/3.0) / ( (10.0 * 1.0) * 0.3 )  ;
    double normal_wrapAngle_mean = (3.1415/2.0)/ ( (10.0 * 1.0 ) * 0.3 ) ;
    double normal_wrapAngle_SDV = (3.1415/2.0) * (1.0/3.0) / ( (10.0 * 1.0) * 0.3 ) ; //0.3 is motorForce_Amplitude ;
    std::default_random_engine turnAngle_seed ;
    std::default_random_engine wrapAngle_seed ;
    std::normal_distribution<double> turnAngle_distribution ;
    std::normal_distribution<double> wrapAngle_distribution ;

    double bacteriaLength = 2 ;
    double bending_Stiffness = 0.5 ;                      // bending constant
    double PI = 3.1415 ;               // prefered angle, Pi
    double linear_Stiffness =  5.0 ;                      // linear spring constant (stiffness)
    double equilibriumLength = bacteriaLength/(nnode-1) ;                      // equilibrium distance of spring
    //parameters to controll Lennard-Jones interaction between bacteria
    double lj_Energy = 0.0014 ;
    double lj_rMin = 0.6 ;
    bool lj_Interaction = true ;    //Bacteria won't interact with one another if this is false
    
    double updatingFrequency = 10.0 ;           //Frequency of updating reversal period with chemotax
    
    double motorForce_Amplitude = (nnode) * 0.3/nnode ;                           //   Ft/(N-1)
    
    // bacterial motor is more efficient in liquid. ( Currently not in use)
    double motorEfficiency = 1.0 ;
    double motorEfficiency_Liquid = 1.0 ;
    double motorEfficiency_Agar = 0.4 ;
    
    double proteinExchangeRate = 0.3 ;       // protein exchange rate, used for Myxo study( Abu's paper)
    long  idum = (-799);    // used for random generator in Gaussian distribution( ranum2 file)
    
    int shiftx = 0 ;        //this value changes in the MinDistance function, call MinDistance before using this value
    int shifty = 0 ;        //this value changes in the MinDistance function, call MinDistance before using this value
    //Size of the domain
    double domainx = 2000.0 ;
    double domainy = 2000.0 ;
    
    //Grid size to calculate liquid/slime values.
    //It should be comparable with the size of bacteria. ( Smaller than grid size for chemical concentration)
    double dx= 0.5 ;
    double dy= 0.5 ;
    int nx = static_cast<int>(domainx/dx) ;               //number of grids in X-axis
    int ny = static_cast<int> (domainy/dy) ;              //number of grids in Y-axis
    int nz = 1 ;
    
    double initialTime = 4.0 ;
    double runTime = 200.0 ;
    double dt=0.000005 ;
    double initialStep = 0.0001 ;
    int index1 = 0 ;    
    
    double slimeSecretionRate = 0.25 ;                      // slime rate of production
    double slimeDecayRate = 0.05 ;                          // slime decay rate
    double slimeEffectiveness = 0.8 ;           // needed for Bacterial_HighwayFollowing, used to be 5
    
    //These constants are used to calculate slime[m][n] and Bacterial_HighwayFollowing
    double liqBackground = 1.0 ;
    double liqLayer = 10.0  ;                  // liquid layer around fungi network
    double liqHyphae = 0.1 ;                    // liquid value for where hyphae is located ( physical barrier)
    int searchAreaForSlime = static_cast<int> (round(( bacteriaLength )/ min(dx , dy))) ;
    double Slime_CutOff = 1.3 ;                 //Threshold to decide bacteria is attached to fungi or not
    vector<vector<double> > slime ;
    vector<vector<double> > viscousDamp ;       //2D grid storing damping coefficient based on slime value
    vector<vector<double> > chemoProfile ;      //Store chemoattractant concentration after diffusion
    vector<vector<double> > sourceChemo ;       //store the source locations in TissueBacteria class
    vector<double> sourceProduction ;           //store the source production rate in TissueBacteria class
    bool sourceAlongHyphae = true ;

    double turnPeriod = 0.1 ;       //duration that bacteria keep changing direction after reversal events
    
    double reversalRate = 1.0/ 5.0 ;    // reversal frequency, needed for uncalibrated durations
    // Bacteria won't reverse or enter to wrap mode until it stays "minimumRunTime" in the run mode
    double minimumRunTime = 0.3 ;             //this has to be larger than turnPeriod and wrapPeriod
    
    // Parameters to control chemotaxis in observational model(classic).
    // Effective range depends on chemoattractant concentration.
    double chemoStrength_run = 5000.0 ;
    double chemoStrength_wrap = 5000.0 ;
    
    //Defining a hard agar around the domain to prevent bacteria from periodic boundary condition(PBC)
    double agarThicknessX = domainx/20 ;
    double agarThicknessY = domainy/20 ;
    
    
    double eta_background = 0.02 ;
    double eta_Barrier = 4.0 ;
    
    //vectors to hold location of the grids. It is used for visualization of slime[m][n]
    vector<double> X ;
    vector<double> Y ;
    vector<double> Z ;
    
    vector<vector<int> > visit ;       //Counts the number time that part of a bacteria was in a grid
    int fNoVisit = 0 ;
    int numOfClasses = 1 ;
    vector<double> frequency ;
    
    //Find the grids with slime secreted in them, used for Myxo study
    vector<vector<int> > surfaceCoverage ;
    double coveragePercentage = 0 ;
    double averageLengthFree = 0.0 ;        //pili average free length
    int nAttachedPili = 0 ;                 //total number of attached pili
    
    //---------------------------- Functions -----------------------------------------------
    TissueBacteria () ;
    
    void Initialze_AllRandomForce () ;
    void UpdateTissue_FromConfigFile () ;
    void Bacteria_Initialization () ;
    // Bacteria_Initialization Call one of the following initialization functions:
    void Initialization () ;
    void Initialization2 () ;
    void CircularInitialization () ;
    void SwarmingInitialization () ;
    void CenterInitialization () ;
    void AlongNetworkInitialization () ;
    
    void Initialize_Pili () ;
    void InitializeMatrix () ;
    
    //Used in Myxo study
    //Following functions calculate how many time a grid had bacteria in it
    // How many of the grids are covered by bacteria
    void VisitsPerGrid () ;
    void Update_SurfaceCoverage () ;
    void Update_SurfaceCoverage (ofstream ProteinLevelFile ,ofstream FrequencyOfVisit ) ;
    int PowerLawExponent () ;
    
    void Update_LJ_NodePositions () ;
    void Initialize_ReversalTimes () ;
    vector<vector<double> > TB_Cal_ChemoDiffusion2D (double xMin, double xMax, double yMin, double yMax,int nGridX, int nGridY ,vector<vector<double> > sources, vector<double> pSource, Chemo_Profile_Type profileType) ;
    void Cal_AllLinearSpring_Forces() ;
    
    //Does not consider periodic boundary effect.
    double Cal_NodeToNodeDistance (int b1 ,int n1, int b2, int n2) ; // bacteria, node, bacteria, node
    void Cal_AllBendingSpring_Forces() ;
    double Cal_Cos0ijk_InteriorNodes(int i,int m) ;
    
    //Calculate point to point distance considering periodic boundary condition
    double Cal_PtoP_Distance_PBC (double x1, double y1, double x2, double y2, int nx , int ny) ;
    //Calculate the minimum distance of 2 points( 2 bacteria or their duplicates)
    double Cal_MinDistance_PBC (double x1 , double y1 , double x2 , double y2) ;
    
    // Update the bacteria-bacteria connection. Used for ProteinExchange() in Myxo study( Abu's paper)
    void Update_BacterialConnection () ;
    double Cal_AllBacteriaLJ_Forces() ;
    void PiliForce () ;
    void TermalFluctiation_Forces() ;
    double Update_LocalFriction (double, double ) ;
    void Update_ViscousDampingCoeff ();
    
    //Apply equal force to each nodes. Direction is toward the next node.
    void Cal_MotorForce () ;
    void PositionUpdating (double t) ;
    //Bacteria follow the highway ( eithier liquid around hyphae of slime produced by other bacteria)
    double Bacterial_HighwayFollowing (int i, double R)  ;      // i th bacteria, radius of the search area
    
    void Reverse_IndividualBacteriaa (int i) ;
    void Check_Perform_AllReversing_andWrapping() ;
    void Update_Bacteria_AllNodes () ;          //Update allNodes with nodes and ljNodes
    void Initialize_BacteriaProteinLevel () ;
    
    void Bacterial_ProteinExchange () ; //Exchange protein between bacteria in Myxo study
    void Update_NodeLevel_ProteinConcentraion () ;
    void BacterialVisualization_ParaView () ;
    void ParaView_Liquid () ;
    
    // Make a duplicate of the bacteria that are leaving the domain
    void Bacteria_CreateDuplicate () ;
    // Merge the bacteria and its duplicate when it is completely in the domain
    void Bacteria_MergeWithDuplicate () ;
    
    //slime[][] would be calculated based on where bacteria were moving. Used for Myxo study
    //Include "bacterial secretion" and "degradation over time"
    void Update_SlimeConcentration () ;
    
    void Handle_BacteriaTurnOrientation () ;   // Bacteria change orientation after reverse or during the wrap mode
    double Cal_BacteriaOrientation (int i) ;
    
    //Calculate the impact of chemo-gradient on reversal based on observational model
    void UpdateReversalFrequency () ;
    
    //Each bacteria has its own max duration in absence of chemoattractant based on distribution observed in experiments
    void Update_BacteriaMaxDuration () ;
    
    //Calculate the change in chemoattractant for i-th bacteria.
    // It is used for changing the bacteria.reversalPeriod
    double Cal_ChemoGradient (int i) ;
    double Cal_ChemoGradient2 (int i) ;     //Current working version
    
    
    void WriteTrajectoryFile () ;
    void WriteNumberReverse () ;
    
    //The entire hyphae width is filled with liquid. Liquid decreass with distance from hyphae
    void Find_FungalNetworkTrace (Fungi tmpFng) ;
    
    //Update the liquid concentration based on relative location with respect to fungi network
    // This function "determine" the highway for bacterial motion near fungi
    void Find_FungalNetworkTrace2 (Fungi tmpFng) ;         // bacteria follows the surrounding.
    
    //Find the coordinates and secretion rates of chemoattractant sources based on our assumption( tips vs unirform along fungi)
    vector<vector<double> > Find_secretion_Coord_Rate (Fungi tmpFng) ;
    
    void Pass_PointSources_To_Bacteria (vector<vector<double> > sourceP) ;
    
    //write information of the baceria art the current time including its physical and chemical environment
    void WriteBacteria_AllStats () ;
    void WriteReversalDataByBacteria(int i);
    void WriteWrapDataByBacteria(int i);
    void WriteWrapDataByBacteria2(int i);
    void WriteReversalForce(int i);
    void WriteReversalForce2(int i);
    
    void Initialize_Distributions_RNG () ;     // Used to calibrate durations to experimental data
    
    //---------------------------- Metabolism Functions -----------------------------------------------
    void Update_MotilityMetabolism (double tmpDt) ;
    void Update_MotilityMetabolism_Only (double tmpDt) ;
    void Update_MotilityMetabolism_Only2 (double tmpDt) ;
    void WriteSwitchProbabilities () ;
    void WriteSwitchProbabilitiesByBacteria();
    void Update_MM_Legand () ;

    
};
