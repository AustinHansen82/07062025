
#ifndef Bacteria_hpp
#define Bacteria_hpp

#include "Nodes.hpp"

enum MovingDirection
{
    backward = 0 ,
    frward = 1
};
class MotilityMetabolism
{
public:
    double a0 = 0.33 ;
    double m0 = 1.0 ;
    double receptorActivity = a0 ;
    double lnA_a = -3.2 ;
    double b = 6.0 ;
    double gamma = 0.5 ;
    double switchProbability ;
    double legand = 20.0 ;              //This is going to be relaced by chemoProfile
    double methylation  = 0.93 ;
    double changeRate = 0.0;
    double km = 1.7 ;
    double kI = 18.2 ;       //18.2 or 1.0
    double kA = 3000.0 ;       //3.0 or 1.0
    double N = 6 ;
    double kR = 0.005 ;
    double kB = 0.010 ;
    double MethylEnergy ;
    double LegandEnergy ;
    bool switchMode = false ;
    
    MotilityMetabolism ();
    double Cal_MethylationEnergy ();
    double Cal_LegandEnergy () ;
    double Cal_ReceptorActivity () ;
    double Cal_MethylationLevel (double tmpDt ) ;
    double Cal_SwitchProbability (double tmpDt ) ;
    void UpdateMotility_FromConfigFile() ;
    
    
};

class bacterium
{   public:
    //---------------------------- Parameters and sub-classes ------------------------------
    MotilityMetabolism motilityMetabolism ;
    vector<node> nodes ;
    //Interact with nodes of other bacteria to avoid volume exclusion
    vector<node> ljnodes ;
    //allnodes is consist of nodes, and ljnodes
    vector<node> allnodes ;
    
    //the duplicates of baceria that are close to the boundary.
    //It is required to handle visualizations and also other interactions
    //There might be multiple duplicated if the bacteria is close to two boundaries
    vector<node> duplicate ;
    //bacteria type IV pili attaching to the substraite or other bacteria depending on species we are modeling
    vector<pilus> pili ;        //Currently not involved
    
    vector<double> connectionToOtherB ;    //Required to tranfer protein between bacteria
    double bhf_Fx = 0;      //Forces to apply impact of Bacterial_HighwayFollowing
    double bhf_Fy = 0;
    
    // Controlling the duration that bacteria stay in the run mode.
    // Changes with chemotaxis. It is also randomly chosen after each reverse
    double reversalPeriod = 7.0 ;   //This value updates during the simulations
    double internalReversalTimer ;  //Keep track of the duration that the bacteria spent on current run mode
    double maxRunDuration ;
    double protein ;
    bool duplicateIsNeeded ;
    double orientation ;
    
    //Holding chemoattractnat source locations.
    //Needed to check if the bacteria is in vicinity of any point sources
    vector<vector<double> > chemoPointSources ;
    
    double locVelocity ;
    double locFriction ;
    
    //Turn angle is randomly chosen from a uniform distribution.
    // To have a normal distribution, defining the distribution is required.
    // Need to calibrate the normal_turnAngle_SDV with experimental data
    //Note: The net turning angle depends on motorForce_Amplitude, turnTimer and chosen angle. Find the correlation
    bool turnStatus = false ;
    double turnTimer = 0.0 ;
    double turnAngle = 0.0 ;        //the overall turn angle depends on motorForce_Amplitude, turnTimer and chosen angle
    double maxTurnAngle = (3.1415)/5.0 ;
    
    double chemotaxisPeriod = 0.1 ;             //How often the bacteria responds to the chemo-gradient
    double chemotaxisTimer = 0.0 ;              // Duration since the last time doing chemotaxis
    
    //oldLocations and oldChem are required to apply chemotaxis using observational model
    vector<double> oldLocation ;
    double oldChem = 0.0 ;
    int numberReverse = 0 ;
    bool sourceWithin = false ;         // Is the bacteria within a source vicinity?
    double sourceVicinity = 5.0 ;
    double timeToSource = 0.0 ;     // The overall duration that bacteria spent to reach the source vicinity
    double timeInSource = 0.0 ;     // the overall duration that bacteria spent in the source vicinity
    
    // We may consider having different velocities for moving forward and backward
    bool directionOfMotion = true ;        //true->forward
    bool attachedToFungi ;          // the bacteria tend to follow the highway if it is attachedToFungi
    
    //---------------------------- Wrapping parameters -----------------------------------------
    //Wrap mode postpone the reverse. Reverse happens when wrap duration is over. ( reverse after unwrap)
    bool wrapMode = false ;
    double wrapPeriod = 1.0 ;   // The duration that bacteria can spend in wrap mode, controlled by chemotaxis
    double wrapTimer = 0.0 ;    //internal timer to keep track of duration spent on current wrap mode
    double wrapAngle = 0.0 ;
    double wrapRate = 1.0 ;     // wrapping frequency, needed for uncalibrated durations
    //double wrapAngularVelocity ;
    double wrapProbability = 1.0 ;      // Probablity of entering to a wrap mode right before reversing
    double wrapSlowDown = 15.0/35.0 ;   // Bacteria moves slower in the wrap mode
    
    // Currently the wrapPeriod follows experimental distribution
    double maxWrapAngle = (3.1415)/6.0 ;    //Needed for uniform distributions
    // Needed for calibration to experimental data.
    
    //---------------------------- Functions -----------------------------------------
    void initialize_RandomForce () ;
    void UpdateBacteria_FromConfigFile() ;
    //Update the maximum durations( run and wrap) based on logNormal distribution
    double LogNormalMaxRunDuration (std::lognormal_distribution<double> &dist, std::default_random_engine &generator , double a, bool calib, double runVal) ;
    bacterium () ;
    bool SourceRegion () ;
};


#endif /* Bacteria_hpp */
