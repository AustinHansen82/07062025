
#include "Bacteria.hpp"

bacterium::bacterium ()
{
    nodes.resize(nnode) ;
    ljnodes.resize(nnode-1) ;
    allnodes.resize(2*nnode-1) ;
    duplicate.resize(nnode) ;
    connectionToOtherB.resize(nbacteria) ;
    oldLocation.resize(2) ;
    
    pili.resize(nPili) ;
}
bool bacterium::SourceRegion()
{
    int numberOfPoints = 0 ;
    bool tmpInsource = false ;
    if (chemoPointSources.size() != 0)
    {
        numberOfPoints = chemoPointSources.at(0).size() ;
    }
    for (int i=0; i< numberOfPoints; i++)
    {
        double pointX = chemoPointSources.at(0).at(i) ;
        double pointY = chemoPointSources.at(1).at(i) ;
        double r2 = sqrt( pow( (nodes[(nnode-1)/2].x - pointX ) ,2 ) + pow( (nodes[(nnode-1)/2].y - pointY) ,2 ) );
        if (r2< sourceVicinity )
        {
            tmpInsource = true ;
            
        }
    }
    //else is not written yet. once it is in source sourceWithin remains true forever
    
    return tmpInsource ;
}

void bacterium::UpdateBacteria_FromConfigFile()
{
    motilityMetabolism.UpdateMotility_FromConfigFile() ;
    
    maxTurnAngle = globalConfigVars.getConfigValue("Bacteria_maxTurnAngle").toDouble() ;
    chemotaxisPeriod = globalConfigVars.getConfigValue("Bacteria_chemotaxisPeriod").toDouble() ;
    wrapPeriod = globalConfigVars.getConfigValue("wrap_duration").toDouble() ;
    wrapProbability = globalConfigVars.getConfigValue("wrap_probability").toDouble() ;
    wrapSlowDown = globalConfigVars.getConfigValue("wrap_slowDown").toDouble() ;
    maxWrapAngle = globalConfigVars.getConfigValue("wrap_maxAngle").toDouble() ;
    maxRunDuration = 1.0 / globalConfigVars.getConfigValue("Bacteria_reversalRate").toDouble() ;
    
}

//-------------------------------MotilityMetabolism----------------------------------------------------

MotilityMetabolism::MotilityMetabolism ()
{
}

double MotilityMetabolism::Cal_MethylationEnergy()
{
    return MethylEnergy = km * (m0 - methylation ) ;
}

double MotilityMetabolism::Cal_LegandEnergy()
{
    return LegandEnergy = log(1.0 + legand/ kI) - log(1.0 + legand/ kA) ;
}

double MotilityMetabolism::Cal_ReceptorActivity()
{
    double totalEnergy = Cal_LegandEnergy() + Cal_MethylationEnergy() ;
    return receptorActivity = 1.0/ ( 1.0 + exp(N * totalEnergy) ) ;
}

double MotilityMetabolism::Cal_MethylationLevel(double tmpDt)
{
    double changeRate = kR * (1.0 - receptorActivity) - kB * receptorActivity ;
    return methylation += tmpDt * changeRate ;
    
}
double MotilityMetabolism::Cal_SwitchProbability(double tmpDt)
{
    
    methylation = Cal_MethylationLevel(tmpDt) ;
    receptorActivity = Cal_ReceptorActivity() ;
    return switchProbability = exp(lnA_a + b * (receptorActivity/ (receptorActivity + gamma ) ) ) ;
}

void bacterium::initialize_RandomForce()
{
    for (int j=0 ; j<nnode ; j++)
    {
        nodes[j].xdev = 0.0 ;
        nodes[j].ydev = 0.0 ;
        
    }
}

void MotilityMetabolism::UpdateMotility_FromConfigFile()
{
    a0 = globalConfigVars.getConfigValue("motility_initialActivity").toDouble() ;
    m0 = globalConfigVars.getConfigValue("motility_initialMethylation").toDouble() ;
    lnA_a = globalConfigVars.getConfigValue("motility_Normalization").toDouble() ;
    b = globalConfigVars.getConfigValue("motility_ScalingFactor").toDouble() ;
    gamma = globalConfigVars.getConfigValue("motility_HalfOccupied").toDouble() ;
    N = globalConfigVars.getConfigValue("motility_receptorNumber").toDouble() ;
    kI = globalConfigVars.getConfigValue("motility_Feedback_kI").toDouble() ;
    kA = globalConfigVars.getConfigValue("motility_Feedback_kA").toDouble() ;
    kR = globalConfigVars.getConfigValue("motility_Feedback_kR").toDouble() ;
    kB = globalConfigVars.getConfigValue("motility_Feedback_kB").toDouble() ;
    km = globalConfigVars.getConfigValue("motility_methylEnergyCoeff_km").toDouble() ;
    
    
}

double bacterium::LogNormalMaxRunDuration(std::lognormal_distribution<double> &dist, std::default_random_engine &generator , double a, bool calib, double runVal)
{
 
    if (calib == true)
    {
        
    maxRunDuration = 1.0/ a * dist(generator) ;
    return maxRunDuration ;
    }
    else
    {
        maxRunDuration = runVal ;
        return maxRunDuration ;
    }
}
