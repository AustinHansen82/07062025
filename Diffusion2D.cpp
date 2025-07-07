#include "Diffusion2D.hpp"



TissueGrid Diffusion2D (double xMin, double xMax, double yMin, double yMax,int nGridX , int nGridY ,vector<vector<double> > sources, vector<double> pSrc, TissueGrid tissue, Chemo_Profile_Type profileType)
{
    vector<vector<double> > pointSource = sources ;
    tissue.DomainBoundaries(xMin, xMax, yMin, yMax, nGridX , nGridY) ;
    
    tissue.xSources = pointSource.at(0) ;
    tissue.ySources = pointSource.at(1) ;
    //Create the 2D grid structure
    tissue.InitializeAllGrids() ;
    //Assign the production rates corresponding to the sources
    tissue.FindProductionPoints(pSrc) ;
    // Let the chemical secrete and diffuse over time
    // Run for a fixed number of iterations unless the profile achive steady-state condition
    if (profileType == production_profile)
    {
        tissue.EulerMethod() ;
    }
    else if (profileType == linear_profile)
    {
        tissue.Create_Linear_Gradient();
        tissue.ParaViewGrids(1);
    }
    else if (profileType == experimental_profile)
    {
        tissue.Create_Experimental_Gradient();
        tissue.ParaViewGrids(1);
    }

    
    return tissue ;
    
    
}
