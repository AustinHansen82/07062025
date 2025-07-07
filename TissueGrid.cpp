
#include "TissueGrid.hpp"

TissueGrid::TissueGrid ()
{
    grid_dx = static_cast<double> ( ( xDomainMax - xDomainMin ) / numberGridsX ) ;
    grid_dy = static_cast<double> ( ( yDomainMax - yDomainMin ) / numberGridsY ) ;
}

void TissueGrid::DomainBoundaries(double xMin, double xMax, double yMin, double yMax, int nGridX , int nGridY)
{
    xDomainMin = xMin ;
    xDomainMax = xMax ;
    yDomainMin = yMin ;
    yDomainMax = yMax ;
    numberGridsX = nGridX ;
    numberGridsY = nGridY ;
    grid_dx = static_cast<double> ( ( xDomainMax - xDomainMin ) / numberGridsX ) ;
    grid_dy = static_cast<double> ( ( yDomainMax - yDomainMin ) / numberGridsY ) ;
    grid_dt *= 1.0 /Diffusion ;
}



void TissueGrid::InitializeAllGrids()
{
    vector<vector <Grid> > tmp (numberGridsY,vector<Grid> ( numberGridsX) ) ;
    grids = tmp ;
    
    return ;
}

void TissueGrid::ClearChanges()
{
    for (int i =0; i < numberGridsY; i++)
    {
        for (int j = 0; j< numberGridsX ; j++)
        {
            
            grids.at(i).at(j).change = 0.0 ;
            
        }
    }
    return ;
}

void TissueGrid::DiffusionChanges()
{
    double tmpChange = 0.0 ;
    
    for (int i =0; i < numberGridsY; i++)
    {
        for (int j = 0; j< numberGridsX -1 ; j++)
        {
            tmpChange = Diffusion * ( grids.at(i).at(j+1).value - grids.at(i).at(j).value) * grid_dt / (grid_dx * grid_dx) ;
            grids.at(i).at(j).change += tmpChange ;
            grids.at(i).at(j+1).change += -tmpChange ;

        }
    }
    for (int i =0; i < numberGridsY -1; i++)
    {
        for (int j = 0; j< numberGridsX ; j++)
        {
            tmpChange = Diffusion * (grids.at(i+1).at(j).value - grids.at(i).at(j).value) * grid_dt / (grid_dy * grid_dy) ;
            grids.at(i).at(j).change += tmpChange ;
            grids.at(i+1).at(j).change += -tmpChange ;
            
        }
    }
    return ;
}

void TissueGrid::FindProductionPoints(vector<double> pSrc)
{
    int tmpIndexX ;
    int tmpIndexY ;
    for (unsigned int i = 0; i < xSources.size(); i++)
    {
        tmpIndexX = static_cast<int>(round ( ( xSources.at(i) - xDomainMin ) / grid_dx ) ) ;
        tmpIndexX = fmod(tmpIndexX, numberGridsX ) ;
        tmpIndexY = static_cast<int>(round ( ( ySources.at(i) - yDomainMin ) / grid_dy ) ) ;
        tmpIndexY = fmod(tmpIndexY, numberGridsY ) ;
        
        //Based on relative distance in segment
        grids.at(tmpIndexY).at(tmpIndexX).productionRate = pSrc.at(i) ;
        
        indexSourceX.push_back(tmpIndexX) ;
        indexSourceY.push_back(tmpIndexY) ;
        
    }
    return ;
}

void TissueGrid::ProductionChanges()
{
    double tmpChange ;
    int tmpIdX ;
    int tmpIdY ;
    for (unsigned int i =0; i < indexSourceX.size() ; i++)
    {
        tmpIdX = indexSourceX.at(i) ;
        tmpIdY = indexSourceY.at(i) ;
        
        tmpChange = ( grids.at(tmpIdY).at(tmpIdX).productionRate) * grid_dt ;
        grids.at(tmpIdY).at(tmpIdX).change += tmpChange ;
        
    }
    return ;
}

void TissueGrid::DegredationChanges()
{
    double tmpChange ;
    for (int i =0; i < numberGridsY ; i++)
    {
        for (int j = 0; j< numberGridsX ; j++)
        {
            tmpChange = - deg * (grids.at(i).at(j).value ) * grid_dt ;
            grids.at(i).at(j).change += tmpChange ;
            
        }
    }
}

void TissueGrid::UpdateChanges()
{
    for (int i =0; i < numberGridsY ; i++)
    {
        for (int j = 0; j< numberGridsX ; j++)
        {
            grids.at(i).at(j).UpdateValue() ;
        }
    }
    return ;
}

void TissueGrid::EulerMethod()
{
    int l = 0 ;
    bool status = false ;
    while (status == false && l < maxIterator) //10000
    {
        status = true ;
        double smallValue = 0.0001 ;
        ClearChanges() ;
        DiffusionChanges() ;
        //DegredationChanges() ;
        ProductionChanges() ;
        
        //Check for steady-state condition
        for (int i=0; i< numberGridsY ; i++)
        {
            for (int j = 0; j< numberGridsX ; j++)
            {
                if (grids.at(i).at(j).change/( (grids.at(i).at(j).value + smallValue) ) > smallValue * grid_dt  )
                {
                    status = false ;
                    break ;
                    
                }
            }
        }
        UpdateChanges() ;
        if (l%100 == 0)
        {
            //ParaViewGrids(l/100) ;
        }
        l++ ;
        
    }
    RoundToZero() ;
    ParaViewGrids(l/100) ;
    cout<<l<<endl ;
    return ;
}
void TissueGrid::Create_Linear_Gradient()
{
 
    double cntrX = numberGridsX/2.0 ;
    double cntrY = numberGridsY/2.0 ;
    for (int i = 0; i<numberGridsY; i++)
    {
        for (int j=0; j< numberGridsX; j++)
        {
            //grids.at(j).at(i).value = grad_scale*(i); // 1D - Linear Function
            //grids.at(j).at(i).value = 180*exp((1500-sqrt(pow((j-cntrX),2.0)+pow((i-cntrY),2.0)))/200); // MWC Profile
            grids.at(j).at(i).value = grad_scale*(sqrt(pow((cntrX),2.0)+pow((cntrY),2.0))-sqrt(pow((j-cntrX),2.0)+pow((i-cntrY),2.0))); // Radial Linear Function
            // grids.at(j).at(i).value = 1; // Constant Function
        }
    }

}

void TissueGrid::Create_Experimental_Gradient()
{
    //std::ifstream file("external_concentrations4_444444444444445_new.txt");
    //std::ifstream file("external_concentrations40.888888888888886_new.txt");
    //std::ifstream file("external_concentrations40.888888888888886_newflipupd.txt");
    //std::ifstream file("external_concentrations40.888888888888886_newtrans.txt");
    //std::ifstream file("external_concentrations40.888888888888886_fliptransup.txt");
    //std::ifstream file("external_concentrations50.67_new.txt");
    //std::ifstream file("external_concentrations50.67_newtrans.txt");
    //std::ifstream file("external_concentrations50.67_fliptransup.txt");
    //std::ifstream file("external_concentrations8.88888888888889_new.txt");
    //std::ifstream file("external_concentrations120.88888888888889_new.txt");
    //std::ifstream file("external_concentrations156.44444444444446_new.txt");
    //std::ifstream file("external_concentrations_t=176000.00_new.txt");
    //std::ifstream file("external_concentrations_t=176000.00_newflipupd.txt");
    //std::ifstream file("/Users/austinhansen/Documents/Research/Bacterial Migration/Pseudomonas_Motility_SCE_Model/external_concentrations_t=160000.00_newtrans2.txt");
    //std::ifstream file("/Users/austinhansen/Documents/Research/Bacterial Migration/Pseudomonas_Motility_SCE_Model/external_concentrations_t=36032.00_newtrans2_Penic_Everywhere.txt");
    //std::ifstream file("/Users/austinhansen/Documents/Research/Bacterial Migration/Pseudomonas_Motility_SCE_Model/external_concentrations_t=36032.00_newtrans2_Penic_Everywhere_4000Grid2.txt");
    //std::ifstream file("/Users/austinhansen/Documents/Research/Bacterial Migration/Pseudomonas_Motility_SCE_Model/external_concentrations_t=36032.00_newtrans2_Penic_Tips.txt");
    //std::ifstream file("/Users/austinhansen/Documents/Research/Bacterial Migration/Pseudomonas_Motility_SCE_Model/external_concentrations_t=36032.00_newtrans2_Penic_Tips_4000Grid2.txt");
    std::ifstream file("/rhome/ahans016/bigdata/Bacterial_Migration/Pseudomonas_Motility_SCE_Model_Redo_5DEF/Diffusion_Point_Source_2000Micron.txt");
    //std::ifstream file("/rhome/ahans016/bigdata/Bacterial_Migration/Pseudomonas_Motility_SCE_Model/external_concentrations_t=36032.00_newtrans2_Penic_Everywhere.txt");
    //std::ifstream file("/rhome/ahans016/bigdata/Bacterial_Migration/Pseudomonas_Motility_SCE_Model/external_concentrations_t=36032.00_newtrans2_Penic_Everywhere_4000Grid2.txt");
    //std::ifstream file("/rhome/ahans016/bigdata/Bacterial_Migration/Pseudomonas_Motility_SCE_Model/external_concentrations_t=36032.00_newtrans2_Penic_Tips.txt");
    //std::ifstream file("/rhome/ahans016/bigdata/Bacterial_Migration/Pseudomonas_Motility_SCE_Model/external_concentrations_t=36032.00_newtrans2_Penic_Tips_4000Grid2.txt");
    //std::ifstream file("/rhome/ahans016/bigdata/Bacterial_Migration/Fungi_In_Liquid_Simulations/external_concentrations_t=160000.00_newtrans2.txt");
    //std::ifstream file("external_concentrations_t=176000.00_fliptransup.txt");

        std::vector<std::vector<double>> External_Concentration;
        std::string line;
        while (std::getline(file, line))
        {
            std::vector<double> row;
            std::stringstream ss(line);
            std::string value;
            while (std::getline(ss, value, ','))
            {
                row.push_back(std::stod(value));
            }
            External_Concentration.push_back(row);
        }

            for (int i = 0; i<numberGridsY; i++)
            {
                for (int j=0; j< numberGridsX; j++)
                {
                    grids.at(j).at(i).value = External_Concentration[j][i];

                }
            }
}

void TissueGrid::ParaViewGrids(int index)
{
    string vtkFileName2 = folderName + "GridChem"+ to_string(index)+ ".vtk" ;
    ofstream SignalOut;
    SignalOut.open(vtkFileName2.c_str());
    SignalOut << "# vtk DataFile Version 2.0" << endl;
    SignalOut << "Result for paraview 2d code" << endl;
    SignalOut << "ASCII" << endl;
    SignalOut << "DATASET RECTILINEAR_GRID" << endl;
    SignalOut << "DIMENSIONS" << " " << numberGridsX  << " " << " " << numberGridsY << " " << 1  << endl;
    
    SignalOut << "X_COORDINATES " << numberGridsX << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int i = 0; i < numberGridsX ; i++) {
        SignalOut << i * grid_dx << endl;
    }
    
    SignalOut << "Y_COORDINATES " << numberGridsY << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int j = 0; j < numberGridsY; j++) {
        SignalOut << j * grid_dy << endl;
    }
    
    SignalOut << "Z_COORDINATES " << 1 << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int k = 0; k < 1 ; k++) {
        SignalOut << 0 << endl;
    }
    
    SignalOut << "POINT_DATA " << (numberGridsX )*( numberGridsY ) << endl;
    SignalOut << "SCALARS trehalose float 1" << endl;
    SignalOut << "LOOKUP_TABLE default" << endl;
    
    for (int k = 0; k < 1 ; k++) {
        for (int j = 0; j < numberGridsY; j++) {
            for (int i = 0; i < numberGridsX ; i++) {
                SignalOut << grids.at(j).at(i).value << endl;
            }
        }
    }
    

}
void TissueGrid::RoundToZero()
{
    for (int i = 0; i<numberGridsY; i++)
    {
        for (int j=0; j< numberGridsX; j++)
        {
            if (grids.at(j).at(i).value < pow(10, -30) )
            {
                grids.at(j).at(i).value = 0.0 ;
            }
        }
    }
    // Code to Save Chemoattractant Profile for Further Use
    /*
    std::string filename = "Chemo_Profile.csv";  // Added .csv extension for clarity
    std::ofstream outFile(filename);
        
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
            return;
        }

        // Set precision for floating-point output
        outFile << std::setprecision(15);

        for (int i = 0; i < numberGridsY; i++)
        {
            for (int j = 0; j < numberGridsX; j++)
            {
                outFile << grids.at(j).at(i).value;
                
                // Add a comma separator if it's not the last element in the row
                if (j < numberGridsX - 1) {
                    outFile << ",";
                }
            }
            outFile << "\n";  // New line after each row
        }

        outFile.close();

        std::cout << "Grid data saved to " << filename << std::endl;
    */
}

void TissueGrid::UpdateTGrid_FromConfigFile()
{
    folderName = globalConfigVars.getConfigValue("AnimationFolder").toString() ;
    statsFolder = globalConfigVars.getConfigValue("StatFolderName").toString() ;
    numberGridsX = globalConfigVars.getConfigValue("grid_NumberX").toDouble() ;
    numberGridsY = globalConfigVars.getConfigValue("grid_NumberY").toDouble() ;
    Diffusion = globalConfigVars.getConfigValue("grid_DiffusionCoeff").toDouble() ;
    deg = globalConfigVars.getConfigValue("grid_degradationRate").toDouble() ;
    pro = globalConfigVars.getConfigValue("grid_productionRate").toDouble() ;
    grid_dt = globalConfigVars.getConfigValue("grid_timeStep").toDouble() ;
    grad_scale = globalConfigVars.getConfigValue("grad_scale").toDouble() ;
    maxIterator = globalConfigVars.getConfigValue("grid_maxIterator").toDouble() ;
    chemo_profile_type =static_cast<Chemo_Profile_Type>( globalConfigVars.getConfigValue("chemo_profile_type").toInt() ) ;

}
