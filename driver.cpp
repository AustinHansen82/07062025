
#include "driver.h"
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
using constants::pi;


// ---------------------------------------------------
//   Main Driver Function
// ---------------------------------------------------
Fungi driver (Fungi fungi){

    // ---------------------------------------------------
    //   Set Up Static Mycelia
    // ---------------------------------------------------
    // Initialize starting parameters
    fungi.UpdateFungi_FromConfigFile() ;
    if (fungi.loading_Network == true)
    {
        fungi = Load_FungalNetwork(fungi) ;
    }
    else
    {
        fungi = Manual_FungalNetwork(fungi) ;
    }
    
    if (fungi.branchIsTip == true)
    {
        fungi.Find_Hyphae_Tips() ;
    }
    else
    {
        fungi.Find_Hyphae_Tips2() ;
    }
    
    fungi.FindCenterPointConnection() ;
    //fungi.FindProductionNetwork() ;

    return fungi ;
    
}
//-----------------------------------------------------------------------------------------------------

Fungi Manual_FungalNetwork ( Fungi fungi)
{
    // NOTE:
    // The i^th hyphae segment is considered a tip if
    //      hy[i].can_extend == true
    // The location of the tip is given by
    //      hy[i].x2 & hy[i].y2

    // Set up initial mycelia
    HyphaeSegment new_hy;
    cout << " " << endl;
    cout << "INITIAL HYPHAE:" << endl;
    cout << "---------------" << endl;
    cout << "initX" << fungi.initX << "initY" << fungi.initY << "hyphaeLength" << new_hy.hyphaeLength << endl;
    double initAngle =  pi ; // pi / 4.0 ;
    for ( int i = 0; i < fungi.init_Count; i++ ) {
        
        // Starting points and angles
        
        new_hy.x1 = fungi.initX + 200.0 ;            //60                                          // One endpoint of i^th segment: x-coordinate
        new_hy.y1 = fungi.initY + 200.0;             //60                                         // One endpoint of i^th segment: y-coordinate
        new_hy.angle = i * ( ( 2 * pi ) / fungi.init_Count ) + initAngle ;  // Angle of i^th segment
        new_hy.x2 = new_hy.x1 + new_hy.hyphaeLength * cos( new_hy.angle );                      // Other endpoint of i^th segment: x-coordinate
        new_hy.y2 = new_hy.y1 + new_hy.hyphaeLength * sin( new_hy.angle );                      // Other endpoint of i^th segment: y-coordinate
        new_hy.can_branch = true;                                            // The i^th segment can extend
        new_hy.can_extend = true;                                            // The i^th segment can branch
        new_hy.from_who = -1 ;
        new_hy.extend_to = -1 ;
        new_hy.branch_to = -1 ;
        
        fungi.hyphaeSegments.push_back( new_hy );                             // Append hy vector for new segment

        // Print out info:
        cout << "Segment " << i << endl;
        cout << "   Endpoints of new segment:       (" << fungi.hyphaeSegments[i].x1 << ", " << fungi.hyphaeSegments[i].y1 << ") to (" << fungi.hyphaeSegments[i].x2 << ", " << fungi.hyphaeSegments[i].y2 << ")" << endl;
        cout << "   Angle of new segment (radians): " << fungi.hyphaeSegments[i].angle << endl;
    }
    
    for ( int i = 0; i < fungi.init_Count ; i++ ) {

        // Starting points and angles

        new_hy.x1 = fungi.initX ;            //60    // One endpoint of i^th segment: x-coordinate
        new_hy.y1 = fungi.initY;             //60    // One endpoint of i^th segment: y-coordinate
        new_hy.angle = i * ( ( 2 * pi ) / fungi.init_Count ) + initAngle ;  // Angle of i^th segment
        new_hy.x2 = new_hy.x1 + new_hy.hyphaeLength * cos( new_hy.angle );                      // Other endpoint of i^th segment: x-coordinate
        new_hy.y2 = new_hy.y1 + new_hy.hyphaeLength * sin( new_hy.angle );  // Other endpoint of i^th segment: y-coordinate
        new_hy.can_branch = true;                                            // The i^th segment can extend
        new_hy.can_extend = true;                                            // The i^th segment can branch
        new_hy.from_who = -1 ;
        new_hy.extend_to = -1 ;
        new_hy.branch_to = -1 ;
        
        fungi.hyphaeSegments.push_back( new_hy );                             // Append hy vector for new segment

        // Print out info:
        cout << "Segment " << i << endl;
        cout << "   Endpoints of new segment:       (" << fungi.hyphaeSegments[i].x1 << ", " << fungi.hyphaeSegments[i].y1 << ") to (" << fungi.hyphaeSegments[i].x2 << ", " << fungi.hyphaeSegments[i].y2 << ")" << endl;
        cout << "   Angle of new segment (radians): " << fungi.hyphaeSegments[i].angle << endl;
    }

    // Define extra segments (make the shape non-uniform/non-circular)
    cout << " " << endl;
    cout << "EXTENDED/BRANCHED HYPHAE:" << endl;
    cout << "-------------------------" << endl;
    //add_new_hyphae( fungi.hyphaeSegments, 0, fungi.extendingAngle, "extend");      //4th
    //add_new_hyphae( fungi.hyphaeSegments, 1, fungi.extendingAngle, "extend");      //4th
    //add_new_hyphae( fungi.hyphaeSegments, 1, ( pi / 3.0 ), "branch");             //7th
    
    add_new_hyphae( fungi.hyphaeSegments, 0, fungi.extendingAngle, "extend");      //4th
    add_new_hyphae( fungi.hyphaeSegments, 2, 0, "extend");                      //5th
    add_new_hyphae( fungi.hyphaeSegments, 5, 0, "extend");                      //6th
    add_new_hyphae( fungi.hyphaeSegments, 5, ( pi / 3.0 ), "branch");             //7th
    
    add_new_hyphae(fungi.hyphaeSegments, 1, 0, "extend") ;                   //8th
    add_new_hyphae(fungi.hyphaeSegments, 1, pi/3.0, "branch") ;                  //9th
    add_new_hyphae(fungi.hyphaeSegments, 9, 0, "extend") ;                      //10th
   // add_new_hyphae(fungi.hyphaeSegments, 10, 0, "extend") ;                      //11th
    
    add_new_hyphae(fungi.hyphaeSegments, 3, -pi/3.0, "branch") ;                  //12th
    
   // add_new_hyphae( fungi.hyphaeSegments, 7, 0, "extend");                      //13th
    
    add_new_hyphae( fungi.hyphaeSegments, 2, -pi/3.0 , "branch");                      //14th
    add_new_hyphae( fungi.hyphaeSegments, 12, 0, "extend");                      //15th
   // add_new_hyphae( fungi.hyphaeSegments, 15, pi/3.0 , "branch");                      //16th
   // add_new_hyphae( fungi.hyphaeSegments, 6, 0, "extend");                      //17th

    cout << " " << endl;

    add_new_hyphae( fungi.hyphaeSegments, 0, 0, "extend");
    add_new_hyphae( fungi.hyphaeSegments, 1, 0, "extend");
    add_new_hyphae( fungi.hyphaeSegments, 2, 0, "extend");
    add_new_hyphae( fungi.hyphaeSegments, 3, 0, "extend");
    
    return fungi ;

}
//-----------------------------------------------------------------------------------------------------

Fungi Load_FungalNetwork ( Fungi fungi)
{
    // open input file
    
    //std::ifstream inputFile("hyphal_coordinates40.888888888888886.txt");
    std::ifstream inputFile("test2_t=36032.00_hyphal_coordinates_run0_Everywhere.txt");
    //std::ifstream inputFile("test2_t=36032.00_hyphal_coordinates_run0_Tips.txt");
    
    // read each line of the file
    std::string line;
    int hyphae_counter = 0;
    HyphaeSegment new_hy;
    while (std::getline(inputFile, line))
    {
        // use stringstream to split line into doubles
        std::stringstream ss(line);
        std::string value;
        int token_counter = 0;
        while (std::getline(ss, value, ','))
        {
            if (token_counter == 0)
            {
                new_hy.x1 = std::stod(value);
                new_hy.x1 = (new_hy.x1) + fungi.initX;
            }
            else if (token_counter == 1)
            {
                new_hy.y1 = std::stod(value);
                new_hy.y1 = (new_hy.y1) + fungi.initY ;
            }
            else if (token_counter == 2)
            {
                new_hy.x2 = std::stod(value);
                new_hy.x2 = (new_hy.x2) + fungi.initX;
            }
            else if (token_counter == 3)
            {
                new_hy.y2 = std::stod(value);
                new_hy.y2 = (new_hy.y2) + fungi.initY;
            }
            token_counter++;
        }
            
        if (hyphae_counter < 4)
        {
            new_hy.can_branch = true;                                            // The i^th segment can extend
            new_hy.can_extend = true;                                            // The i^th segment can branch
            new_hy.from_who = -1 ;
            new_hy.extend_to = -1 ;
            new_hy.branch_to = -1 ;
            new_hy.hyphae_ID = hyphae_counter ;
            fungi.hyphaeSegments.push_back( new_hy );
        }
        else
        {
            for (int k = 0; k < hyphae_counter; k++)
            {
                if (new_hy.x1 == fungi.hyphaeSegments[k].x2 && new_hy.y1 == fungi.hyphaeSegments[k].y2)
                {
                    new_hy.can_branch = true ;                                            // The i^th segment can extend
                    new_hy.can_extend = true ;                                            // The i^th segment can branch
                    new_hy.from_who = k ;
                    new_hy.extend_to = -1 ;
                    new_hy.branch_to = -1 ;
                    new_hy.hyphae_ID = hyphae_counter ;
                    fungi.hyphaeSegments[k].can_extend = false ;
                    fungi.hyphaeSegments[k].extend_to = hyphae_counter ;
                }
                else if (std::abs(new_hy.x1 - (fungi.hyphaeSegments[k].x1 + fungi.hyphaeSegments[k].x2)/2.0)<.001 && std::abs(new_hy.y1 - (fungi.hyphaeSegments[k].y1 + fungi.hyphaeSegments[k].y2)/2.0) < .001)
                {
                    new_hy.can_branch = true ;                                            // The i^th segment can extend
                    new_hy.can_extend = true ;                                            // The i^th segment can branch
                    new_hy.from_who = k ;
                    new_hy.extend_to = -1 ;
                    new_hy.branch_to = -1 ;
                    new_hy.hyphae_ID = hyphae_counter ;
                    fungi.hyphaeSegments[k].can_branch = false ;
                    fungi.hyphaeSegments[k].branch_to = hyphae_counter ;
                }
            }
            fungi.hyphaeSegments.push_back( new_hy );
        }
        hyphae_counter ++;
    }

    // close input file
    inputFile.close();
    
    for (const auto& obj : fungi.hyphaeSegments)
    {
        std::cout << "x1 is -- " << obj.x1 << "x2 is -- " << obj.x2 << std::endl;
        std::cout << "y1 is -- " << obj.y1 << "y2 is -- " << obj.y2 << std::endl;
        std::cout << "can_branch is --  " << obj.can_branch << std::endl;
        std::cout << "can_extend is -- " << obj.can_extend << std::endl;
        std::cout << "from_who is -- " << obj.from_who << std::endl;
        std::cout << "extend_to is -- " << obj.extend_to << std::endl;
        std::cout << "branch_to is -- " << obj.branch_to << std::endl;
        std::cout << "hyphae_ID is -- " << obj.hyphae_ID << std::endl;
    }
    
    return fungi ;
}
//-----------------------------------------------------------------------------------------------------

