

#ifndef Grid_hpp
#define Grid_hpp

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <vector>
#include <numeric>
#include <stdio.h>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include "driver.h"
using namespace std ;

class Grid
{
public:
    double value = 0.0 ;
    double change = 0.0 ;
    double productionRate = 0.0 ;
    
    
    void UpdateValue () ;
    
};
#endif /* Grid_hpp */
