//
//  driver.h
//  Myxobacteria
//
//  Created by Alireza Ramezani on 11/12/20.
//  Copyright © 2020 Alireza Ramezani. All rights reserved.
//

#ifndef driver_h
#define driver_h


#include "Fungi.hpp"


Fungi driver (Fungi fungi) ;
Fungi Manual_FungalNetwork ( Fungi fungi) ;    //This is just a template. Modify according to your plan
//Load the fungal network based on coordinates of hyphae segments.
//The coordinates are found from simulating the fungal growth
Fungi Load_FungalNetwork ( Fungi fungi) ;


#endif /* driver_h */
