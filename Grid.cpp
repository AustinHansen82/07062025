#include "Grid.hpp"



void Grid::UpdateValue()
{
    value += change ;
    change = 0.0 ;
    
    return ;
}
