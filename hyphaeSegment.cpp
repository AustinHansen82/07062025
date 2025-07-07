#include "hyphaeSegment.h"




HyphaeSegment:: HyphaeSegment()
{
    hyphaeLength = globalConfigVars.getConfigValue("hyphae_length").toDouble() ;
}

