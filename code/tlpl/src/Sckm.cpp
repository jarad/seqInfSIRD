/* 
Sckm class definition. 
*/

#include "Sckm.h"

using namespace std;

Sckm::Sckm() 
{
    allocatedPre    = false;
    allocatedPost   = false;
    allocatedStoich = false;
    allocatedRates  = false;
    allocatedState  = false;

    s = 0;
    r = 0;
};

Sckm::Sckm(unsigned int _s, unsigned int _r)
{
    s = _s;
    r = _r;
}

Sckm::Sckm(unsigned int _s, unsigned int _r, 
           unsigned int *_pre, unsigned int *_post, int *_stoich, 
           double *_rates, unsigned int *_state)
{
    s = _s;
    r = _r;

    unsigned int *pre = malloc(s*r * sizeof(unsigned int));
    if (NULL == pre) 
    {
        error("pre could not be allocated")
    } else 
    {
        allocatedPre = true;
    }
}


Sckm::~Sckm() 
{
  if (allocatedPre   ) delete [] pre;
  if (allocatedPost  ) delete [] post;
  if (allocatedState ) delete [] state;
  if (allocatedStoich) delete [] stoich;
  if (allocatedRates ) delete [] rates;
}


int main()
{
}

