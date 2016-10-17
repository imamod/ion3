#include <cstdio>
#include "src/elements.h"
#include "src/sahasolver.h"
#include <cmath>

int main()
{
    SahaSolver c(elements::Cu);
    //TElement elem(29);
    //SahaSolver c(elem);

    for(double lgT = 4.6; lgT >= 0.499; lgT -= 0.05)
    {
        for(double lgV = -3; lgV <= 6.01; lgV += 0.05)
        {
            SahaPoint res = c.Calculate_lgTeV_lgVae(lgT,lgV);
            //printf("lgT = %g lgV = %g xe = %g P = %g E = %g Mu = %g\n",lgT, lgV, res.Xe,pow(10,res.lgP), pow(10,res.lgE), res.M);
            printf("%.3f ",res.Xe);
        }
        printf("\n");
    }
    printf("\n");

    return 0;
}
