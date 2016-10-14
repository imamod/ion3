#include <cstdio>
#include "src/elements.h"
#include "src/sahasolver.h"

int main()
{
    SahaSolver c(elements::Cu);
    c.calculate_lgTeV_lgVae(1,3);

    printf("xe = %g",c.xe());
    return 0;
}
