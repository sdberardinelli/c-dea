#include <stdlib.h>
#include <stdio.h>
#include "util.h"
#include "lp_glpk.h"
#include "lp_lp_solve.h"

#define _GLPK 1
#define _LPSOLVE 1

int main ( int argc, char * argv[] ) {
    dea_obj dea;
    output_obj output;

    dea = process_file("../examples/GLPKdata-CCR-example.txt");

    //display_dea_obj(&dea);

#if _GLPK
    output = lp_glpk(&dea);
    printf("GLPK: Runtime 1 = %lf\n", output.runtimes[0]);
    printf("GLPK: Runtime 2 = %lf\n", output.runtimes[1]);
    write_output_obj(&output,"glpk-EFFscores.csv");
    free_output_obj(&output);
#endif

#if _LPSOLVE
    output = lp_lp_solve(&dea);
    printf("LP_SOLVE: Runtime 1 = %lf\n", output.runtimes[0]);
    printf("LP_SOLVE: Runtime 2 = %lf\n", output.runtimes[1]);
    write_output_obj(&output,"lpsolve-EFFscores.csv");
    free_output_obj(&output);
#endif

    free_dea_obj(&dea);

    exit(EXIT_SUCCESS);
}