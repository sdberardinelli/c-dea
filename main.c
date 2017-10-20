#include <stdlib.h>
#include <stdio.h>
#include "lp_glpk.h"


int main ( int argc, char * argv[] ) {
    dea_obj dea;
    output_obj output;

    dea = process_file("../examples/GLPKdata.txt");

    output = lp_glpk(&dea);

    printf("Runtime 1 = %lf\n", output.runtimes[0]);
    printf("Runtime 2 = %lf\n", output.runtimes[1]);

    write_output_obj(&output,"EFFscores.csv");

    free_dea_obj(&dea);
    free_output_obj(&output);

    exit(EXIT_SUCCESS);
}