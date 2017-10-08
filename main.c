#include <stdlib.h>
#include "lp_glpk.h"
#include "util.h"


int main ( int argc, char * argv[] ) {
    dea_obj obj;

    obj = process_file("../examples/GLPKdata-CCR-example.txt");

    lp_glpk(&obj);

    free_dea_obj(&obj);

    exit(EXIT_SUCCESS);
}