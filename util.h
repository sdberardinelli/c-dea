#ifndef UTIL_H
#define UTIL_H
#include <stdlib.h>

typedef struct {
    int nDMU, nY, nX, nobjvals, ncons, objtype;
    float *rhs;//int *rhs_values;
    int *rest;//int *const_signs;
    float *objvals; //int *obj_coeffs;
    float **A;
} dea_obj;

typedef struct {
    double obj;
    int nruntimes;
    double *runtimes;
    int nsoln;
    double *soln;
    int nEFFscores;
    int mEFFscores;
    double **EFFscores;
} output_obj;


void init_dea_obj(dea_obj *);
void init_output_obj(output_obj *);
void free_dea_obj(dea_obj *);
void free_output_obj(output_obj *);
void display_dea_obj(dea_obj *);
void display_output_obj(output_obj *);
void write_output_obj( output_obj *, const char * );

char** str_split(char*, char, int*);
void* _alloc(size_t, size_t);
void _free(void *);
void** _alloc2d(size_t, size_t, size_t);
void _free2d(void**, size_t);
void rm_newline(char*, ssize_t);
dea_obj process_file(const char *);
void* copy_array(void*,size_t,size_t);
void _init_int_array2d(int **, size_t, size_t, int);
void _init_float_array2d(float **, size_t, size_t, float);
void _init_double_array2d(double **, size_t, size_t, double);
void _init_int_array(int *, size_t, int);
void _init_float_array(float *, size_t, float);
void _init_double_array(double *, size_t, double);

int * alloc_init_int_array(size_t, int);
int ** alloc_init_int_array2d(size_t, size_t, int);
float * alloc_init_float_array(size_t, float);
float ** alloc_init_float_array2d(size_t, size_t, float);
double * alloc_init_double_array(size_t, double);
double ** alloc_init_double_array2d(size_t, size_t, double);

void* _resize(void*, size_t, size_t);

#endif