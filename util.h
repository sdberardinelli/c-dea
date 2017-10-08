#ifndef UTIL_H
#define UTIL_H
#include <stdlib.h>

typedef struct {
    int nDMU, nY, nX, nobjvals, ncons, objtype;
    int *rhs_values;
    int *const_signs;
    int *obj_coeffs;
    float** A;
} dea_obj;

void init_dea_obj(dea_obj *);
void free_dea_obj(dea_obj *);
void display_dea_obj(dea_obj *);
char** str_split(char*, char, int*);
void* _alloc(size_t, size_t);
void _free(void *);
void** _alloc2d(size_t, size_t, size_t);
void _free2d(void**, size_t);
void rm_newline(char*, ssize_t);
dea_obj process_file(const char *);

#endif