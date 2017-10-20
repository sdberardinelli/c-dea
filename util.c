#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "def.h"

char** str_split( char* str, char delim, int* num_splits ) {
    char** ret;
    int ret_len;
    char* c;

    if ( ( str == NULL ) ||
         ( delim == '\0' ) )
    {
        /* Either of those will cause problems */
        ret = NULL;
        ret_len = -1;
    }
    else
    {
        ret_len = 0;
        c = str;

        /* Pre-calculate number of elements */
        do
        {
            if ( *c == delim )
            {
                ret_len++;
            }

            c++;
        } while ( *c != '\0' );

        ret = malloc( ( ret_len + 1 ) * sizeof( *ret ) );
        ret[ret_len] = NULL;

        c = str;
        ret_len = 1;
        ret[0] = str;

        do
        {
            if ( *c == delim )
            {
                ret[ret_len++] = &c[1];
                *c = '\0';
            }

            c++;
        } while ( *c != '\0' );
    }

    if ( num_splits != NULL )
    {
        *num_splits = ret_len;
    }

    return ret;
}

void* _alloc ( size_t n, size_t t) {
    return malloc(n*t);
}

int * alloc_init_int_array(size_t n, int value) {
    int * ptr;
    ptr = (int*)_alloc(n,sizeof(int));
    _init_int_array(ptr,n,value);
    return ptr;
}

int ** alloc_init_int_array2d(size_t n, size_t m, int value) {
    int ** ptr;
    ptr = (int**)_alloc2d(n,m,sizeof(int));
    _init_int_array2d(ptr,n,m,value);
    return ptr;
}

float * alloc_init_float_array(size_t n, float value) {
    float * ptr;
    ptr = (float*)_alloc(n,sizeof(float));
    _init_float_array(ptr,n,value);
    return ptr;
}

float ** alloc_init_float_array2d(size_t n, size_t m, float value) {
    float ** ptr;
    ptr = (float**)_alloc2d(n,m,sizeof(float));
    _init_float_array2d(ptr,n,m,value);
    return ptr;
}

double * alloc_init_double_array(size_t n, double value) {
    double * ptr;
    ptr = (double*)_alloc(n,sizeof(double));
    _init_double_array(ptr,n,value);
    return ptr;
}

double ** alloc_init_double_array2d(size_t n, size_t m, double value) {
    double ** ptr;
    ptr = (double**)_alloc2d(n,m,sizeof(double));
    _init_double_array2d(ptr,n,m,value);
    return ptr;
}

void* _resize(void *a, size_t n, size_t t) {
    void * ptr;
    ptr = _alloc(n,t);
    memcpy(ptr,a,n*t);
    _free(a);
    return ptr;
}

void _init_int_array(int *a, size_t n, int value) {
    int i;
    for ( i = 0; i < n; i++ ) {
        a[i] = value;
    }
}

void _init_float_array(float *a, size_t n, float value) {
    int i;
    for ( i = 0; i < n; i++ ) {
        a[i] = value;
    }
}

void _init_double_array(double *a, size_t n, double value) {
    int i;
    for ( i = 0; i < n; i++ ) {
        a[i] = value;
    }
}

void _init_float_array2d(float **a, size_t n, size_t m, float value) {
    int i;
    for ( i = 0; i < n; i++ ) {
        _init_float_array(a[i],m,value);
    }
}

void _init_double_array2d(double **a, size_t n, size_t m, double value) {
    int i;
    for ( i = 0; i < n; i++ ) {
        _init_double_array(a[i],m,value);
    }
}

void _init_int_array2d(int **a, size_t n, size_t m, int value) {
    int i;
    for ( i = 0; i < n; i++ ) {
        _init_int_array(a[i],m,value);
    }
}

void _free( void * a ) {
    if (a)
    {
        free(a);
    }
}

void** _alloc2d ( size_t n, size_t m, size_t t) {
    int i;
    void ** arr;
    arr = malloc(n*m*t);
    for ( i = 0; i < n; i++ ) {
        arr[i] = _alloc(m,t);
    }
    return arr;
}

void _free2d ( void ** a, size_t n) {
    int i;
    for ( i = 0; i < n; i++ ) {
        _free(a[i]);
    }
    _free((void*)a);
}

void rm_newline(char* line, ssize_t len) {
    if (*line && line[len-1] == '\n'){
        line[len-1] = '\0';
    }
}

void init_dea_obj( dea_obj * obj ) {
    obj->nDMU = 0;
    obj->nY = 0;
    obj->nX = 0;
    obj->nobjvals = 0;
    obj->ncons = 0;
    obj->objtype = 0;
    obj->rhs = 0;
    obj->rest = NULL;
    obj->objvals = NULL;
    obj->A = NULL;
}

void free_dea_obj( dea_obj * obj ) {
    if (obj->rhs) {
        _free(obj->rhs);
    }
    if ( obj->rest ) {
        _free(obj->rest);
    }
    if ( obj->objvals ) {
        _free(obj->objvals);
    }
    if ( obj->A ) {
        _free2d((void**)obj->A,(size_t)obj->ncons);
    }
}

void init_output_obj( output_obj * obj ) {
    obj->nruntimes = 0;
    obj->runtimes = NULL;
    obj->nsoln = 0;
    obj->soln = NULL;
    obj->nEFFscores = 0;
    obj->mEFFscores = 0;
    obj->EFFscores = NULL;
}

void free_output_obj( output_obj * obj ) {
    if (obj->runtimes) {
        _free(obj->runtimes);
    }
    if ( obj->soln ) {
        _free(obj->soln);
    }
    if ( obj->EFFscores ) {
        _free2d((void**)obj->EFFscores,(size_t)obj->nEFFscores);
    }
}

void display_dea_obj( dea_obj * obj ) {

    int i, j;
    printf("-----------\n");
    printf("%d\n", obj->objtype);
    printf("%d %d %d\n", obj->nDMU, obj->nY, obj->nX);
    printf("%d\n", obj->nobjvals);
    if (obj->objvals) {
        for (i = 0; i < obj->nobjvals; i++) {
            printf("%f ", obj->objvals[i]);
        }
        printf("\n");
    }
    printf("%d\n", obj->ncons);
    if (obj->rest) {
        for (i = 0; i < obj->nobjvals; i++) {
            switch ( obj->rest[i] )
            {
                case LTE: { printf("<="); } break;
                case GTE: { printf(">="); } break;
                case LT: { printf("<"); } break;
                case GT: { printf(">"); } break;
                case EQ: { printf("=="); } break;
                default: { printf("??"); } break;
            }
            printf(" ");
        }
        printf("\n");
    }
    if (obj->rhs) {
        for (i = 0; i < obj->ncons; i++) {
            printf("%f ", obj->rhs[i]);
        }
        printf("\n");
    }
    if (obj->A) {
        for (i = 0; i < obj->ncons; i++) {
            for (j = 0; j < obj->nobjvals; j++) {
                printf("%f ", obj->A[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

void display_output_obj( output_obj * obj ) {
    int i, j;
    printf("-----------\n");
    printf("%lf\n", obj->obj);
    printf("%d\n", obj->nruntimes);
    if (obj->runtimes) {
        for (i = 0; i < obj->nruntimes; i++) {
            printf("%lf ", obj->runtimes[i]);
        }
        printf("\n");
    }
    printf("%d\n", obj->nsoln);
    if (obj->soln) {
        for (i = 0; i < obj->nsoln; i++) {
            printf("%lf ", obj->soln[i]);
        }
        printf("\n");
    }

    printf("%d %d\n", obj->nEFFscores, obj->mEFFscores);
    if (obj->EFFscores) {
        for (i = 0; i < obj->nEFFscores; i++) {
            for (j = 0; j < obj->mEFFscores; j++) {
                printf("%f ", obj->EFFscores[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

void write_output_obj( output_obj * obj, const char * filename ) {
    FILE * file;

    file = fopen(filename,"w");

    fprintf(file, "DEA,QDEA1,QDEA2,\n");

    int i, j;
    if (obj->EFFscores) {
        for (i = 0; i < obj->nEFFscores; i++) {
            for (j = 0; j < obj->mEFFscores; j++) {
                fprintf(file, "%lf,",obj->EFFscores[i][j]);
            }
            fprintf(file, "\n");
        }
    }

    fclose(file);
}

void* copy_array(void* ar,size_t n, size_t t) {
    void * ptr;
    ptr = _alloc(n,t);
    memcpy(ptr,ar,n*t);
    return ptr;
}

dea_obj process_file(const char* filename) {
    int i, j, k;
    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        exit(EXIT_FAILURE);
    }

    dea_obj obj;

    init_dea_obj(&obj);

    char **arr;

    arr = NULL;

    while ((read = getline(&line, &len, fp)) != -1) {
        rm_newline(line, read);
        if (strcmp(line, "\"LPlist.objtype\"") == 0) {
            read = getline(&line, &len, fp);
            rm_newline(line, read);
            if (strcmp(line, "\"max\"") == 0) {
                obj.objtype = 1;
            } else if (strcmp(line, "\"min\"") == 0) {
                obj.objtype = -1;
            } else {
                obj.objtype = 0;
            }
        }
        if (strcmp(line, "\"nDMU\",\"nY\",\"nX\"") == 0) {
            read = getline(&line, &len, fp);
            rm_newline(line, read);
            arr = str_split(line, ',', &i);
            if (i == 3) {
                obj.nDMU = (int) strtol(arr[0], NULL, 10);
                obj.nY = (int) strtol(arr[1], NULL, 10);
                obj.nX = (int) strtol(arr[2], NULL, 10);
            }
            _free(arr);
        }
        if (strcmp(line, "\"nobjvals\"") == 0) {
            read = getline(&line, &len, fp);
            rm_newline(line, read);
            obj.nobjvals = (int) strtol(line, NULL, 10);
        }
        if (strcmp(line, "\"obj coeffs\"") == 0) {
            read = getline(&line, &len, fp);
            rm_newline(line, read);
            if (obj.nobjvals) {
                arr = str_split(line,',',&i);
                if ( i == obj.nobjvals ) {
                    obj.objvals = (float *)_alloc((size_t)obj.nobjvals, sizeof(float));
                    for ( j = 0; j < i; j++ ) {
                        obj.objvals[j] = strtof(arr[j], NULL);
                    }
                }
                _free(arr);
            }
        }
        if (strcmp(line, "\"ncons\"") == 0) {
            read = getline(&line, &len, fp);
            rm_newline(line, read);
            obj.ncons = (int) strtol(line, NULL, 10);
        }
        if (strcmp(line, "\"const signs\"") == 0) {
            read = getline(&line, &len, fp);
            rm_newline(line, read);
            if (obj.nobjvals) {
                arr = str_split(line,',',&i);
                if (i == obj.nobjvals || i == obj.ncons) {
                    obj.rest = (int *)_alloc((size_t)i, sizeof(int));
                    for ( j = 0; j < i; j++ ) {
                        if (strcmp(arr[j], "\"<=\"") == 0) {
                            obj.rest[j] = LTE;
                        }
                        else if (strcmp(arr[j], "\">=\"") == 0) {
                            obj.rest[j] = GTE;
                        }
                        else if (strcmp(arr[j], "\">\"") == 0) {
                            obj.rest[j] = GT;
                        }
                        else if (strcmp(arr[j], "\"<\"") == 0) {
                            obj.rest[j] = LT;
                        }
                        else if (strcmp(arr[j], "\"==\"") == 0) {
                            obj.rest[j] = EQ;
                        }
                    }
                }
                else {
                    printf("---> ERROR: ");
                    printf("\"nobjvals\" did not match number of \"const signs\"\n");
                    printf("nobjvals = %d\n\nncos = %d\nfound %d\n",obj.nobjvals,obj.ncons,i);
                    exit(-1);
                }
                _free(arr);
            }
        }
        if (strcmp(line, "\"rhs values\"") == 0) {
            read = getline(&line, &len, fp);
            rm_newline(line, read);
            arr = str_split(line,',',&i);
            if (i == obj.ncons) {
                obj.rhs = (float *) _alloc((size_t)i, sizeof(float));
                for (j = 0; j < i; j++) {
                    obj.rhs[j] = strtof(arr[j], (char **)NULL);
                }
            }
            _free(arr);
        }
        if (strcmp(line, "\"A coefficients\"") == 0) {
            obj.A = (float**)_alloc2d((size_t)obj.ncons,(size_t)obj.nobjvals,sizeof(float));
            for ( i = 0; i < obj.ncons; i++ )
            {
                read = getline(&line, &len, fp);
                if ( read == -1 ){
                    exit(EXIT_FAILURE);
                }
                rm_newline(line, read);
                arr = str_split(line,',',&k);
                if ( k == obj.nobjvals ) {
                    for ( j = 0; j < k; j++ ) {
                        obj.A[i][j] = strtof (arr[j], NULL);
                    }
                }
                _free(arr);
            }
        }
    }

    fclose(fp);
    if (line) {
        _free(line);
    }

    return obj;
}
