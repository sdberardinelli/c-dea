#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "def.h"

char** str_split( char* str, char delim, int* num_splits )
{
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

void* _alloc ( size_t n, size_t t)
{
    return malloc(n*t);
}

void _free( void * a )
{
    if (a)
    {
        free(a);
    }
}

void** _alloc2d ( size_t n, size_t m, size_t t)
{
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
    obj->rhs_values = 0;
    obj->const_signs = NULL;
    obj->obj_coeffs = NULL;
    obj->A = NULL;
}

void free_dea_obj( dea_obj * obj ) {
    if (obj->rhs_values) {
        _free(obj->rhs_values);
    }
    if ( obj->const_signs ) {
        _free(obj->const_signs);
    }
    if ( obj->obj_coeffs ) {
        _free(obj->obj_coeffs);
    }
    if ( obj->A ) {
        _free2d((void**)obj->A,(size_t)obj->ncons);
    }
}

void display_dea_obj( dea_obj * obj ){

    int i, j;
    printf("-----------\n");
    printf("%d\n", obj->objtype);
    printf("%d %d %d\n", obj->nDMU, obj->nY, obj->nX);
    printf("%d\n", obj->nobjvals);
    if (obj->obj_coeffs) {
        for (i = 0; i < obj->nobjvals; i++) {
            printf("%d ", obj->obj_coeffs[i]);
        }
        printf("\n");
    }
    printf("%d\n", obj->ncons);
    if (obj->const_signs) {
        for (i = 0; i < obj->nobjvals; i++) {
            switch ( obj->const_signs[i] )
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
    if (obj->rhs_values) {
        for (i = 0; i < obj->ncons; i++) {
            printf("%d ", obj->rhs_values[i]);
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
                    obj.obj_coeffs = (int *)_alloc((size_t)obj.nobjvals, sizeof(int));
                    for ( j = 0; j < i; j++ ) {
                        obj.obj_coeffs[j] = (int) strtol(arr[j], NULL, 10);
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
                if ( i == obj.nobjvals ) {
                    obj.const_signs = (int *)_alloc((size_t)obj.nobjvals, sizeof(int));
                    for ( j = 0; j < i; j++ ) {
                        if (strcmp(arr[j], "\"<=\"") == 0) {
                            obj.const_signs[j] = LTE;
                        }
                        else if (strcmp(arr[j], "\">=\"") == 0) {
                            obj.const_signs[j] = GTE;
                        }
                        else if (strcmp(arr[j], "\">\"") == 0) {
                            obj.const_signs[j] = GT;
                        }
                        else if (strcmp(arr[j], "\"<\"") == 0) {
                            obj.const_signs[j] = LT;
                        }
                        else if (strcmp(arr[j], "\"==\"") == 0) {
                            obj.const_signs[j] = EQ;
                        }
                    }
                }
                _free(arr);
            }
        }
        if (strcmp(line, "\"rhs values\"") == 0) {
            read = getline(&line, &len, fp);
            rm_newline(line, read);
            arr = str_split(line,',',&i);
            if (i == obj.ncons) {
                obj.rhs_values = (int *) _alloc((size_t)i, sizeof(int));
                for (j = 0; j < i; j++) {
                    obj.rhs_values[j] = (int) strtol(arr[j], (char **)NULL, 10);
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

void test_ (void )
{
    int i, j;
    int ** arr;

    arr = (int**)_alloc2d(10,10,sizeof(int));

    for ( i = 0; i < 10; i++ ) {
        for ( j = 0; j < 10; j++ ) {
            arr[i][j] = i+1+j;
        }
    }

    for ( i = 0; i < 10; i++ ) {
        for ( j = 0; j < 10; j++ ) {
            printf("%d ", arr[i][j]);
        }
        printf("\n");
    }

    _free2d((void**)arr,10);


    int * array;

    array = (int*)_alloc(10,sizeof(int));

    for ( i = 0; i < 10; i++ )
    {
        array[i] = i+1;
    }
    _free(array);
}