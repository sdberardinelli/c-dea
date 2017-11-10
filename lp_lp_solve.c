#include <stdio.h>
#include <math.h>
#include <time.h>
#include <lpsolve/lp_lib.h>
#include "lp_lp_solve.h"
#include "def.h"

#define _DEA 1
#define _QDEA 1
#define _DEA2 1

output_obj lp_lp_solve ( dea_obj * object ) {
    int nc, nr, nDMU, nY, nX;
    int *rest;
    float *objvals, *rhs;
    float **A;
    int objtype;

    int nYX, nDMUp1;

    float *objvalsD, *rhsD, **AD;
    double *arD;
    int *restD, *ijconrowD, *iconrowD, *jconrowD;
    int nrD, ncD, kD;
    lprec * dea;

    clock_t t1, t2, t3;
    FILE *file;

    double *ar;
    int *ijconrowQ, *iconrowQ, *jconrowQ;
    int kQ;
    lprec * qdea;

    float *rhsD2, **AD2;
    double *arD2;
    int *restD2, *iaD2, *jaD2;
    int nrD2, kD2;
    lprec * dea2;

    int i, j, k, itmp, ret, jDMU, jcount, jtmp, jtmp1, jtmp2, dstart, dend, dcount, ndea2;

    double obj, objD, objD2;

    double *dmult,*dvals;


    A = object->A;
    nc = object->nobjvals;
    nr = object->ncons;
    nY = object->nY;
    nX = object->nX;
    nDMU = object->nDMU;
    objvals = object->objvals;
    rhs = object->rhs;
    objtype = object->objtype;
    rest = object->rest;

    output_obj output;
    init_output_obj(&output);
    output.nEFFscores = nDMU;
    output.mEFFscores = 3;
    output.EFFscores = alloc_init_double_array2d((size_t)output.nEFFscores, (size_t)output.mEFFscores, 0.0);
    output.nruntimes = 2;
    output.runtimes =  alloc_init_double_array((size_t)output.nruntimes, 0.0);
    output.nsoln = nc;
    output.soln = alloc_init_double_array((size_t)output.nsoln*2, 0.0);


    double ** EFFscores;
    double * soln;
    double * runtimes;
    EFFscores = output.EFFscores;
    soln = output.soln;
    runtimes = output.runtimes;

    nYX = nY + nX;
    nDMUp1 = nDMU+1;
    nrD = nDMUp1;
    ncD = nYX;

    dstart=nYX+2;
    dend=dstart+nDMU-1;

    t1 = clock();

    restD = (int*)copy_array(rest,(size_t)nDMUp1,sizeof(int));
    objvalsD = (float*)copy_array(objvals,(size_t)nYX,sizeof(float));
    rhsD = (float*)copy_array(rhs,(size_t)nDMUp1,sizeof(float));
    AD = (float**)copy_array(A,((size_t)nDMUp1*nYX),sizeof(float));

    dea = make_lp(nrD+1,ncD+1);
    set_verbose(dea, IMPORTANT);
#if _DEA
    {
        switch (objtype) {
            case LP_MAXIMIZATION: {
                set_maxim(dea);
            }
                break;
            case LP_MINIZATION: {
                set_minim(dea);
            }
                break;
            case LP_OPTIMIZATION: { ;
            }
                break;
            default: { ;
            }
                break;
        }

        for (i = 0; i < nrD; i++) {
            switch (restD[i]) {
                case LTE: {
                    ret = set_constr_type(dea,i+1,LE);
                    if ( ret != TRUE ) {
                        printf("dea: set_constr_type failed : %d\n", ret);
                    }
                }
                    break;
                case GTE: {
                    ret = set_constr_type(dea,i+1,GE);
                    if ( ret != TRUE ) {
                        printf("dea: set_constr_type failed : %d\n", ret);
                    }
                }
                    break;
                default: { ;
                }
                    break;
            }
            ret = set_rh(dea,i+1,rhsD[i]);
            if ( ret != TRUE ) {
                printf("dea: set_rh failed : %d\n", ret);
            }
        }

        k = (nDMU + 1) * (nY + nX) + 1;
        arD = alloc_init_double_array((size_t) k, 0.0);
        ijconrowD = alloc_init_int_array((size_t) nrD * ncD, 0);
        iconrowD = alloc_init_int_array((size_t) nrD * ncD, 0);
        jconrowD = alloc_init_int_array((size_t) nrD * ncD, 0);

        itmp = 0;
        kD = 0;
        for (i = 0; i < nrD; i++) {
            for (j = 0; j < ncD; j++) {
                if (fabsf(AD[i][j]) > EPS) {
                    kD++;
                    arD[kD] = (double) AD[i][j];
                    if (i + 1 == nDMUp1 && j + 1 > nY && j + 1 <= (nY + nX)) {
                        ijconrowD[itmp] = kD;
                        iconrowD[itmp] = i+1;
                        jconrowD[itmp] = j+1;
                        itmp++;
                    }
                    ret = set_mat(dea,i+1,j+1,arD[kD]);
                    if ( ret != TRUE ) {
                        printf("dea: set_mat failed : %d\n", ret);
                    }
                }
            }
        }
        ijconrowD = (int *) _resize(ijconrowD, (size_t) itmp, sizeof(int));
        iconrowD = (int *) _resize(iconrowD, (size_t) itmp, sizeof(int));
        jconrowD = (int *) _resize(jconrowD, (size_t) itmp, sizeof(int));

        for ( j = 0; j < ncD; j++ ) {
            ret = set_obj(dea, j+1, objvalsD[j]);
            if ( ret != TRUE ) {
                printf("dea: set_obj failed : %d\n", ret);
            }
        }
    }
#endif

    qdea = make_lp(nr+1,nc+1);
    set_verbose(qdea,IMPORTANT);
#if _QDEA
    {
        switch (objtype) {
            case LP_MAXIMIZATION: {
                set_maxim(qdea);
            }
                break;
            case LP_MINIZATION: {
                set_minim(qdea);
            }
                break;
            case LP_OPTIMIZATION: { ;
            }
                break;
            default: { ;
            }
                break;
        }

        for (i = 0; i < nr; i++) {
            switch (rest[i]) {
                case LTE: {
                    ret = set_constr_type(qdea,i+1,LE);
                    if ( ret != TRUE ) {
                        printf("qdea: set_constr_type failed : %d\n", ret);
                    }
                }
                    break;
                case GTE: {
                    ret = set_constr_type(qdea,i+1,GE);
                    if ( ret != TRUE ) {
                        printf("qdea: set_constr_type failed : %d\n", ret);
                    }
                }
                    break;
                default: { ;
                }
                    break;
            }
            ret = set_rh(qdea,i+1,rhs[i]);
            if ( ret != TRUE ) {
                printf("qdea: set_rh failed : %d\n", ret);
            }
        }

        for ( j = 0; j < nc; j++ ) {
            ret = set_obj(qdea, j+1, objvals[j]);
            if ( ret != TRUE ) {
                printf("qdea: set_obj failed : %d\n", ret);
            }
        }

        k = (nr*nc) + 1;
        ar = alloc_init_double_array((size_t) k, 0.0);
        ijconrowQ = alloc_init_int_array((size_t) nr * nc, 0);
        iconrowQ = alloc_init_int_array((size_t) nr * nc, 0);
        jconrowQ = alloc_init_int_array((size_t) nr * nc, 0);

        itmp = 0;
        kQ = 0;
        for (i = 0; i < nr; i++) {
            for (j = 0; j < nc; j++) {
                if (fabsf(A[i][j]) > EPS) {
                    kQ++;
                    ar[kQ] = (double) A[i][j];
                    if (i + 1 == nDMUp1 && j + 1  > nY && j + 1 <= (nY + nX)) {
                        ijconrowQ[itmp] = kQ;
                        iconrowQ[itmp] = i+1;
                        jconrowQ[itmp] = j+1;
                        itmp++;
                    }
                    ret = set_mat(qdea,i+1,j+1,ar[kQ]);
                    if ( ret != TRUE ) {
                        printf("qdea: set_mat failed : %d\n", ret);
                    }
                }
            }
        }
        ijconrowQ = (int *) _resize(ijconrowQ, (size_t) itmp, sizeof(int));
        iconrowQ = (int *) _resize(iconrowQ, (size_t) itmp, sizeof(int));
        jconrowQ = (int *) _resize(jconrowQ, (size_t) itmp, sizeof(int));
    }
#endif
    jcount = 0;

    dmult = (double*)_alloc((size_t)nDMU, sizeof(double));
    dvals = (double*)_alloc((size_t)nDMU, sizeof(double));

    ndea2 = (nDMU+1)*(nY+nX) + 1;
    AD2 = alloc_init_float_array2d((size_t)nDMU+1,(size_t)nY+nX, 0.0);
    arD2 = alloc_init_double_array((size_t)ndea2, 0.0);
    iaD2 = alloc_init_int_array((size_t)ndea2, 0);
    jaD2 = alloc_init_int_array((size_t)ndea2, 0);
    rhsD2 = alloc_init_float_array((size_t)nDMU+1,0.0);

    restD2 = (int*)_alloc((size_t)nDMU+1, sizeof(int));

    t2 = clock();

    for ( jDMU = 0; jDMU < nDMU; jDMU++ )
    {
        for ( i = 0; i < nY; i++) {
            objvalsD[i] = AD[jDMU][i];
            objvals[i] = A[jDMU][i];
        }

        jtmp1 = nY; //?
        jtmp2 = nY + nX;

        for ( i = jtmp1; i < jtmp2; i++ ) {
            AD[nrD-1][i] = AD[jDMU][i] * (float)-1.0;
        }

        for ( j = 0; j < nY+nX; j++ ) {
            ret = set_obj(dea,j+1,objvalsD[j]);
            if ( ret != TRUE ) {
                printf("dea dmu: set_obj failed : %d\n", ret);
            }
            ret = set_obj(qdea,j+1,objvals[j]);
            if ( ret != TRUE ) {
                printf("dea dmu: set_obj failed : %d\n", ret);
            }
        }

        for ( j = 0; j < nX; j++ ) {
            jtmp= nY + j;
            itmp = ijconrowD[j];
            arD[itmp] = AD[jDMU][jtmp] * (float)-1.0;
            ret = set_mat(dea,iconrowD[j],jconrowD[j],arD[itmp]);
            if ( ret != TRUE ) {
                printf("dea dmu: set_mat failed : %d\n", ret);
            }
        }

        for ( j = 0; j < nX; j++ ) {
            jtmp= nY + j;
            itmp = ijconrowQ[j]+1;
            ar[itmp] = A[jDMU][jtmp] * (float)-1.0;

            ret = set_mat(qdea,iconrowQ[j],jconrowQ[j],ar[itmp]);
            if ( ret != TRUE ) {
                printf("qdea dmu: set_mat failed : %d\n", ret);
            }
        }

        ret = solve(dea);
        if ( ret != OPTIMAL ) {
            printf("dea solve: not OPTIMAL : %d\n", ret);
        }
        objD = get_objective(dea);

        EFFscores[jDMU][0] = objD;

        ret = solve(qdea);
        obj = get_objective(qdea);
        //if ( ret != OPTIMAL ) {
        //    printf("qdea solve: not OPTIMAL : %d\n", ret);
        //}

        EFFscores[jDMU][1] = obj;

        ret = get_variables(qdea,soln);
        if ( ret != TRUE ) {
            printf("qdea dmu: get_variables failed : %d\n", ret);
        }
        _init_double_array(dmult,(size_t)nDMU,0.0);
        _init_double_array(dvals,(size_t)nDMU,0.0);
        for ( i = 0, j = dstart; j < dend; j++, i++ ) {
            dvals[i+1] = soln[j];
        }
        dcount = 0;
        for ( i = 0; i < nDMU; i++ ) {
            if ( dvals[i] > EPS ) {
                dmult[i] = 1.0000000005; //?
            }
            dcount += (int)dmult[i];
        }

        nrD2 = nrD - dcount;

        _init_float_array2d(AD2,(size_t)nDMU+1, (size_t)nY+nX, 0.0);
        _init_float_array(rhsD2,(size_t)nDMU+1, 0.0);
        _init_int_array(restD2, (size_t)nDMU+1, LTE);
        itmp = 0;
        for ( i = 0; i < nrD; i++ ) {
            if ( i < nDMU ) {
                if (dmult[i] < 0.01) {
                    for ( j = 0; j < ncD; j++ ) {
                        AD2[itmp][j] = AD[i][j];
                    }
                    rhsD2[itmp]=rhsD[i];
                    restD2[itmp]=restD[i];
                    itmp++;
                }
            }

            if ( i >= nDMU ) {

                for ( j = 0; j < ncD; j++ ) {
                    AD2[itmp][j] = AD[i][j];
                }
                rhsD2[itmp]=rhsD[i];
                restD2[itmp]=restD[i];
                itmp++;
            }
        }

        dea2 = make_lp(nrD2+1,ncD+1);
        set_verbose(dea2,IMPORTANT);
#if _DEA2
        {
            switch (objtype) {
                case LP_MAXIMIZATION: {
                    set_maxim(dea2);
                }
                    break;
                case LP_MINIZATION: {
                    set_minim(dea2);
                }
                    break;
                case LP_OPTIMIZATION: { ;
                }
                    break;
                default: { ;
                }
                    break;
            }

            for (i = 0; i < nrD2; i++) {
                switch (restD2[i]) {
                    case LTE: {
                        ret = set_constr_type(dea2,i+1,LE);
                        if ( ret != TRUE ) {
                            printf("dea2 dmu: set_constr_type failed : %d\n", ret);
                        }
                    }
                        break;
                    case GTE: {
                        ret = set_constr_type(dea2,i+1,GE);
                        if ( ret != TRUE ) {
                            printf("dea2 dmu: set_constr_type failed : %d\n", ret);
                        }
                    }
                        break;
                    default: { ;
                    }
                        break;
                }
                ret = set_rh(dea2,i+1,rhsD2[i]);
                if ( ret != TRUE ) {
                    printf("dea2 dmu: set_rh failed : %d\n", ret);
                }
            }

            for ( j = 0; j < ncD; j++ ) {
                ret = set_obj(dea2, j+1, objvalsD[j]);
                if ( ret != TRUE ) {
                    printf("dea2 dmu: set_obj failed : %d\n", ret);
                }
            }

            _init_double_array(arD2,(size_t)ndea2, 0.0);
            _init_int_array(iaD2,(size_t)ndea2, 0);
            _init_int_array(jaD2,(size_t)ndea2, 0);
            kD2 = 0;
            for (i = 0; i < nrD2; i++) {
                for (j = 0; j < ncD; j++) {
                    if (fabsf(AD2[i][j]) > EPS) {
                        kD2++;
                        iaD2[kD2] = i + 1;
                        jaD2[kD2] = j + 1;
                        arD2[kD2] = (double) AD2[i][j];
                        ret = set_mat(dea2,i+1,j+1,arD2[kD2]);
                        if ( ret != TRUE ) {
                            printf("dea2 dmu: set_mat failed : %d\n", ret);
                        }
                    }
                }
            }
        }
#endif
        ret = solve(dea2);
        if ( ret != OPTIMAL ) {
            printf("dea2 solve: not OPTIMAL : %d\n", ret);
        }
        objD2=get_objective(dea2);
        EFFscores[jDMU][2] = objD2;

        delete_lp(dea2);

        jcount=jcount+1;

        if (jcount == 50) {
            t3 = clock();
            file = fopen("fortran-report.txt", "a");
            fprintf(file, "lpsolve jDMU, runtimes = %d %lf %lf\n", jDMU,(double)(t3 - t2) / CLOCKS_PER_SEC,(double)(t3 - t1) / CLOCKS_PER_SEC);
            fclose(file);
        }
    }

    t3 = clock();

    runtimes[0] = (double)(t3 - t1) / CLOCKS_PER_SEC;
    runtimes[1] = (double)(t3 - t2) / CLOCKS_PER_SEC;

    file = fopen("fortran-report.txt", "a");
    fprintf(file, "lpsolve jDMU, runtimes = %d %lf %lf\n", jDMU,(double)(t3 - t2) / CLOCKS_PER_SEC,(double)(t3 - t1) / CLOCKS_PER_SEC);
    fclose(file);

    if (restD) {
        _free(restD);
    }
    if (objvalsD) {
        _free(objvalsD);
    }
    if (rhsD) {
        _free(rhsD);
    }
    if (AD) {
        _free(AD);
    }
    if (arD) {
        _free(arD);
    }
    if (ijconrowD) {
        _free(ijconrowD);
    }
    if (ar) {
        _free(ar);
    }
    if (dmult) {
        _free(dmult);
    }
    if (dvals) {
        _free(dvals);
    }
    if (AD2) {
        _free2d((void**)AD2,(size_t)nDMU+1);
    }
    if (arD2) {
        _free(arD2);
    }
    if (iaD2) {
        _free(iaD2);
    }
    if (jaD2) {
        _free(jaD2);
    }
    if (rhsD2) {
        _free(rhsD2);
    }
    if (restD2) {
        _free(restD2);
    }
    if (ijconrowQ) {
        _free(ijconrowQ);
    }
    if (iconrowQ) {
        _free(iconrowQ);
    }
    if (jconrowQ) {
        _free(jconrowQ);
    }
    if (iconrowD) {
        _free(iconrowD);
    }
    if (jconrowD) {
        _free(jconrowD);
    }
    if(dea != NULL) {
        delete_lp(dea);
    }
    if(qdea != NULL) {
        delete_lp(qdea);
    }

    return output;
}