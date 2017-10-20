#include <stdio.h>
#include <math.h>
#include <time.h>
#include <glpk.h>
#include "lp_glpk.h"
#include "util.h"
#include "def.h"

#define _DEA 1
#define _QDEA 1
#define _DEA2 1

output_obj lp_glpk ( dea_obj * object ) {
    int nc, nr, nDMU, nY, nX;
    int *rest;
    float *objvals, *rhs;
    float **A;
    int objtype;

    int nYX, nDMUp1;

    float *objvalsD, *rhsD, **AD;
    double *arD;
    int *restD, *iaD, *jaD, *ijconrowD;
    int nrD, ncD, kD;
    glp_prob * dea;
    glp_smcp parmD;

    clock_t t1, t2, t3;
    FILE *file;

    double *ar;
    int *ia, *ja, *ijconrowQ;
    int kQ;
    glp_prob * qdea;
    glp_smcp parmQ;

    float *rhsD2, **AD2;
    double *arD2;
    int *restD2, *iaD2, *jaD2;
    int nrD2, kD2;
    glp_prob * dea2;
    glp_smcp parmD2;

    int i, j, k, itmp, ret, jDMU, jcount, jtmp, jtmp1, jtmp2, dstart, dend, dcount, ndea2;

    float z8;

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
    output.soln = alloc_init_double_array((size_t)output.nsoln, 0.0);


    double ** EFFscores;
    double * soln;
    double * runtimes;
    EFFscores = output.EFFscores;
    soln = output.soln;
    runtimes = output.runtimes;

    z8 = 0.0;

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

    ret = glp_term_out(GLP_MSG_OFF);

    dea = glp_create_prob();
#if _DEA
    {
        switch (objtype) {
            case LP_MAXIMIZATION: {
                glp_set_obj_dir(dea, GLP_MAX);
            }
                break;
            case LP_MINIZATION: {
                glp_set_obj_dir(dea, GLP_MIN);
            }
                break;
            case LP_OPTIMIZATION: { ;
            }
                break;
            default: { ;
            }
                break;
        }
        ret = glp_add_rows(dea, nrD);

        for (i = 0; i < nrD; i++) {
            switch (rest[i]) {
                case LTE: {
                    glp_set_row_bnds(dea, i + 1, GLP_UP, z8, rhsD[i]);
                }
                    break;
                case GTE: {
                    glp_set_row_bnds(dea, i + 1, GLP_LO, z8, rhsD[i]);
                }
                    break;
                default: { ;
                }
                    break;
            }
        }

        ret = glp_add_cols(dea, ncD);

        for (j = 0; j < ncD; j++) {
            glp_set_obj_coef(dea, j + 1, objvalsD[j]);
            glp_set_col_bnds(dea, j + 1, GLP_LO, z8, z8);
        }

        k = (nDMU + 1) * (nY + nX) + 1;
        arD = alloc_init_double_array((size_t) k, 0.0);
        iaD = alloc_init_int_array((size_t) k, 0);
        jaD = alloc_init_int_array((size_t) k, 0);
        ijconrowD = alloc_init_int_array((size_t) nrD * ncD, 0);

        itmp = 0;
        kD = 0;
        for (i = 0; i < nrD; i++) {
            for (j = 0; j < ncD; j++) {
                if (fabsf(AD[i][j]) > EPS) {
                    kD++;
                    iaD[kD] = i + 1;
                    jaD[kD] = j + 1;
                    arD[kD] = (double) AD[i][j];
                    if (i + 1 == nDMUp1 && j + 1 > nY && j + 1 <= (nY + nX)) {
                        ijconrowD[itmp] = kD;
                        itmp++;
                    }
                }
            }
        }
        ijconrowD = (int *) _resize(ijconrowD, (size_t) itmp, sizeof(int));

        glp_load_matrix(dea, kD, iaD, jaD, arD);
        glp_init_smcp(&parmD);
    }
#endif

    qdea = glp_create_prob();
#if _QDEA
    {
        switch (objtype) {
            case LP_MAXIMIZATION: {
                glp_set_obj_dir(qdea, GLP_MAX);
            }
                break;
            case LP_MINIZATION: {
                glp_set_obj_dir(qdea, GLP_MIN);
            }
                break;
            case LP_OPTIMIZATION: { ;
            }
                break;
            default: { ;
            }
                break;
        }
        ret = glp_add_rows(qdea, nr);

        for (i = 0; i < nr; i++) {
            switch (rest[i]) {
                case LTE: {
                    glp_set_row_bnds(qdea, i + 1, GLP_UP, z8, rhs[i]);
                }
                    break;
                case GTE: {
                    glp_set_row_bnds(qdea, i + 1, GLP_LO, z8, rhs[i]);
                }
                    break;
                default: { ;
                }
                    break;
            }
        }

        ret = glp_add_cols(qdea, nc);

        for (j = 0; j < nc; j++) {
            glp_set_obj_coef(qdea, j + 1, objvals[j]);
            glp_set_col_bnds(qdea, j + 1, GLP_LO, z8, z8);
        }

        k = (nr*nc) + 1;
        ar = alloc_init_double_array((size_t) k, 0.0);
        ia = alloc_init_int_array((size_t) k, 0);
        ja = alloc_init_int_array((size_t) k, 0);
        ijconrowQ = alloc_init_int_array((size_t) nr * nc, 0);

        itmp = 0;
        kQ = 0;
        for (i = 0; i < nr; i++) {
            for (j = 0; j < nc; j++) {
                if (fabsf(A[i][j]) > EPS) {
                    kQ++;
                    ia[kQ] = i + 1;
                    ja[kQ] = j + 1;
                    ar[kQ] = (double) A[i][j];
                    if (i + 1 == nDMUp1 && j + 1 > nY && j + 1 <= (nY + nX)) {
                        ijconrowQ[itmp] = kQ;
                        itmp++;
                    }
                }
            }
        }
        ijconrowQ = (int *) _resize(ijconrowQ, (size_t) itmp, sizeof(int));

        glp_load_matrix(qdea, kQ, ia, ja, ar);
        glp_init_smcp(&parmQ);
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
            glp_set_obj_coef(dea, j + 1, objvalsD[j]);
            glp_set_obj_coef(qdea, j + 1, objvals[j]);
        }

        for ( j = 0; j < nX; j++ ) {
            jtmp= nY + j ;
            itmp = ijconrowD[j];
            arD[itmp] = AD[jDMU][jtmp] * (float)-1.0;
            itmp = ijconrowQ[j];
            ar[itmp] = A[jDMU][jtmp] * (float)-1.0;
        }
        glp_load_matrix(dea, kD, iaD, jaD, arD);
        glp_load_matrix(qdea, kQ, ia, ja, ar);

        ret = glp_simplex(dea, &parmD);
        objD = glp_get_obj_val(dea);

        EFFscores[jDMU][0] = objD;

        ret = glp_simplex(qdea, &parmQ);
        obj = glp_get_obj_val(qdea);

        EFFscores[jDMU][1] = obj;

        for ( j = 0; j < nc; j++ ) {
            soln[j] = glp_get_col_prim(qdea,j+1);
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

        dea2 = glp_create_prob();
#if _DEA2
        {
            switch (objtype) {
                case LP_MAXIMIZATION: {
                    glp_set_obj_dir(dea2, GLP_MAX);
                }
                    break;
                case LP_MINIZATION: {
                    glp_set_obj_dir(dea2, GLP_MIN);
                }
                    break;
                case LP_OPTIMIZATION: { ;
                }
                    break;
                default: { ;
                }
                    break;
            }
            ret = glp_add_rows(dea2, nrD2);

            for (i = 0; i < nrD2; i++) {
                switch (restD2[i]) {
                    case LTE: {
                        glp_set_row_bnds(dea2, i + 1, GLP_UP, z8, rhsD2[i]);
                    }
                        break;
                    case GTE: {
                        glp_set_row_bnds(dea2, i + 1, GLP_LO, z8, rhsD2[i]);
                    }
                        break;
                    default: { ;
                    }
                        break;
                }
            }

            ret = glp_add_cols(dea2, ncD);

            for (j = 0; j < ncD; j++) {
                glp_set_obj_coef(dea2, j + 1, objvalsD[j]);
                glp_set_col_bnds(dea2, j + 1, GLP_LO, z8, z8);
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

                    }
                }
            }

            glp_load_matrix(dea2, kD2, iaD2, jaD2, arD2);
            glp_init_smcp(&parmD2);
        }
#endif
        ret = glp_simplex(dea2, &parmD2);
        objD2 = glp_get_obj_val(dea2);
        EFFscores[jDMU][2] = objD2;

        glp_delete_prob(dea2);

        jcount=jcount+1;

        if (jcount == 50) {
            t3 = clock();
            file = fopen("fortran-report.txt", "a");
            fprintf(file, "jDMU, runtimes = %d %lf %lf", jDMU,(double)(t3 - t2) / CLOCKS_PER_SEC,(double)(t3 - t1) / CLOCKS_PER_SEC);
            fclose(file);
        }
    }

    t3 = clock();

    runtimes[0] = (double)(t3 - t1) / CLOCKS_PER_SEC;
    runtimes[1] = (double)(t3 - t2) / CLOCKS_PER_SEC;

    file = fopen("fortran-report.txt", "a");
    fprintf(file, "jDMU, runtimes = %d %lf %lf\n", jDMU,(double)(t3 - t2) / CLOCKS_PER_SEC,(double)(t3 - t1) / CLOCKS_PER_SEC);
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
    if (iaD) {
        _free(iaD);
    }
    if (jaD) {
        _free(jaD);
    }
    if (ijconrowD) {
        _free(ijconrowD);
    }
    if (ar) {
        _free(ar);
    }
    if (ia) {
        _free(ia);
    }
    if (ja) {
        _free(ja);
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
    glp_delete_prob(dea);
    glp_delete_prob(qdea);
    glp_free_env();

    return output;
}