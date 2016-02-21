/*----------------------------------------------------------------------------*/
/* UFPR – Bacharelado em Ciência da Computação                                */
/* CI164: Introducao a Computacao Cientifica, 1sem/2014                       */
/* Primeiro Trabalho Pratico                                                  */
/* Alunos:  Marfan Fragoso Veras Junior - GRR20113754                         */
/*          Leonardo C Robaskievicz     - GRR20117797                         */
/*----------------------------------------------------------------------------*/
#include "cholesky.h"
/*----------------------------------------------------------------------------*/
double *cholesky(double *A, int n) {
    int i, j, k;
    double *L = (double*)calloc(n * n, sizeof(double));
    if (L == NULL)
        fprintf (stderr, "Falta de memória.\n");
 
    for (i = 0; i < n; i++) {
        for (j = 0; j < (i+1); j++) {
            double s = 0;
            for (k = 0; k < j; k++) {
                s += L[i * n + k] * L[j * n + k];
            }
            if (i == j) {
                L[i * n + j] = sqrt(A[i * n + i] - s);
            }
            else {
                L[i * n + j] = (1.0 / L[j * n + j] * (A[i * n + j] - s));
            }
        }
    }
    return L;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double* forward_subs (double *a, int n) {
    int i = 0, j = 0, k;
    double *inv = (double*)calloc(n * n, sizeof(double));
    if (inv == NULL)
        fprintf (stderr, "Falta de memória.\n");
    for (k = 0; k < n; k++) {  //colunas de inv
        for (i = 0; i < n; i++) {
            if (i == k)
                inv[i * n + k] = 1;
            else
                inv[i * n + k] = 0;
            for (j = 0; j < i; j++) {
                inv[i * n + k] -= a[i * n + j] * inv[j * n + k];
            }
            inv[i * n + k] /= a[i * n + i];
        }
    }
    return inv;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double* mult_at_a (double *a, int n) { //retorna a^t*a
    int h, i, j;
    double *p = (double*)calloc(n * n, sizeof(double));
    if (p == NULL)
        fprintf (stderr, "Falta de memória.\n");
    for (h = 0; h < n; h++) {
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                p[h * n + j] += a[i * n + h] * a[i * n + j];
            }
        }
    }       
    return p;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
int Check_Matrix(double *A, int n) {
    int i,j,k; double sum;
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = A[i * n + j];
            for (k = i - 1; k >= 0; k--)
                sum -= A[i * n + k] * A[j * n + k];
            if ( (i == j) && (sum <= 0.0) )
                    return 0;
        }
    }
    return 1;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double *inverse (double *A, int n, double *tempo) {
    double t1 = timestamp();

    double *L = cholesky (A, n);

    double *L_inv = forward_subs (L, n);

    double *A_inv = mult_at_a (L_inv, n);
    
    free(L);
    free(L_inv);

    *tempo = timestamp()-t1;
    return A_inv;
}
/*----------------------------------------------------------------------------*/
double *mult (double *A, double *B, int n) {
    int h, i, j;
    double *p = (double*)calloc(n * n, sizeof(double));
    if (p == NULL)
        fprintf (stderr, "Falta de memória.\n");
    for (h = 0; h < n; h++) {
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                p[h * n + j] += A[h * n + i] * B[i * n + j];
            }
        }
    }
    return p;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double *subtrai_da_identidade(double *A, int n) {
    int i, j;
    double *a = (double*)malloc(n * n * sizeof(double));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j)
                a[i * n + j] = 1 - A[i * n + j];
            else
                a[i * n + j] = 0 - A[i * n + j];
        }
    }
    return a;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double residuo (double *A, double *A_inv, int n, double *tempo_err) {
    double t1 = timestamp();
    double res = 0.0;
    double *aux = mult (A, A_inv, n);
    int i;
    
    double *R = subtrai_da_identidade(aux, n);

    for (i = 0; i < (n*n); i++) {
        res += R[i]*R[i];
    }
    res = sqrt(res);
    
    free (aux);
    free(R);
    *tempo_err = timestamp()-t1;

    return res;
}
/*----------------------------------------------------------------------------*/
