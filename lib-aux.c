/*----------------------------------------------------------------------------*/
/* UFPR – Bacharelado em Ciência da Computação                                */
/* CI164: Introducao a Computacao Cientifica, 1sem/2014                       */
/* Primeiro Trabalho Pratico                                                  */
/* Alunos:  Marfan Fragoso Veras Junior - GRR20113754                         */
/*          Leonardo C Robaskievicz     - GRR20117797                         */
/*----------------------------------------------------------------------------*/
#include "lib-aux.h"

/*----------------------------------------------------------------------------*/
double *generateSquareRandomPositiveDefiniteMatrix( unsigned int n )
{
  double *mat = NULL;

  /* return NULL if memory allocation fails */
  if ( ! (mat = (double *) malloc(n*n*sizeof(double))) )
    return (NULL);

  /* generate a randomly initialized matrix in row-major order */
  double *ptr = mat;
  double *end = mat + n*n;

  double invRandMax = 1.0/(double)RAND_MAX;

  while( ptr != end ) {
    *ptr++ = (double)rand() * invRandMax;
  }

  /* Now we want to make this matrix positive definite. Since all
     values are positive and <1.0, we have to ensure it is symmetric
     and diagonally dominant.

     A = A + transpose(A)
     A = A + I*n                        */
  unsigned int i,j;
  for (i=0; i<n; ++i)
    for (j=i+1; j<n; ++j)
    {
      double aux = mat[i*n+j];
      mat[i*n+j] += mat[j*n+i];
      mat[j*n+i] += aux;
    }

  for (i=0; i<n; ++i)
    mat[i*n+i] += mat[i*n+i] + n;

  return (mat);
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void imprime (double *matriz, unsigned int n) {
    int i;
    for (i = 0; i < n*n; i++) {
        printf("%g  ", matriz[i]);
        //printf("%lf  ", matriz[i]);
        if ((i + 1) % n == 0)
            printf("\n");
    }
    printf("\n");
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double* le_entrada (char *arquivo, int *n) {
    FILE *input;
    if (strcmp (arquivo, "stdin") != 0) 
       input = fopen (arquivo, "r");
    else
        input = stdin;
    int i = 0;
    if (input == NULL) {
        fprintf (stderr, "%s\nNome de arquivo invalido\n", arquivo);
        return NULL;
    }
    else {
        fscanf(input, "%d", n);
        double *A = (double*) malloc ((*n) * (*n) * sizeof(double));
        if (A == NULL)
            fprintf (stderr, "Falta de memória.\n");
        while ( !feof(input) && (i < (*n) * (*n)) ) {
            fscanf(input, "%lf", &(A[i]));
            i++;
        }
        if (input != stdin)
            fclose(input);
        return A;
    }
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void imprime_saida (double *matriz, unsigned int n, double tempo_inv,
                     double tempo_err, double residuo, FILE* output) {
    int i;
    fprintf(output, "#\n");
    fprintf(output, "# Erro: %.17g\n", residuo);
    fprintf(output, "# Tempo Invert: %.17g\n", tempo_inv);
    fprintf(output, "# Tempo Erro: %.17g\n", tempo_err);
    fprintf(output, "#\n");
    fprintf(output, "%d\n", n);
    for (i = 0; i < n*n; i++) {
        fprintf(output, "%.17g  ", matriz[i]);
        if ((i + 1) % n == 0)
            fprintf(output, "\n");
    }
    if (output != stdout)
        fclose(output);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double timestamp(void){
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec + tp.tv_usec/1000000.0));
}
/*----------------------------------------------------------------------------*/
