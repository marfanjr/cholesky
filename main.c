/*----------------------------------------------------------------------------*/
/* UFPR – Bacharelado em Ciência da Computação                                */
/* CI164: Introducao a Computacao Cientifica, 1sem/2014                       */
/* Primeiro Trabalho Pratico                                                  */
/* Alunos:  Marfan Fragoso Veras Junior - GRR20113754                         */
/*          Leonardo C Robaskievicz     - GRR20117797                         */
/*----------------------------------------------------------------------------*/
#include "lib-aux.h"
#include "cholesky.h"
#include <unistd.h>

int main (int argc, char **argv) {
    srand( 20141 );
    
    int c, //guarda a opçao (-i, -o, -r)
        n = -1; //guarta a dimensao da matriz (n x n)
        
    char *arq_out; //string p/ guardar o caminho do arquivo de saida
    
    FILE *input = stdin,
         *output = stdout;
    
    double *A = NULL,
            tempo_inv,
            tempo_err;
    
    while ((c = getopt(argc, argv, "i:o:r:")) != -1) {
        switch (c) {
            case 'i':
                A = le_entrada (optarg, &n);
            break;
            case 'o':
                output = fopen (optarg, "w");
            break;
            case 'r':
                n = atoi(optarg);
            break;
            default:
                //err = 1;
            break;
        }
    }
    if ( (A == NULL) && (n == -1) )  // se n foi dado um arquivo como entrada
        A = le_entrada ("stdin", &n);// le da entrada padrao

    if ( (A == NULL) && (n != -1) )  // se foi dado apenas a dimensao da matriz
        A = generateSquareRandomPositiveDefiniteMatrix( n );

//    if (!Check_Matrix(A, n)) {
//        fprintf (stderr, "Matriz não definida positiva.\n");
//        return -1;
//    }

 
    double *A_inv = inverse (A, n, &tempo_inv);
    double res = residuo (A, A_inv, n, &tempo_err);
    
    imprime_saida (A_inv, n, tempo_inv, tempo_err, res, output);
    free(A);
    free(A_inv);
    return 0;
}
