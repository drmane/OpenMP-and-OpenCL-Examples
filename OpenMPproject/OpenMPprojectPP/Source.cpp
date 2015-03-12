#define MAX(A,B) (((A) > (B)) ? (A) : (B)) //declarar la operación de maximo
/*
PRACTICA OPENMP PP
##################
Version Final
Daniel Ruiz Manero DNI: 71299479P

Version 0.1
*/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP //definir las librerias correspondientes a open 
#include <omp.h>
#endif

/*
El programa resuelve el problema de la ecuación 2-D de Laplace Mediante métodos numéricos.
Concretando, la resolución se realiza mediante el método de relajación.
Se ha de discretizar la ecuación correspondiente.
Se implementó el programa sin la utilización de OpenMP y tras ello, se añadieron las directivas y pragmas.
*/

void laplace2DNormal();
void laplace2DOpenMp();


int main() {

	int datos_leidos, opcion;

	do{

		printf("\nIntroduce la opcion deseada: \n");
		printf("############################ \n\n");
		printf("1. Resolucion Ecuacion Laplace 2D Normal: \n");
		printf("2. Resolucion Ecuacion Laplace 2D OpenMP:  \n");
		printf("0. Salir del Programa:  \n");

		datos_leidos = scanf("%u" , &opcion);

		if (datos_leidos!=1){
				  fprintf(stderr, "ERROR!! Introduzca un solo valor!\n");
		}

		if(opcion == 1){
			printf("Se ha elegido la opcion 1 \n\n");
			laplace2DNormal();
		}
		else if(opcion == 2){
			printf("Se ha elegido la opcion 2 \n\n");
			laplace2DOpenMp(); 
		}
	}while(opcion != 0);


   return 0;
}

void laplace2DNormal(){
	double *T, *Tnew, *Tmp;
   double tolerancia, maximo_valor_double = DBL_MAX, cota_superior = 100.0; //la tolerancia es el error admisible

   //DBL_MAX es el mayor valor admisible por un double
   unsigned tamano_grilla, n2, iteraciones_maximas, iteraciones = 0;
   unsigned i, j; //para recorrer los bucles for

   int datos_leidos;

   printf("Introduce el tamaño de la grilla, numero maximo de iteraciones y tolerancia: ");

   //pide tres variables de entrada
   datos_leidos = scanf("%u" , &tamano_grilla);
   datos_leidos += scanf("%u", &iteraciones_maximas);
   datos_leidos += scanf("%lf",&tolerancia);

   if (datos_leidos!=3) { //si no son tres datos leidos da error
      fprintf(stderr, "ERROR!! fallo en la entrada de datos!\n");
      exit(-1);
   }

   time_t tiempo_inicio = clock();

   n2 = tamano_grilla + 2;

   //inicializa los arrays de dimension n2*n2
   T = (double *)calloc(n2*n2, sizeof(*T));
   Tnew = (double *) calloc(n2*n2, sizeof(*T));

   if (T == NULL || Tnew == NULL) {
      fprintf(stderr, "ERROR!! no hay memoria suficiente\n");
      exit(EXIT_FAILURE);
   }

   // establecer las condiciones limite

   for (i = 1; i <= tamano_grilla; i++) {

      T[(tamano_grilla+1)*n2+i] = Tnew[(tamano_grilla+1)*n2+i] = i * cota_superior / (tamano_grilla + 1);

      T[i*n2+tamano_grilla+1] = Tnew[i*n2+tamano_grilla+1] = i * cota_superior / (tamano_grilla + 1);

   }

   while(maximo_valor_double > tolerancia && iteraciones <= iteraciones_maximas) { 

      ++iteraciones;

      maximo_valor_double = 0.0;

      for (i=1; i<=tamano_grilla; ++i) {

         for (j=1; j<=tamano_grilla; ++j) {   

            Tnew[i*n2+j] = 0.25*( T[(i-1)*n2+j] + T[(i+1)*n2+j] + T[i*n2+(j-1)] + T[i*n2+(j+1)] );

			//convergencia
                                
            maximo_valor_double = MAX(maximo_valor_double, fabs(Tnew[i*n2+j] - T[i*n2+j]));
         }
      }
    
    Tmp = T; 
	T = Tnew; 
	Tnew=Tmp;

    if (iteraciones % 100 == 0){
         printf("Iteracion: %8u, variacion = %12.4lE\n", iteraciones, maximo_valor_double);
	 }
   }

   double tiempo_final = (clock() - tiempo_inicio) / (double) CLOCKS_PER_SEC; 
	//calcula el tiempo restando el final menos el del principio y dividiendolo entre los clocks por segundo   

   //muestra el resultado por pantalla

   printf("Tiempo Transcurrido (s) = %.2lf\n", tiempo_final);
   printf("Tamano malla: %u\n", tamano_grilla);
   printf("Parado en la iteracion: %u\n", iteraciones);
   printf("Error maximo: %lE\n", maximo_valor_double);

   free(T);
   free(Tnew);

   //para poder visualizar los datos por pantalla

   system("PAUSE");

}

/**
VERSION OPENMP
##############
Modificación 1: Cambio en el comienzo del reloj
Modificación 2: Paralelizacion del bucle for
Modificación 3: Cambio en el final del reloj


**/
void laplace2DOpenMp(){
	double *T, *Tnew, *Tmp;
   double tolerancia, maximo_valor_double = DBL_MAX, cota_superior = 100.0; //la tolerancia es el error admisible

   //DBL_MAX es el mayor valor admisible por un double
   unsigned tamano_grilla, n2, iteraciones_maximas, iteraciones = 0;
   unsigned i, j; //para recorrer los bucles for

   int datos_leidos;
   FILE *fichero_salida;

   printf("Introduce el tamaño de la grilla, numero maximo de iteraciones y tolerancia: ");

   //pide tres variables de entrada
   datos_leidos = scanf("%u" , &tamano_grilla);
   datos_leidos += scanf("%u", &iteraciones_maximas);
   datos_leidos += scanf("%lf",&tolerancia);

   if (datos_leidos!=3) { //si no son tres datos leidos da error
      fprintf(stderr, "ERROR!! fallo en la entrada de datos!\n");
      exit(-1);
   }


   //Modificación 1
	#ifdef _OPENMP 
		double tiempo_inicio = omp_get_wtime();
	#else
		time_t tiempo_inicio = clock();
	#endif

   n2 = tamano_grilla + 2;

   //inicializa los arrays de dimension n2*n2
   T = (double *)calloc(n2*n2, sizeof(*T));
   Tnew = (double *) calloc(n2*n2, sizeof(*T));

   if (T == NULL || Tnew == NULL) {
      fprintf(stderr, "ERROR!! no hay memoria suficiente\n");
      exit(EXIT_FAILURE);
   }

   // establecer las condiciones limite

   for (i = 1; i <= tamano_grilla; i++) {

      T[(tamano_grilla+1)*n2+i] = Tnew[(tamano_grilla+1)*n2+i] = i * cota_superior / (tamano_grilla + 1);

      T[i*n2+tamano_grilla+1] = Tnew[i*n2+tamano_grilla+1] = i * cota_superior / (tamano_grilla + 1);

   }

   while(maximo_valor_double > tolerancia && iteraciones <= iteraciones_maximas) { 

      ++iteraciones;

      maximo_valor_double = 0.0;

	  //Modificación 2
	  #pragma omp parallel for private(j) 
		  for (int i=1; i<=tamano_grilla; ++i) {

			 for (j=1; j<=tamano_grilla; ++j) {   

				Tnew[i*n2+j] = 0.25*( T[(i-1)*n2+j] + T[(i+1)*n2+j] + T[i*n2+(j-1)] + T[i*n2+(j+1)] );

				//convergencia
                                
				maximo_valor_double = MAX(maximo_valor_double, fabs(Tnew[i*n2+j] - T[i*n2+j]));
			 }
		  }
	  
    
    Tmp = T; 
	T = Tnew; 
	Tnew=Tmp;

    if (iteraciones % 100 == 0){
         printf("Iteracion: %8u, variacion = %12.4lE\n", iteraciones, maximo_valor_double);
	 }
   }


   //Modificación 3
	#ifdef _OPENMP
		double tiempo_final = omp_get_wtime() - tiempo_inicio;
	#else
		double tiempo_final = (clock() - tiempo_inicio) / (double) CLOCKS_PER_SEC;
	#endif

   printf("Tiempo Transcurrido (s) = %.2lf\n", tiempo_final);
   printf("Tamano malla: %u\n", tamano_grilla);
   printf("Parado en la iteracion: %u\n", iteraciones);
   printf("Error maximo: %lE\n", maximo_valor_double);

   free(T);
   free(Tnew);

   //para poder visualizar los datos por pantalla

   system("PAUSE");
}