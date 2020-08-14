#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

const int n = 20;//2000; 	//wielkosc mnozonych macierzy
const int PP = 2; 			//pierwiastek z liczby procesow
const int P = 4;				// liczba procesow

float a[n / PP][n / PP], b[n / PP][n / PP], c[n / PP][n / PP]; // macierze
float aa[n / PP][n / PP], bb[n / PP][n / PP];
//float(*psa)[n / PP], (*psb)[n / PP], (*pra)[n / PP], (*prb)[n / PP]; 

double startwtime1,startwtime2, endwtime;

int main(int argc, char **argv)
{

	FILE *plik;	//plik wejsciowy
	FILE *plik_out; // plik wyjsciowy

	int my_rank, ncpus; //nr procesu, liczba procesów z linii komend  (-np)
	int row, col, mod = 0;
	int data_received = -1;
	int tag = 101;
	int koniec;
	
	MPI_Status  statRecv[2];
	MPI_Request reqSend[2], reqRecv[2];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ncpus); 

	if (my_rank == 0)
	{
		printf("obliczenia metoda Cannona dla tablicy %d x %d elementów \n",n,n);
		startwtime1 = MPI_Wtime();//czas w sekundach
	}

	//wczytanie danych przez proces rank=0
	if (my_rank == 0)
	{
  		plik = fopen("macierz.txt","r");
   		if (plik == NULL) 
   		{
      		printf("Blad otwarcia pliku \"liczby.txt\"\n");
      		koniec=1;
      		MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
      		MPI_Finalize();
      		exit (0);
   		}
		else 
		{
 			koniec=0;
  			MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (koniec) 
			{
				MPI_Finalize();
				exit(0);
			}
	}

	
	if (ncpus != P)
	{ 
       	if (my_rank == 0) 
       	{
       		printf("Blednie wywolano obliczenia iloczynu macierzy metoda cannona na %d procesach - uruchom mpiexec -n %d matrixmult\n",ncpus, P);
       	}
		MPI_Finalize();
		exit (0);
	}

	if (my_rank == 0)
	{
		for(int kk=1; kk<PP*PP; kk++)
		{
			for (int i = 0; i < n / PP; i++)
			{
				for (int j = 0; j < n / PP; j++)
   				{
      				fscanf(plik,"%f",&a[i][j]);
      			}
      		}
			MPI_Isend(a, n*n / PP / PP, MPI_FLOAT, kk, tag, MPI_COMM_WORLD, reqSend);
			//test konca komunikacji
		}
		// czyta dla siebie
		for (int i = 0; i < n / PP; i++)
		{
			for (int j = 0; j < n / PP; j++)
   			{
      			fscanf(plik,"%f",&aa[i][j]);
         	}
         }
		//czyta b do wyslania do innych procesow

		for (int kk =1; kk< PP*PP; kk++)//kolejne identyfikatory procesow
		{
			for (int i = 0; i < n / PP; i++)
			{
				for (int j = 0; j < n / PP; j++)
   				{
      				fscanf(plik,"%f",&bb[i][j]);
      			}
      		}
			MPI_Isend(b, n*n / PP / PP, MPI_FLOAT, kk, tag, MPI_COMM_WORLD, reqSend);
			//test konca komunikacji			
		}

		// czyta dla siebie
		for (int i = 0; i < n / PP; i++)
		{
			for (int j = 0; j < n / PP; j++)
   			{
     			fscanf(plik,"%f",&b[i][j]);
        	}
    	}
    	fclose(plik);
	}
	else
	{
		MPI_Irecv(aa, n*n / PP / PP, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, reqRecv);
		MPI_Wait(reqRecv, statRecv);
		MPI_Irecv(bb, n*n / PP / PP, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqRecv[1]);
		MPI_Wait(&reqRecv[1], &statRecv[1]);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//przygotowanie tablicy wynikowej

	for (int i = 0; i < n / PP; i++)
	{
		for (int j = 0; j < n / PP; j++)
   		{
      		c[i][j]=0;
        }
    }

    row = my_rank / PP; 
	col = my_rank % PP;

	int right = row * PP + (col+1)%PP;
	int left = row * PP + (col-1 + PP)%PP;
	int down = ((row+1)%PP)*PP + col;
	int up = ((row+PP-1)%PP) * PP + col;

	if (my_rank == 0)
	{
		startwtime2 = MPI_Wtime();//czas w sekundach
	}

	//obliczenia iloczynu macierzy zgodnie z algorytmem Cannona 

    for(int i = 1; i < PP; i++)
    {
        if(row >= i)
        {
            MPI_Isend(aa, n*n / PP / PP, MPI_FLOAT, left, tag, MPI_COMM_WORLD, reqSend);
            MPI_Irecv(aa, n*n / PP / PP, MPI_FLOAT, right, tag, MPI_COMM_WORLD, reqRecv);
            MPI_Wait(reqRecv, statRecv);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    for(int i = 1; i < PP; i++)
    {
        if(col >= i)
        {
            MPI_Isend(bb, n*n / PP / PP, MPI_FLOAT, up, tag, MPI_COMM_WORLD, &reqSend[1]);
            MPI_Irecv(bb, n*n / PP / PP, MPI_FLOAT, down, tag, MPI_COMM_WORLD, &reqRecv[1]);
            MPI_Wait(&reqRecv[1], &statRecv[1]);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }


    for (int step = 0; step < PP; step++) {
        for (int i = 0; i < n/PP; i++) {
                for (int j = 0; j < n/PP; j++) {
                        for (int x = 0; x < n/PP; x++) {
                                c[i][j] += aa[i][x] * bb[j][x];
                        }
                }
        }

        if(step < PP - 1)
        {
            MPI_Isend(aa, n*n / PP/ PP, MPI_FLOAT, left, tag, MPI_COMM_WORLD, reqSend);
            MPI_Irecv(aa, n*n / PP / PP, MPI_FLOAT, right, tag, MPI_COMM_WORLD, reqRecv);
            MPI_Wait(reqRecv, statRecv);
            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Isend(aa, n*n / PP / PP, MPI_FLOAT, up, tag, MPI_COMM_WORLD, &reqSend[1]);
            MPI_Irecv(bb, n*n / PP / PP, MPI_FLOAT, down, tag, MPI_COMM_WORLD, &reqRecv[1]);
            MPI_Wait(&reqRecv[1], &statRecv[1]);
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    /*for (int step = 0; step < PP; step++) 
    {
        for (int i = 0; i < n/PP; i++) 
        {
                for (int j = 0; j < n/PP; j++) 
                {
                        for (int x = 0; x < n/PP; x++) 
                        {
                                c[i][j] += aa[i][x] * bb[j][x];
                        }
                }
        }
        MPI_Irecv(aa, n*n / PP / PP, MPI_FLOAT, left, tag, MPI_COMM_WORLD, reqRecv);
        MPI_Irecv(bb, n*n / PP / PP, MPI_FLOAT, down, tag, MPI_COMM_WORLD, &reqRecv[1]);
        MPI_Isend(aa, n*n / PP/ PP, MPI_FLOAT, right, tag, MPI_COMM_WORLD, reqSend);
        MPI_Isend(aa, n*n / PP / PP, MPI_FLOAT, up, tag, MPI_COMM_WORLD, &reqSend[1]);
        MPI_Wait(&reqRecv[1], &statRecv[1]);
    }*/


	//

	if (my_rank == 0)
	{
		endwtime = MPI_Wtime();
		printf("Calkowity czas przetwarzania wynosi %f sekund\n",endwtime - startwtime1);
		printf("Calkowity czas obliczen wynosi %f sekund\n", endwtime - startwtime2);

	}

	// test poprawnosci wyniku - wynik do pliku lub inny sposob

	 MPI_Barrier(MPI_COMM_WORLD);

	if (my_rank == 0)
	{
		plik_out = fopen("wynik.txt", "w");
		if (plik_out == NULL)
		{
			printf("Blad otwarcia pliku \"wynik.txt\"\n");
			exit(0);
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				fprintf(plik_out, "%6.1f", c[i][j]);
			}
			fprintf(plik_out, "\n");
		}
		fclose(plik_out);
	}
	
	MPI_Finalize();
	return 0;
}