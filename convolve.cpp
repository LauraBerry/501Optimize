/*  Include files  */
#include <stdio.h>
#include <iostream>
#include <math.h>
#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr
using namespace std;

char * inFile;
char * InputResponse;
char * outFile;


/*  Function prototypes  */
void convolve(float x[], int N, float h[], int M, float y[], int P);
void print_vector(char *title, float x[], int N);


/*****************************************************************************
*
*    Function:     main
*
*    Description:  Tests the convolve function with various input signals
*
*****************************************************************************/
/*this code is taken from code provided by Prof. Manzara*/
void method(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) 
	{
		if (j > i) 
		{
			SWAP(data[j], data[i]);
			SWAP(data[j+1], data[i+1]);
		}
		m = nn;
		while (m >= 2 && j > m) 
		{
			j -= m;
			m >>= 1;
		}
		j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}


/*****************************************************************************
*
*    Function:     convolve
*
*    Description:  Convolves two signals, producing an output signal.
*                  The convolution is done in the time domain using the
*                  "Input Side Algorithm" (see Smith, p. 112-115).
*
*    Parameters:   x[] is the signal to be convolved
*                  N is the number of samples in the vector x[]
*                  h[] is the impulse response, which is convolved with x[]
*                  M is the number of samples in the vector h[]
*                  y[] is the output signal, the result of the convolution
*                  P is the number of samples in the vector y[].  P must
*                       equal N + M - 1
*
*****************************************************************************/

void convolve(float x[], int N, float h[], int M, float y[], int P)
{
  int n, m;

  /*  Make sure the output buffer is the right size: P = N + M - 1  */
  if (P != (N + M - 1)) 
  {
    printf("Output signal vector is the wrong size\n");
    printf("It is %-d, but should be %-d\n", P, (N + M - 1));
    printf("Aborting convolution\n");
    return;
  }

  /*  Clear the output buffer y[] to all zero values  */  
  for (n = 0; n < P; n++)
    y[n] = 0.0;

  /*  Do the convolution  */
  /*  Outer loop:  process each input value x[n] in turn  */
  for (n = 0; n < N; n++) 
  {
    /*  Inner loop:  process x[n] with each sample of h[]  */
    for (m = 0; m < M; m++)
      y[n+m] += x[n] * h[m];
  }
}



/*****************************************************************************
*
*    Function:     print_vector
*
*    Description:  Prints the vector out to the screen
*
*    Parameters:   title is a string naming the vector
*                  x[] is the vector to be printed out
*                  N is the number of samples in the vector x[]
*
*****************************************************************************/

void print_vector(char *title, float x[], int N)
{
  int i;

  printf("\n%s\n", title);
  printf("Vector size:  %-d\n", N);
  printf("Sample Number \tSample Value\n");
  for (i = 0; i < N; i++)
    printf("%-d\t\t%f\n", i, x[i]);
}

/*
* function reader
* description: reads the input file
*/
void reader(char * in, char * inResponse)
{
	FILE * input=fopen(in, "rb"); //reading the file in binary
}


int main(int argc, char ** argv)
{
	if(argc<4)
	{
		return 0;
	}
	inFile=argv[1];
	InputResponse=argv[2];
	outFile=argv[3];
	
	//method for reading Input
	return 0;
}