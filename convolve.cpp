//important note everywhere it says SizeOfChunk2 is wrong Laura


/*  Include files  */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <iostream>
#include <string.h>
#include <math.h>
#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr
using namespace std;

char * inFile;
char * InputResponse;
char * outFile;
int inFileSize=0;
int sizeOfSignal=0;


typedef struct header_file
{
    char chunk_id[4];
    int chunk_size;
    char format[4];
    char subchunk1_id[4];
    int subchunk1_size;
    short int audio_format;
    short int num_channels;		// number of channels i.e 1-mono 2-stereo
    int sample_rate;			// sample_rate denotes the sampling rate.
    int byte_rate;
    short int block_align;
    short int bits_per_sample;		// number of bits per sample
} header;
typedef struct header_file* header_p;

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
	cout<<"reader"<<endl;
	FILE * input=fopen(in, "rb"); 						//reading the file in binary
	
	header_p head= (header_p)malloc(sizeof(header));
	
	if(in == NULL)
	{
		//do something
	}
	else
	{
		fseek(input, 0, SEEK_END);
		inFileSize=ftell(input);						//length of input file
		rewind(input);
		
		char * fileContent= (char *) malloc (inFileSize);
		fread(fileContent, inFileSize, 1, input);			//read input file
		rewind(input);
		
		fread(head, 1, sizeof(head), input);			//get header info
		int a=0;
		for(int i=0;i<inFileSize-5;i++)
		{
			if(fileContent[i]=='d' && fileContent[i+1]=='a' && fileContent[i+2]=='t' && fileContent[i+3]==fileContent[i+1])
			{
				a=i+4;
			}
		}
		fseek(input, a, SEEK_SET);
		char * SizeOfChunk2;
		fread(SizeOfChunk2, 4,1, input);
		
		fseek(input, a+4, SEEK_SET);
		char * data= (char*)malloc(sizeof(SizeOfChunk2));
		fread(data,sizeof(SizeOfChunk2),1,input);
		fclose(input);
		
		short * smokeSignal = (short*)malloc(sizeof(SizeOfChunk2));
		int lengthSmokeSignal=sizeof(SizeOfChunk2)/2;
		for(int j=0; j<sizeOfSignal;j++)
		{
			short half= (short) ((unsigned char) sizeof(SizeOfChunk2));
			short otherHalf=(short)((unsigned char)sizeof(SizeofChunk2));
		}			
	}
}


int main(int argc, char ** argv)
{
	cout<<"main"<<endl;
	if(argc<4)
	{
		return 0;
	}
	inFile=argv[1];
	outFile=argv[2];
	InputResponse=argv[3];
	
	reader(inFile, InputResponse);
	
	return 0;
}