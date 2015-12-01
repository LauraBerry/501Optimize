/*  Include files  */
#include <stdio.h>
#include <iostream>
#include <math.h>
#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr
using namespace std;
/*
typedef struct HEADER
{
	char RIFF[4];
	unsigned long 	chunkSize;
	char	WAVE[4];
	char	fmt[4];
	unsigned long	subchunk1Size;
	unsigned short	AudioFormat;
	unsigned short Channels;
	unsigned long  SamplesPerSec;
	unsigned long bytesperSec;
	unsigned short blockAlign;
	unsigned short bitsPerSample;
	char	subChunk2ID[4];
	unsigned long subChunk2Size;
}wav_hdr;*/


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

int main(void)
{
  float input_signal[100], impulse_response[20], output_signal[120];
  int input_size, impulse_size, output_size;
  Cwave input, Holder;

/*	wav_hdr wavHeader;
	FILE *wavfile; 
	int headerSize = sizeof(wav_hdr);
	int filelength=0;*/
	cout<<"what is the file name of the input sound?"<<endl;
	cin >> input;
	cin.get();
	cout<<input<<endl;
	/*wavfile=fopen(input);
		Holder.Mix(input);
		Holder.Play();
		
	*/



  /*  Create an example input signal  */
  input_signal[0] = 1.0;
  input_signal[1] = 0.5;
  input_signal[2] = 0.25;
  input_signal[3] = 0.125;
  input_size = 4;
    
  /*  Print out the input signal to the screen  */
  print_vector("Original input signal", input_signal, input_size);

  /*  Create an "identity" impulse response.  The output should be
      the same as the input when convolved with this  */
  impulse_response[0] = 1.0;
  impulse_size = 1;

  /*  Set the expected size of the output signal  */
  output_size = input_size + impulse_size - 1;

  /*  Do the convolution, and print the output signal  */
  convolve(input_signal, input_size, impulse_response, impulse_size,
	   output_signal, output_size);
  print_vector("Output signal using identity IR", output_signal, output_size);


  /*  Create an "inverse" impulse response.  The output should be
      inverted when convolved with this  */
  impulse_response[0] = -1.0;
  impulse_size = 1;

  /*  Set the expected size of the output signal  */
  output_size = input_size + impulse_size - 1;

  /*  Do the convolution, and print the output signal  */
  convolve(input_signal, input_size, impulse_response, impulse_size,
	   output_signal, output_size);
  print_vector("Output signal using inverse IR", output_signal, output_size);


  /*  Create a "scaling" impulse response.  The output should be
      1/2 the amplitude when convolved with this  */
  impulse_response[0] = 0.5;
  impulse_size = 1;

  /*  Set the expected size of the output signal  */
  output_size = input_size + impulse_size - 1;

  /*  Do the convolution, and print the output signal  */
  convolve(input_signal, input_size, impulse_response, impulse_size,
	   output_signal, output_size);
  print_vector("Output signal scaled by 1/2", output_signal, output_size);


  /*  Create a "delay" impulse response.  The output should be
      delayed by 2 samples  */
  impulse_response[0] = 0.0;
  impulse_response[1] = 0.0;
  impulse_response[2] = 1.0;
  impulse_size = 3;

  /*  Set the expected size of the output signal  */
  output_size = input_size + impulse_size - 1;

  /*  Do the convolution, and print the output signal  */
  convolve(input_signal, input_size, impulse_response, impulse_size,
	   output_signal, output_size);
  print_vector("Output delayed 2 samples", output_signal, output_size);


  /*  Create a "delay and scaling" impulse response.  The output should be
      delayed by 2 samples and be 1/2 the amplitude  */
  impulse_response[0] = 0.0;
  impulse_response[1] = 0.0;
  impulse_response[2] = 0.5;
  impulse_size = 3;

  /*  Set the expected size of the output signal  */
  output_size = input_size + impulse_size - 1;

  /*  Do the convolution, and print the output signal  */
  convolve(input_signal, input_size, impulse_response, impulse_size,
	   output_signal, output_size);
  print_vector("Delayed 2 samples, 1/2 amplitude", output_signal, output_size);


  /*  Create an "echo effect".  The output will contain the original signal
      plus a copy delayed by 2 samples and 1/2 the amplitude.  The original
      and copy will overlap starting at the 3rd sample  */
  impulse_response[0] = 1.0;
  impulse_response[1] = 0.0;
  impulse_response[2] = 0.5;
  impulse_size = 3;

  /*  Set the expected size of the output signal  */
  output_size = input_size + impulse_size - 1;

  /*  Do the convolution, and print the output signal  */
  convolve(input_signal, input_size, impulse_response, impulse_size,
	   output_signal, output_size);
  print_vector("Overlapping echo", output_signal, output_size);


  /*  Create an "echo effect" that doesn't overlap.  The output will
      contain the original signal  plus a copy delayed by 5 samples
      and 1/2 the amplitude.  */
  impulse_response[0] = 1.0;
  impulse_response[1] = 0.0;
  impulse_response[2] = 0.0;
  impulse_response[3] = 0.0;
  impulse_response[4] = 0.0;
  impulse_response[5] = 0.5;
  impulse_size = 6;

  /*  Set the expected size of the output signal  */
  output_size = input_size + impulse_size - 1;

  /*  Do the convolution, and print the output signal  */
  convolve(input_signal, input_size, impulse_response, impulse_size,
	   output_signal, output_size);
  print_vector("Non-overlapping echo", output_signal, output_size);


  /*  Interchange the input signal and impulse response.  Since
      convolution is commutative, you should get the same output  */
  convolve(impulse_response, impulse_size, input_signal, input_size,
	   output_signal, output_size);
  print_vector("Same as above, but with interchanged h[] and x[]",
	       output_signal, output_size);

  
  /*  End of program  */
  return 0;
}

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
  if (P != (N + M - 1)) {
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
  for (n = 0; n < N; n++) {
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
