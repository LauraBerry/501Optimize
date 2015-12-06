/*
*	name: Laura Berry
* 	ID: 10111166
*	Class: 501 L01
*	this code is taken in part from code given in lecture and tutorial 
*/


#include <sys/time.h>
#include<ctime>
#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<iostream>


using namespace std;

char chunkID[5];
int chunkSize;
char format[5];
char subChunk1ID[5];
int subChunk1Size;
short audioFormat;
short numChannels;
int sampleRate;
int byteRate;
short blockAlign;
short bitsPerSample;
int channelSize;
char subChunk2ID[5];
int subChunk2Size;
short* data;
int numSamples;

char chunkIDIR[5];
int chunkSizeIR;
char formatIR[5];
char subChunk1IDIR[5];
int subChunk1SizeIR;
short audioFormatIR;
short numChannelsIR;
int sampleRateIR;
int byteRateIR;
short blockAlignIR;
short bitsPerSampleIR;
int channelSizeIR;
char subChunk2IDIR[5];
int subChunk2SizeIR;
short* dataIR;
int numSamplesIR;

void print()
{
	chunkID[5] = '\0';
	format[5] = '\0';
	subChunk1ID[5] = '\0';
	subChunk2ID[5] = '\0';
	
	printf("\n============= HEADER INFO =============\n");
	printf(" chunkID:%s\n", chunkID);
	printf(" chunkSize:%d\n", chunkSize);
	printf(" format:%s\n", format);
	printf(" subChunk1ID:%s\n", subChunk1ID);
	printf(" subChunk1Size:%d\n", subChunk1Size);
	printf(" audioFormat:%d\n", audioFormat);
	printf(" numChannels:%d\n", numChannels);
	printf(" sampleRate:%d\n", sampleRate);
	printf(" byteRate:%d\n", byteRate);
	printf(" blockAlign:%d\n", blockAlign);
	printf(" bitsPerSample:%d\n", bitsPerSample);
	printf(" subChunk2ID:%s\n", subChunk2ID);
	printf(" subChunk2Size:%d\n", subChunk2Size);
	cout<<" "<<endl;
}

int loadWave(char* filename, char* response)
{
	FILE* in = fopen(filename, "rb");
	FILE* inResp = fopen(response, "rb");

	if (in != NULL)
	{		
		printf("Reading %s...\n",filename);

		fread(chunkID, 1, 4, in);
		fread(&chunkSize, 1, 4, in);
		fread(format, 1, 4, in);

		//sub chunk 1
		fread(subChunk1ID, 1, 4, in);
		fread(&subChunk1Size, 1, 4, in);
		fread(&audioFormat, 1, 2, in);
		fread(&numChannels, 1, 2, in);
		fread(&sampleRate, 1, 4, in);
		fread(&byteRate, 1, 4, in);
		fread(&blockAlign, 1, 2, in);
		fread(&bitsPerSample, 1, 2, in);		
		
		//read extra bytes
		if(subChunk1Size == 18)
		{
			short empty;
			fread(&empty, 1, 2, in);		
		}
		
		//sub chunk 2
		fread(subChunk2ID, 1, 4, in);
		fread(&subChunk2Size, 1, 4, in);

		//read data		
		int bytesPerSample = bitsPerSample/8;
		numSamples = subChunk2Size / bytesPerSample;
		data = (short*) malloc(sizeof(short) * numSamples);
		
		//fread(data, 1, bytesPerSample*numSamples, in);
		
		int i=0;
		short sample=0;
		while(fread(&sample, 1, bytesPerSample, in) == bytesPerSample)
		{		
			data[i++] = sample;
			sample = 0;			
		}
		
		fclose(in);
		printf("Closing %s...\n",filename);
	}
	else
	{
		printf("Can't open file\n");
		return 0;
	}
	if (inResp != NULL)
	{		
		printf("Reading %s...\n",filename);

		fread(chunkIDIR, 1, 4, inResp);
		fread(&chunkSizeIR, 1, 4, inResp);
		fread(formatIR, 1, 4, inResp);

		//sub chunk 1
		fread(subChunk1IDIR, 1, 4, inResp);
		fread(&subChunk1SizeIR, 1, 4, inResp);
		fread(&audioFormatIR, 1, 2, inResp);
		fread(&numChannelsIR, 1, 2, inResp);
		fread(&sampleRateIR, 1, 4, inResp);
		fread(&byteRateIR, 1, 4, inResp);
		fread(&blockAlignIR, 1, 2, inResp);
		fread(&bitsPerSampleIR, 1, 2, inResp);		
		
		//read extra bytes
		if(subChunk1SizeIR == 18)
		{
			short emptyIR;
			fread(&emptyIR, 1, 2, inResp);		
		}
		
		//sub chunk 2
		fread(subChunk2IDIR, 1, 4, inResp);
		fread(&subChunk2SizeIR, 1, 4, inResp);

		//read data		
		int bytesPerSampleIR = bitsPerSampleIR/8;
		numSamplesIR = subChunk2SizeIR / bytesPerSampleIR;
		dataIR = (short*) malloc(sizeof(short) * numSamplesIR);
		
		//fread(data, 1, bytesPerSample*numSamples, in);
		
		int i=0;
		short sampleIR=0;
		while(fread(&sampleIR, 1, bytesPerSampleIR, inResp) == bytesPerSampleIR)
		{		
			dataIR[i++] = sampleIR;
			sampleIR = 0;			
		}
		
		fclose(inResp);
		printf("Closing %s...\n",inResp);
	}
	else
	{
		printf("Can't open response file\n");
		return 0;
	}
	return 1;
}

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;
    mmax = 2;
    while (n > mmax) 
	{
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / mmax);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		if (j > i) 
		{
			double holder=data[j];
			data[j]=data[i];
			data[i]=data[j];
			holder=data[j+1];
			data[j+1]=data[i+1];
			data[i+1]=data[j+1];
		}
		m = nn;
		while (m >= 2 && j > m) 
		{
			j -= m;
			m >>= 1;
		}
		j += m;
		for (i = m; i <= n; i += istep) 
		{
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
		mmax = istep;
    }
}

double * convolve_d(double signal[], int numOfSamples, double responseArr[], int sizeResponseArr)
{

  /*  Do the convolution  */
	 double * convolvedSignal;
	 double * convolvedResponse;
	 double * output;
	 
	 int i_power = (int)pow(2, (int) log2(numOfSamples) + 1);
	 int i_power_2 = 2 * i_power; 
	 convolvedResponse = (double *) malloc(i_power_2 + i_power_2 + i_power_2 + i_power_2); 
	 convolvedSignal = (double *) malloc(i_power_2 + i_power_2 + i_power_2 + i_power_2); 
	 output = (double *) malloc(i_power_2 + i_power_2 + i_power_2 + i_power_2); 
	 
	for (int i=0; i<i_power_2; i+=2)
	{
		convolvedResponse[i] = 0.0;
		convolvedResponse[i+1] = 0.0;
		convolvedSignal[i] = 0.0;
		convolvedSignal[i+1] = 0.0;
	} 

	 for(int j = 0; j < numOfSamples; j++)
	 {
		convolvedSignal[ j << 1 ] = signal[j];
	 }
	 	 
	 for(int k = 0; k < sizeResponseArr; k++)
	 {
		convolvedResponse[ k << 1 ] = responseArr[k];
		
	 }
	
	 four1(convolvedResponse - 1, i_power, 1);
	 four1(convolvedSignal - 1, i_power, 1);
	
	for(int m = 0 ; m < i_power; m ++)
	{
		output[ m << 1] = (convolvedSignal[m] * convolvedResponse[m]) - (convolvedSignal[m+1] * convolvedResponse[m + 1]);
		output[ (m << 1) + 1] = (convolvedSignal[m+1] * convolvedResponse[m]) - (convolvedSignal[m+1] * convolvedResponse[m + 1]);	
	}
	
	// Revert back
	four1(output - 1, i_power, -1);
	return output;
}

int saveWave(char* filename)
{
	FILE* out = fopen(filename, "wb");

	if (out != NULL)
	{		
		printf("\nWriting %s...\n",filename);

		fwrite(chunkID, 1, 4, out);
		fwrite(&chunkSize, 1, 4, out);
		fwrite(format, 1, 4, out);

		//sub chunk 1
		fwrite(subChunk1ID, 1, 4, out);
		fwrite(&subChunk1Size, 1, 4, out);
		fwrite(&audioFormat, 1, 2, out);
		fwrite(&numChannels, 1, 2, out);
		fwrite(&sampleRate, 1, 4, out);
		fwrite(&byteRate, 1, 4, out);
		fwrite(&blockAlign, 1, 2, out);
		fwrite(&bitsPerSample, 1, 2, out);		
		
		//read extra bytes
		if(subChunk1Size == 18)
		{
			short empty = 0;
			fwrite(&empty, 1, 2, out);		
		}
		
		//sub chunk 2
		fwrite(subChunk2ID, 1, 4, out);
		fwrite(&subChunk2Size, 1, 4, out);

		//read data		
		int bytesPerSample = bitsPerSample / 8;
		int sampleCount =  subChunk2Size / bytesPerSample;
		
		//write the data
		double* newData = (double*) malloc(sizeof(double) * (numSamples+numSamplesIR-1));
		double maxSample = -1;
		double MAX_VAL = 32767.f;	//FIXME: find based on bits per sample
		clock_t start=clock();
		four1(newData, (int)numSamplesIR, (int)numSamples);	
		double* convolvedData= convolve_d((double*)data, (int)numSamples, (double*)dataIR, sizeof(dataIR)); 
		/*for(int i=0; i<numSamples; ++i)
		{			
			//convolve
			for(int j=0; j<numSamplesIR; ++j)
				newData[i+j] += ((float)data[i] / MAX_VAL) * dataIR[j];
			
			//Keep track of max value for scaling
			if(i==0)
				maxSample = newData[0];
			else if(newData[i] > maxSample)
				maxSample = newData[i];
		}*/
		clock_t convolveTime=clock()-start;
		cout<<"\ntime to convolve: ";
		cout<<convolveTime;
		cout<<" cycles\n"<<endl;
		//scale and re write the data
		for(int i=0; i<sampleCount + chunkSize; ++i)
		{
			newData[i] = (newData[i] / maxSample) ;
			short sample = (short) (newData[i] * MAX_VAL);
			fwrite(&sample, 1, bytesPerSample, out);
		}
		
		//clean up
		free(newData);
		fclose(out);
		printf("Closing %s...\n",filename);
	}
	else
	{
		printf("Can't open file\n");
		return 0;
	}
	return 1;
}

int main(int argc, char* argv[])
{
	char* filename = argv[1];
	char* response= argv[2];
	char* outputFile=argv[3];
	clock_t start=clock();
	if(loadWave(filename, response))
		print();
	clock_t loadedFileTime=clock()-start;
	cout <<"time to load file: ";
	cout<<loadedFileTime;
	cout<<" cycles"<<endl;
	saveWave("out.wav");
	free(data);
}
