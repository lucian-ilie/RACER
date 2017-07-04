#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <cstring>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <omp.h>
#include <assert.h>
using namespace std;

//#define VERBOSE

//*********************************************************************************
//   U(w,l,n,p) gives the total number of reads uncorrectable with witLength = w
//   w = witLength, l = readLength, n = numberOfReads, p = error
//*********************************************************************************/

double U(int64_t witLength, int64_t readLength, int64_t numberOfReads) 
{
	int64_t i,j,k;
	int64_t** fw = new int64_t*[readLength+1];
	double sum=0, temp=0;
	double error = 0.01;
	for (i=0; i<readLength+1; ++i) {
		fw[i] = new int64_t[readLength+1];

		for (j=0; j<=readLength; ++j)
		fw[i][j] = 0;
	}

	for (j=1; (j<=witLength) && (j<=readLength); ++j)
		fw[1][j] = j;

	for (j=witLength+1; (j<=2*witLength-1) && (j<=readLength); ++j)
		fw[1][j] = 2*witLength-j;

	for (j=2*witLength; j<=readLength; ++j) 
		fw[1][j] = 0;

	for (i=2; i<=readLength; ++i) {
		fw[i][1] = 0;

		for (j=2; (j<=witLength+1) && (j<=readLength); ++j) {
			sum += fw[i-1][j-1];
			fw[i][j] = static_cast<int64_t>(sum);
		}

		for (j=witLength+2; j<=readLength; ++j) {
			sum -= fw[i-1][j-witLength-1];
			sum += fw[i-1][j-1];
			fw[i][j] = static_cast<int64_t>(sum);
		}
	}

	sum = 0; temp = 2; k = 1;

	while (k<=readLength) {

		temp = (double)fw[k][readLength] * pow(error,k) * pow(1-error,readLength-k) * (double)numberOfReads;
		sum += temp;
		k++;

		if ((temp>0) && (temp < 1))

		k = readLength+1;

	}
	for (i=0; i<readLength+1; ++i)
		delete [] fw[i];
	delete [] fw;
	
	return(sum);
}

//*********************************************************************************
//   D(w,l,n,L,p) = expected number of destructible reads, that is, correct reads
//   that are turned wrong because they contain va with v containing errors,
//   a correct but v appearing elsewhere as correct vb
//*********************************************************************************/

double D(int64_t witLength, int64_t readLength, int64_t numberOfReads, uint64_t genomeLength) 
{
	// ** v has some erros: (1-(1-p)^w)
	// ** a is correct: (1-p)
	// ** v appears elsewhere: (1-(1-1/4^w)^L)
	// ** followed by b<>a: 3/4
	double error = 0.01;
	double qw = (1 - pow(1-error, witLength)) * (1 - error) * (1 - pow(1 - pow(0.25, witLength), genomeLength)) * 0.75; //** prob 1 wit gone wrong

	double expCorrReads = pow(1 - error, readLength) * (double)numberOfReads;

	//  ** 1 - (1 - qw)^(l-w) = prob 1 read gone wrong

	return ((1 - pow(1 - qw, readLength - witLength)) * expCorrReads);
}

//*********************************************************************************
//   compute the paremeters needed for correction
//   witLength_small = wm and witLength_large = wM
//   no larger than 28 so that it fits into 7 bytes with 2 bits per base
//*********************************************************************************/

void computeWitLength(int64_t readLength, int64_t numberOfReads, uint64_t genomeLength, int64_t &witLength_small, int64_t &witLength_large) 
{
	//** compute wM, wm = witness lengths
	//** wM = min (w | D(w) < 0.0001 * expErrReads)
	//** wm = arg min (U(w) + D(w))
	double error = 0.01;
	double expErrReads = (1 - pow((1 - error), readLength)) * (double)numberOfReads;
	double temp = 0, temp1 = 0, temp2 = 0; 
	witLength_large = 1;

	while (D(witLength_large, readLength, numberOfReads, genomeLength) > 0.0001 * expErrReads) {
		witLength_large = witLength_large + 1;
	}

	if (witLength_large > 25) 
		witLength_large = 25;

	witLength_small = witLength_large;

	temp = U(witLength_small, readLength, numberOfReads) + D(witLength_small, readLength, numberOfReads, genomeLength);

	temp1 = U(witLength_small - 1, readLength, numberOfReads) + D(witLength_small - 1, readLength, numberOfReads, genomeLength);

	temp2 = U(witLength_small + 1, readLength, numberOfReads) + D(witLength_small + 1, readLength, numberOfReads, genomeLength);

	if (temp1 < temp) 
		while (temp1 < temp) {
			witLength_small--;
			temp = U(witLength_small, readLength, numberOfReads) + D(witLength_small, readLength, numberOfReads, genomeLength);
			temp1 = U(witLength_small - 1, readLength, numberOfReads) + D(witLength_small - 1, readLength, numberOfReads, genomeLength);
		}
	else 
		while (temp2 < temp) {
			witLength_small++;
			temp = U(witLength_small, readLength, numberOfReads) + D(witLength_small, readLength, numberOfReads, genomeLength);
			temp2 = U(witLength_small + 1, readLength, numberOfReads) + D(witLength_small + 1, readLength, numberOfReads, genomeLength);
		}

	if (witLength_small > 24) 
		witLength_small = 24;

	if (witLength_small >= witLength_large) 
		witLength_large = witLength_small + 1;
}

//*********************************************************************************
//   compute Threshold = T
//   counts T or larger are correct, smaller are wrong
//*********************************************************************************/

void computeThreshold(int64_t readLength, int64_t numberOfReads, uint64_t genomeLength, int64_t witLength, int64_t &threshold) 
{
	// compute T (done for each iteration)
	// T = min (k | Wc(k) > We(k) + 1)
	double error = 0.01;
	double qc = (readLength - witLength) * pow((1 - error), (witLength + 1)) / (double) genomeLength;
	double qe = (readLength - witLength) * pow((1 - error), witLength) * error / (double) (3 * genomeLength);
	// while (dbinom(T,numberOfReads,qc) <= dbinom(T,numberOfReads,qe)) {
	threshold = 1;

	while (pow(1-error, threshold) * pow(1-qc, numberOfReads-threshold) <= pow(error/3, threshold) * pow(1-qe, numberOfReads-threshold))
		threshold++;
	
	threshold += 2;
}  

//*********************************************************************************
//	Read input from the file provided by the user.  The reads are stored in a
//	character array and then converted to binary and stored in the binReads array.
//*********************************************************************************/
void readInput(uint8_t* &binReads, char* argv[], char &fileFormatFlag, uint64_t &totalBases, int64_t &numberOfReads, int64_t &maxReadLength, int* &goodRead, uint64_t* &readLocation, uint64_t &currentMemory, uint64_t &peakMemory)
{
	time_t t_start, t_end;
	time(&t_start);

	/* Create and initialize variables. */
	uint64_t baseCounter=0, readsSize=0, chunkSize=0, binReadSize=0, offset=0, k=0, i=0, j=0, totalNs=0, badReads=0, binReadsLocation=0, previousPeak=0;
	int64_t currReadLength=0, arraySize=5000000, oldArraySize=0;
	int chooseLetter=0;
	
	/* Open the read file provided by the user. */
	FILE *inReadfile = fopen(argv[1], "r");
	if (!inReadfile) 
	{
		cout << "ERROR: Unable to open text file. " << endl << flush;
		exit( EXIT_FAILURE );
	}
	
	/* Find the size of the input and store it in the readsSize variable. */
	fseek(inReadfile, 0, SEEK_END);
	readsSize = ftell(inReadfile);
	chunkSize = static_cast<uint64_t>(0.1*readsSize);
	rewind(inReadfile);
		
	/* Store the reads file into a character array. */
	unsigned char * charReads;
	charReads = new unsigned char[chunkSize];
	fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile); 


	currentMemory += sizeof(char)*chunkSize;
	if (currentMemory > peakMemory)
	{	previousPeak = peakMemory;
		peakMemory = currentMemory;
#ifdef VERBOSE
		cout << "New peak Memory = " << (peakMemory/1048576) << " MB. charReads array added " << (peakMemory-previousPeak)/1048576 << " MB." << endl << flush;
#endif
	}

    
	/* Determine the file format and read length. */
	while ((i < readsSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
		i++; //blank lines
	if (charReads[i] == '>')
		fileFormatFlag = '>';
	else if (charReads[i] == '@')
		fileFormatFlag = '@';	//fileFormatFlag is '>' for fasta and '@' for fastq files.
	else
	{
		cerr << "Wrong input file format (must be fasta or fastq)" << endl << flush;
		exit(1);
	}
	
	/* Create and initialize the binReads array to store the reads in binary. */
	if (fileFormatFlag == '@')
		binReadSize = static_cast<uint64_t>(floor(readsSize/8));
	else
		binReadSize = static_cast<uint64_t>(floor(readsSize/4));
	
	binReads = new uint8_t[binReadSize];
	for (k=0; k<binReadSize; ++k) 
		binReads[k]=0;
    
	currentMemory += sizeof(uint8_t)*binReadSize;
	if (currentMemory > peakMemory)
	{	previousPeak = peakMemory;
		peakMemory = currentMemory;
#ifdef VERBOSE
		cout << "New peak Memory = " << (peakMemory/1048576) << " MB. binReads array added " << (peakMemory-previousPeak)/1048576 << " MB." << endl << flush;
#endif
    }

    cout << "reading input file..." << flush;

	/* Processing the reads by encoding the base pairs into 2 bits. */
	for (j=0, i=0; j < readsSize; ++i, ++j)
	{
		totalNs = 0;
		currReadLength = 0;
		baseCounter = 0;
		
		if (i == chunkSize)
		{	fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
			i = 0;
		}
		
		while ((i < chunkSize) && (charReads[i] != '\r') && (charReads[i] != '\n'))
		{	i++; //read header
			j++;
			if (i == chunkSize)
			{	fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
				i = 0;
			}
		}
		while ((i < chunkSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
		{	i++; //blank lines
			j++;
			if (i == chunkSize)
			{	fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
				i = 0;
			}
		}
		while ((i < chunkSize) && (charReads[i] != '\r') && (charReads[i] != '\n')) 
		{
			switch(charReads[i]) 
			{
				case 'A':
					binReads[binReadsLocation] <<= 2;
					currReadLength++;
					break;
				case 'C': 
					binReads[binReadsLocation] <<= 2;
					binReads[binReadsLocation] |= 1;
					currReadLength++;
					break;
				case 'G':
					binReads[binReadsLocation] <<= 2;
					binReads[binReadsLocation] |= 2;
					currReadLength++;
					break;
				case 'T':
					binReads[binReadsLocation] <<= 2;
					binReads[binReadsLocation] |= 3;
					currReadLength++;
					break;
				default:
					chooseLetter = rand() % 4;
					if (chooseLetter == 0)
						binReads[binReadsLocation] <<= 2;
					else if (chooseLetter == 1)
					{	binReads[binReadsLocation] <<= 2;
						binReads[binReadsLocation] |= 1;
					}
					else if (chooseLetter == 2)
					{	binReads[binReadsLocation] <<= 2;
						binReads[binReadsLocation] |= 2;
					}
					else
					{	binReads[binReadsLocation] <<= 2;
						binReads[binReadsLocation] |= 3;
					}
					totalNs++;
					currReadLength++;
			}
			totalBases++;
			baseCounter++;
			if ((baseCounter%4)==0)			
				binReadsLocation++;
			i++;
			j++;
			if (i == chunkSize)
			{	fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
				i = 0;
			}
		}
		
		/* Increment the count once 4 bases are encoded which constitutes one byte. */
		if ((baseCounter%4)!=0) 
		{	
			offset = (currReadLength*2)%8;
			binReads[binReadsLocation] <<= (8-offset);
			binReadsLocation++;
			binReads[binReadsLocation]=0;
		}
				
		if(maxReadLength < currReadLength)
			maxReadLength = currReadLength;
		
		/* Check if the read has more than 25% N's. */
		if (totalNs > 0.5*currReadLength)
		{
			goodRead[numberOfReads] = currReadLength*(-1);
			readLocation[numberOfReads] = binReadsLocation-1;
			badReads++;
		}
		else
		{
			goodRead[numberOfReads] = currReadLength;
			readLocation[numberOfReads] = binReadsLocation-1;
		}

		/* Resize the arrays if they are full. */
		if (arraySize == numberOfReads+1)
		{	oldArraySize = arraySize;
			arraySize = 2*arraySize;
			
			uint64_t* readLocationOld;
			int* goodReadOld;
			
			readLocationOld = readLocation;
			goodReadOld = goodRead;
						
			assert(readLocation = new uint64_t[arraySize]);
            

			currentMemory += sizeof(uint64_t)*arraySize;
			if (currentMemory > peakMemory)
			{	previousPeak = peakMemory;
				peakMemory = currentMemory;
#ifdef VERBOSE
				cout << "New peak Memory = " << (peakMemory/1048576) << " MB. readLocation array added " << (peakMemory-previousPeak)/1048576 << " MB." << endl << flush;
#endif
			}
            
			assert(goodRead = new int[arraySize]);
			
			currentMemory += sizeof(int)*arraySize;
			if (currentMemory > peakMemory)
			{	previousPeak = peakMemory;
				peakMemory = currentMemory;
#ifdef VERBOSE
				cout << "New peak Memory = " << (peakMemory/1048576) << " MB. goodRead array added " << (peakMemory-previousPeak)/1048576 << " MB." << endl << flush;
#endif
			}
            
			for(int64_t i=0; i<oldArraySize; ++i)
			{	readLocation[i] = readLocationOld[i];
				goodRead[i] = goodReadOld[i];
			}
			
			delete [] readLocationOld;
			delete [] goodReadOld;
			
			currentMemory -= sizeof(int)*oldArraySize;
			currentMemory -= sizeof(uint64_t)*oldArraySize;
		}
		
		/* Skip over blank lines. */
		while ((i < chunkSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
		{	i++; //blank lines
			j++;
			if (i == chunkSize)
			{	fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
				i = 0;
			}
		}
		
		/* The if block is to process the quality scores in the fastq files. */
		if (fileFormatFlag == '@')
		{
			while ((i < chunkSize) && (charReads[i] != '\r') && (charReads[i] != '\n'))
			{	i++; // next line
				j++;
				if (i == chunkSize)
				{	fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
					i = 0;
				}
			}	
			while ((i < chunkSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
			{	i++; // blank lines
				j++;
				if (i == chunkSize)
				{	fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
					i = 0;
				}
			}
			while ((i < chunkSize) && (charReads[i] != '\r') && (charReads[i] != '\n'))
			{	i++; // quality score
				j++;
				if (i == chunkSize)
				{	fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
					i = 0;
				}
			}
			while ((i < chunkSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
			{	i++; // blank lines
				j++;
				if (i == chunkSize)
				{	fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
					i = 0;
				}
			}
		}
		
		numberOfReads++;
	}
	
	/* Delete the array used to store the input. */ 
	delete [] charReads;
	currentMemory -= sizeof(char)*chunkSize;
	
	fclose(inReadfile);

	time(&t_end);
    
    cout << "done reading input file" << endl << flush;
    
#ifdef VERBOSE
	cout << "Reads file processed. Successfully converted into binary file in " << difftime (t_end,t_start) << " seconds." << endl << flush;
	cout << "Total reads with > 50% N's = " << badReads << endl << flush;
#endif
    
}

//*********************************************************************************
//	Caluclating all the possible values for the reverse complement of an 8 bit 
//	value.  The results are stored in the allRevCompOf8Bit array.
//*********************************************************************************/
void revComplOf8Bit(uint8_t* allRevCompOf8Bit)
{
	uint8_t tempRevComp=0, complementOf_i=0, m=0;
	for (m=0; m<255; ++m) 
	{
		allRevCompOf8Bit[m] = 0;
		complementOf_i = ~ m;
		tempRevComp = 3 & complementOf_i;
		allRevCompOf8Bit[m] |= (tempRevComp << 6);
		complementOf_i = complementOf_i >> 2;
		tempRevComp = 3 & complementOf_i;
		allRevCompOf8Bit[m] |= (tempRevComp << 4);
		complementOf_i = complementOf_i >> 2;
		tempRevComp = 3 & complementOf_i;
		allRevCompOf8Bit[m] |= (tempRevComp << 2);
		complementOf_i = complementOf_i >> 2;
		tempRevComp = 3 & complementOf_i;
		allRevCompOf8Bit[m] |= tempRevComp;
	}
	allRevCompOf8Bit[255]=0;
}

//*********************************************************************************
//  Double the size of the hash table if the number of witnesses is half the size
//	of the current hash table.
//*********************************************************************************/
void reHash(uint64_t* &witness, uint8_t* &counters, uint64_t* hashTableSize, int &hashTableSizeIndex, uint64_t &currHashTabSize)
{
	time_t t_start, t_end;
	time(&t_start);
	
	int y=0, d=0;
	uint64_t oldHashTableSize=0, witnessHashValue=0, tempCurrHashTabSize=0, j=0, x=0, n=0;
	
	uint64_t * oldWitness;
	uint8_t * oldCounters;
	
	/* Find first index with a hash table size larger than twice the currect hash table size. */
	oldHashTableSize = hashTableSize[hashTableSizeIndex];
	cout << endl << "The hash table has filled half its size. Starting to increase its size." << endl << flush;
	
	cout << "Current hash table size is " << oldHashTableSize << flush;

	for (n=hashTableSizeIndex; n<450; n++) 
		if (hashTableSize[n] > 2*oldHashTableSize) 
		{
			hashTableSizeIndex = n;
			break;
		}
	
	/* Transfering the current hash table values to the old hash table array */
	oldWitness = witness;
	oldCounters = counters;
	
	/* Create and initialize the new witness array and counter array */
	tempCurrHashTabSize = hashTableSize[hashTableSizeIndex];
	witness = new uint64_t[tempCurrHashTabSize];
	counters = new uint8_t[8*tempCurrHashTabSize];
	
	cout << ", increasing to " << tempCurrHashTabSize << endl;
	
	for (x=0; x < tempCurrHashTabSize; x++) 
	{
		witness[x] = 0;
		for (y=0; y<8; y++) 
			counters[8*x+y]=0;
	}
	
	/* Creating the hash values of the witness based on the new hash table size. */
	for (j=0; j<oldHashTableSize; j++) 
	{	
		if (oldWitness[j] != 0)
		{		
			witnessHashValue = oldWitness[j] % tempCurrHashTabSize;
			
			/* Use linear probing if the hash value is already taken in the hash table. */
			while (witness[witnessHashValue] != 0)
				if (witnessHashValue < (tempCurrHashTabSize-1)) 
					witnessHashValue++;
				else 
					witnessHashValue = 0;
			
			for (d=0; d<8; d++) 
				counters[8*witnessHashValue+d]=oldCounters[8*j+d];
			
			witness[witnessHashValue] = oldWitness[j];
		}
	}
	
	/* Delete the old witness and counter arrays from memory. */
	delete [] oldWitness;
	delete [] oldCounters;
	
	time(&t_end);
	cout << "Successfully increased the size of the hash table in " << difftime (t_end,t_start) << " seconds." << endl;
	cout << "New hash table size is " << tempCurrHashTabSize << endl << endl;
	
	currHashTabSize = tempCurrHashTabSize;
}

//*********************************************************************************
//  Find all the witnesses and the base before and after the witness.  Stores the 
//	witnesses in the witness array, and the number of occurrences of each base 
//  before and after the witness in the counter array.
//*********************************************************************************/
void buildWitnessesAndCounters(uint8_t* &binReads, uint64_t* &witness, uint8_t* &counters, int64_t startLocation, int64_t endLocation, uint8_t* allRevCompOf8Bit, uint64_t witnessLengthInBits, uint64_t* hashTableSize, int &hashTableSizeIndex, uint64_t &currHashTabSize, int maxReadLength, int* &goodRead, uint64_t* &readLocation, uint64_t &currentMemory, uint64_t &peakMemory)
{   
	/* Initialize the witness and counter arrays. */
	int y=0, reHashing=0;
	for (uint64_t x=0; x < currHashTabSize; ++x) 
	{
		witness[x] = 0;
		for (y=0; y<8; ++y) 
			counters[8*x+y]=0;
	}
	
	uint64_t numOfWitness=0, previousPeak=0;
	
	#pragma omp parallel
	{
		/* Create variables needed in the function. */
		int d=0, f=0, witnessLength = witnessLengthInBits/2;
		uint64_t currReadShare=0, currReadShare_RevComp=0, tempReadShare=0, tempReadShare_RevComp=0, numOfBytesRead=0, numOfReadsProcessed=0, after=0, before=0, witnessValue=0, witnessMinValue=0, witnessHashValue=0, witnessValue_RevCom=0, temp=0, extraRead=0, extraRead_RevComp=0, witnessMask=0, counterForRead=0, counterForRead_RevComp=0, readLengthInBits=0, readLengthInBytes=0, offset=0, k=0, numCollisionsThread=0, witnessHashValueTemp=0, increment=0;

		witnessMask = static_cast<uint64_t>(pow(2, witnessLengthInBits) - 1);
		
		/* Create and initialize arrays to hold a single read and its reverse complement. */
		offset = (maxReadLength*2)%8;
		readLengthInBytes = (offset == 0) ? (maxReadLength*2)/8 : (((maxReadLength*2)/8)+1);
		uint8_t* currRead;                                        
		uint8_t* currRead_RevComp; 
		currRead = new uint8_t[readLengthInBytes];
		currRead_RevComp = new uint8_t[readLengthInBytes];
		
		/* Caluclating the witness pool and the counters. */
		#pragma omp for
		for (int64_t counter=startLocation; counter>=endLocation; --counter)	
		{	
			if (goodRead[counter] > witnessLength)
			{
				/* Initialize variables for every read. */
				readLengthInBits = goodRead[counter]*2;
				offset = (readLengthInBits % 8);
				readLengthInBytes = (offset == 0) ? readLengthInBits/8 : ((readLengthInBits/8)+1);
				currReadShare = 0;
				currReadShare_RevComp = 0;
				tempReadShare = 0;
				tempReadShare_RevComp = 0;
				numOfBytesRead = 0;
				counterForRead = readLengthInBytes - 1;
				counterForRead_RevComp = 0;
				numOfReadsProcessed = 0;
						
				/* Put a single read and its reverse complement into arrays */
				for (d=readLengthInBytes-1,f=0; d >=0 ; --d, ++f) 
				{	
					currRead[d] = binReads[readLocation[counter]-f];
					currRead_RevComp[f] = allRevCompOf8Bit[currRead[d]];
				}
				
				/* Loading the last 64 bits of the read into the currReadShare variable */
				for (k=0; k<=readLengthInBytes && k<8;) 
				{
					tempReadShare = currRead[counterForRead];
					tempReadShare_RevComp = currRead_RevComp[counterForRead_RevComp];
					
					tempReadShare <<= 8*k;
					tempReadShare_RevComp <<= (56 - 8*k);
					
					currReadShare |= tempReadShare;
					currReadShare_RevComp |= tempReadShare_RevComp;
					
					counterForRead--;
					counterForRead_RevComp++;
					
					numOfBytesRead++;
					k++;
				}
				
				/* If readlength/8 is not a whole number then to remove the extra zeroes at the end. */
				if (offset!=0) 
				{
					currReadShare >>= (8-offset);	
					currReadShare_RevComp <<= (8-offset);
					numOfReadsProcessed = 8 - offset;
				}
				
				/* Find each witness and the bases before and after each witness. */
				for (k=0; k<=readLengthInBits-witnessLengthInBits; k+=2,numOfReadsProcessed+=2) 
				{
					/* Find the base after the current witness. */
					if(k!=0)
					{
						after = currReadShare & 3;
						currReadShare >>= 2;
						currReadShare_RevComp <<= 2;
					}
					
					/* Find the current witness and its reverse complement.  Only use the smaller witness value. */
					witnessValue = currReadShare & witnessMask;
					witnessValue_RevCom = (currReadShare_RevComp & (witnessMask << (64-witnessLengthInBits)))>>(64-witnessLengthInBits);
					witnessMinValue = (witnessValue < witnessValue_RevCom) ? witnessValue : witnessValue_RevCom;	
					witnessMinValue = (witnessMinValue == 0) ? 1 : witnessMinValue;
					
					/* Find the base before the current witness. */
					if(k != readLengthInBits-witnessLengthInBits)
						before = (currReadShare >> witnessLengthInBits) & 3;
					
					/* Find the hash function value for the current witness. */
					witnessHashValue = witnessMinValue % currHashTabSize;
					
					/* Use new probing if the hash value is already taken in the hash table. */
					numCollisionsThread=0;
					increment = 1 + (witnessMinValue % hashTableSize[hashTableSizeIndex-1]);
					witnessHashValueTemp = witnessHashValue;
		
					while (witness[witnessHashValueTemp] != witnessMinValue && witness[witnessHashValueTemp] != 0)
					{	
						numCollisionsThread++;
						witnessHashValueTemp = witnessHashValue+(numCollisionsThread * increment);
						while(witnessHashValueTemp>=currHashTabSize)
							witnessHashValueTemp-=currHashTabSize;
					}
					witnessHashValue = witnessHashValueTemp;
					
					if (witness[witnessHashValue] == 0) 
					{
						#pragma omp atomic
						numOfWitness++;
					}
					
					/* Incrementing the appropriate counts for the bases before and after the current witness. */
					witness[witnessHashValue] = witnessMinValue;
					if (witnessMinValue == witnessValue || (witnessMinValue == 1 && witnessValue == 0)) 
					{	
						if(k!=0)
							if (counters[8*witnessHashValue+after] < 255) 
								counters[8*witnessHashValue+after]++;
						if(k != readLengthInBits-witnessLengthInBits)
							if (counters[8*witnessHashValue+before+4] < 255)
								counters[8*witnessHashValue+before+4]++;
					}
					else 
					{
						temp = after;
						after = (~ before) & 3;
						before = (~ temp) & 3;
						if(k!=0)
							if (counters[8*witnessHashValue+before+4] < 255)
								counters[8*witnessHashValue+before+4]++;
						if(k != readLengthInBits-witnessLengthInBits)
							if (counters[8*witnessHashValue+after] < 255) 
								counters[8*witnessHashValue+after]++;
					}
					
					/* Add the next byte from left to the 64 bit variable after processing the right 8 bits. */
					if (numOfReadsProcessed%8==0 && numOfReadsProcessed!=0) 
					{
						if (numOfBytesRead <= readLengthInBytes) 
						{
							extraRead = currRead[counterForRead];
							extraRead_RevComp = currRead_RevComp[counterForRead_RevComp];
							
							counterForRead_RevComp++;
							counterForRead--;
							
							extraRead <<= 56;
							
							currReadShare_RevComp |= extraRead_RevComp;
							currReadShare |= extraRead;
							
							numOfBytesRead++;
						}
					}
				}
			}
			
			/* If the hash table is half full then double the size of the hash table. */
			if (numOfWitness > currHashTabSize*0.9)
			{		
				reHashing = 1;
				counter=endLocation-1;
			}
		}
		
		/* Delete the arrays needed for the function/ */
		delete [] currRead;
		delete [] currRead_RevComp;
	}
	
	if (reHashing == 1)
	{	
		/* Delete the old witness and counter arrays from memory. */
		delete [] witness;
		delete [] counters;
		
		currentMemory -= sizeof(uint8_t)*8*currHashTabSize;
		currentMemory -= sizeof(uint64_t)*currHashTabSize;
        
#ifdef VERBOSE
		cout << "Doubling hash table size." << endl << flush;
#endif
        
		for (int n=hashTableSizeIndex; n<450; n++)
			if (hashTableSize[n] > 2*currHashTabSize) 
			{
				hashTableSizeIndex = n;
				break;
			}
				
		/* Create the new witness and counter arrays from memory. */
		currHashTabSize = hashTableSize[hashTableSizeIndex];
		assert(witness = new uint64_t[currHashTabSize]);
		assert(counters = new uint8_t[8*currHashTabSize]);

		currentMemory += sizeof(uint8_t)*8*currHashTabSize;
		if (currentMemory > peakMemory)
		{	previousPeak = peakMemory;
			peakMemory = currentMemory;
#ifdef VERBOSE
			cout << "New peak Memory = " << (peakMemory/1048576) << " MB. counters array added " << (peakMemory-previousPeak)/1048576 << " MB." << endl << flush;
#endif
		}
		
		currentMemory += sizeof(uint64_t)*currHashTabSize;
		if (currentMemory > peakMemory)
		{	previousPeak = peakMemory;
			peakMemory = currentMemory;
#ifdef VERBOSE
			cout << "New peak Memory = " << (peakMemory/1048576) << " MB. witness array added " << (peakMemory-previousPeak)/1048576 << " MB." << endl << flush;
#endif
		}
        
		buildWitnessesAndCounters(binReads, witness, counters, startLocation, endLocation, allRevCompOf8Bit, witnessLengthInBits, hashTableSize, hashTableSizeIndex, currHashTabSize, maxReadLength, goodRead, readLocation, currentMemory, peakMemory);
	}
	else
    {
#ifdef VERBOSE
		cout << "Total number of witnesses = " << numOfWitness << endl << flush;
#endif
    }
    
}

//*********************************************************************************
//  Correct the errors using the witness and counter arrays.  The base before and
//  after each witness is checked to see if it appears more times than the threshold.
//  If it is below the threshold it is changed to the base that is above the
//  threshold.
//*********************************************************************************/
void correctErrors(uint8_t* binReads, uint64_t* &witness, uint8_t* &counters, int64_t startLocation, int64_t endLocation, uint64_t witnessLengthInBits, uint8_t* allRevCompOf8Bit, int64_t threshold, uint64_t &numOfPosChanged, uint64_t &currHashTabSize, int maxReadLength, int* &goodRead, uint64_t* &readLocation, uint64_t* hashTableSize, int &hashTableSizeIndex)
{	
	uint64_t numOfAmbiguousCorrectionsTotal=0, numOfPosChangedAfterTotal=0, numOfPosChangedBeforeTotal=0;
	numOfPosChanged=0;
	
	#pragma omp parallel
	{
		/* Create variables needed in the function. */
		int f=0, d=0, j=0, corrected=0, witnessLength = witnessLengthInBits/2;
		
		uint64_t currReadShare=0, currReadShare_RevComp=0, tempReadShare=0, tempReadShare_RevComp=0, numOfReadsProcessed=0, after=0, before=0, witnessValue=0, witnessMinValue=0, witnessHashValue=0, witnessValue_RevCom=0, extraRead=0, extraRead_RevComp=0, after_count=0, before_count=0, maskForBaseInsertion=0, maskForBaseInsertion_RevComp=0, numThreshold=0, numOfAmbiguousCorrections=0, numOfPosChangedAfter=0, numOfPosChangedBefore=0, numOfBytesRead=0, k=0, readLengthInBytes=0, readLengthInBits=0, offset=0, counterForRead=0, counterForRead_RevComp=0, numCollisionsThread=0, witnessHashValueTemp=0, increment=0;
		
		/* Set masks used for correcting the witnesses. */
		maskForBaseInsertion = 3;
		maskForBaseInsertion_RevComp = 3;
		maskForBaseInsertion = ~ (maskForBaseInsertion << witnessLengthInBits);
		maskForBaseInsertion_RevComp = ~ (maskForBaseInsertion_RevComp << (62 - witnessLengthInBits));
		uint64_t witnessMask = static_cast<uint64_t>(pow(2, witnessLengthInBits) - 1);
		
		/* Create and initialize arrays to hold a single read and its reverse complement. */
		offset = (maxReadLength*2)%8;
		readLengthInBytes = (offset == 0) ? (maxReadLength*2)/8 : (((maxReadLength*2)/8)+1);
		uint8_t* currRead;
		uint8_t* currRead_RevComp; 
		currRead = new uint8_t[readLengthInBytes];
		currRead_RevComp = new uint8_t[readLengthInBytes];

		/* Initialize variables at the start of each iteration. */
		numOfPosChangedAfter=0;
		numOfPosChangedBefore=0;
		numOfAmbiguousCorrections=0;
        
		/* Correct the reads for this split. */
		#pragma omp for
		for (int64_t counter=startLocation; counter>=endLocation; --counter)
		{		
			if (goodRead[counter] > witnessLength)
			{
				/* Initialize variables for every read. */
				readLengthInBits = goodRead[counter]*2;
				offset = (readLengthInBits % 8);
				readLengthInBytes = (offset == 0) ? readLengthInBits/8 : ((readLengthInBits/8)+1);
				currReadShare = 0;
				currReadShare_RevComp = 0;
				tempReadShare = 0;
				tempReadShare_RevComp = 0;
				numOfBytesRead = 0;
				counterForRead = readLengthInBytes - 1;
				counterForRead_RevComp = 0;
				numOfReadsProcessed = 0;
				
				/* Put a single read and its reverse complement into arrays. */
				for (d=readLengthInBytes-1,f=0; d >=0 ; --d, ++f) 
				{
					currRead[d] = binReads[readLocation[counter]-f];
					currRead_RevComp[f] = allRevCompOf8Bit[currRead[d]];
				}
				
				/* Loading the last 64 bits of the read into the currReadShare variable. */
				for (k=0; k<=readLengthInBytes && k<8;) 
				{
					tempReadShare = currRead[counterForRead];
					tempReadShare_RevComp = currRead_RevComp[counterForRead_RevComp];
					
					tempReadShare <<= 8*k;
					tempReadShare_RevComp <<= (56 - 8*k);
					
					currReadShare |= tempReadShare;
					currReadShare_RevComp |= tempReadShare_RevComp;
					
					counterForRead--;
					counterForRead_RevComp++;
					
					numOfBytesRead++;
					k++;
				}
				
				/* If readlength/8 is not a whole number then to remove the extra zeroes at the end. */
				if (offset!=0) 
				{
					currReadShare >>= (8-offset);
					currReadShare_RevComp <<= (8-offset);
					numOfReadsProcessed = 8 - offset;
				}
				
				
				/* Find and correct each witness and the bases before and after each witness. */
				for (k=0; k<=readLengthInBits-witnessLengthInBits; k+=2,numOfReadsProcessed+=2) 
				{
					/* Find the base after the current witness. */
					if(k!=0)
					{
						after = currReadShare & 3;
						currReadShare >>= 2;
						currReadShare_RevComp <<= 2;
					}
					
					/* Find the current witness and its reverse complement.  Only use the smaller witness value. */
					witnessValue = currReadShare & witnessMask;
					witnessValue_RevCom = (currReadShare_RevComp & (witnessMask << (64-witnessLengthInBits)))>>(64-witnessLengthInBits);
					witnessMinValue = (witnessValue < witnessValue_RevCom) ? witnessValue : witnessValue_RevCom;
					witnessMinValue = (witnessMinValue == 0) ? 1 : witnessMinValue;
					
					/* Find the base before the current witness. */
					if(k != readLengthInBits-witnessLengthInBits)
						before = (currReadShare >> witnessLengthInBits) & 3;
					
					/* Find the hash function value for the current witness. */
					witnessHashValue = witnessMinValue % currHashTabSize;
					
					/* Use new probing if the hash value is already taken in the hash table. */
					numCollisionsThread=0;
					increment = 1 + (witnessMinValue % hashTableSize[hashTableSizeIndex-1]);
					witnessHashValueTemp = witnessHashValue;
		
					while (witness[witnessHashValueTemp] != witnessMinValue && witness[witnessHashValueTemp] != 0)
					{	
						numCollisionsThread++;
						witnessHashValueTemp = witnessHashValue+(numCollisionsThread * increment);
						while(witnessHashValueTemp>=currHashTabSize)
							witnessHashValueTemp-=currHashTabSize;
					}
					witnessHashValue = witnessHashValueTemp;
					
					/* Find and correct any errors in the current witness. */
					if (witnessMinValue == witnessValue || witnessMinValue == 1)
					{
						if(k != readLengthInBits-witnessLengthInBits)
						{
							/* Set the default base to the existing base. */
							before_count = counters[8*witnessHashValue+before+4];
							corrected = 0; // Flag to see if the witness has been corrected.
							
							/* Check to see if the base before the witness should be corrected. */
							for (j=4; j<8; ++j) 
							{
								if (counters[8*witnessHashValue+j] > before_count && counters[8*witnessHashValue+j] > threshold) 
								{
									before = j-4;
									before_count = counters[8*witnessHashValue+j];
									corrected = 1;
								}
							}
							
							/* If corrected then check if there is more than one base above the threshold. */
							if (corrected == 1) 
							{
								numThreshold = 0;
								for (j=4; j<8; ++j) 
								{
									if (counters[8*witnessHashValue+j] > threshold) 
									{
										numThreshold++;
									}
								}
								if (numThreshold > 1) 
								{
									numOfAmbiguousCorrections++;
								}
							}
							
							/* If corrected then increment the value of the corrected base. */
							if (corrected == 1 && numThreshold == 1) 
							{
								
								currRead[readLengthInBytes - 1 - (static_cast<uint64_t>(floor((numOfReadsProcessed+witnessLengthInBits)/8)))] = ((currRead[readLengthInBytes - 1 - (static_cast<uint64_t>(floor((numOfReadsProcessed+witnessLengthInBits)/8)))] & (~ (3 << (((numOfReadsProcessed+witnessLengthInBits))%8)))) | (before << (((numOfReadsProcessed+witnessLengthInBits))%8)));
								
								/* Adding the corrected base to the currReadShare array. */
								tempReadShare = currReadShare & maskForBaseInsertion;
								currReadShare = tempReadShare | (before << witnessLengthInBits);
								tempReadShare = currReadShare_RevComp & maskForBaseInsertion_RevComp;
								before = (~ before) & 3;
								currReadShare_RevComp = tempReadShare | (before << (62-witnessLengthInBits));
								numOfPosChangedBefore++;
							}
						}
						
						/* Get corrected value of the after and before base value. */
						if (k!=0) 
						{
							/* Set the default base to the existing base. */
							after_count = counters[8*witnessHashValue+after];
							corrected = 0; // Flag to see if the witness has been corrected.
							
							/* Check to see if the base after the witness should be corrected. */
							for (j=0; j<4; ++j) 
							{
								if (counters[8*witnessHashValue+j] > after_count && counters[8*witnessHashValue+j] > threshold) 
								{
									after = j;
									after_count = counters[8*witnessHashValue+j];
									corrected = 1;
								}
							}
							
							/* If corrected then check if there is more than one base above the threshold. */
							if (corrected == 1) 
							{
								numThreshold = 0;
								for (j=0; j<4; ++j) 
								{
									if (counters[8*witnessHashValue+j] > threshold) 
									{
										numThreshold++;
									}
								}
								if (numThreshold > 1) 
								{
									numOfAmbiguousCorrections++;
								}
							}
							
							/* If corrected then increment the value of the corrected base. */
							if (corrected == 1 && numThreshold == 1) 
							{
								currRead[readLengthInBytes - 1 - (static_cast<uint64_t>(floor((numOfReadsProcessed - 2)/8)))] = ((currRead[readLengthInBytes - 1 - (static_cast<uint64_t>(floor((numOfReadsProcessed - 2)/8)))] & (~ (3 << ((numOfReadsProcessed-2)%8)))) | (after << ((numOfReadsProcessed-2)%8)));
								numOfPosChangedAfter++;
							}
							
						}
					}
					else 
					{
						if(k != readLengthInBits-witnessLengthInBits)
						{
							/* Set the default base to the existing base. */
							after_count = counters[8*witnessHashValue+((~ before) & 3)];
							corrected = 0; // Flag to see if the witness has been corrected.
							
							/* Check to see if the base after the witness should be corrected. */
							for (j=0; j<4; ++j) 
							{
								if (counters[8*witnessHashValue+j] > after_count && counters[8*witnessHashValue+j] > threshold) 
								{
									before = j;
									after_count = counters[8*witnessHashValue+j];
									corrected = 1;
								}
							}
							
							/* If corrected then check if there is more than one base above the threshold. */
							if (corrected == 1) 
							{
								numThreshold = 0;
								for (j=0; j<4; ++j) 
								{
									if (counters[8*witnessHashValue+j] > threshold) 
									{
										numThreshold++;
									}
								}
								if (numThreshold > 1) 
								{
									numOfAmbiguousCorrections++;
								}
							}
							
							/* If corrected then increment the value of the corrected base. */
							if (corrected == 1 && numThreshold == 1) 
							{
								currRead[readLengthInBytes - 1 - (static_cast<uint64_t>(floor((numOfReadsProcessed+witnessLengthInBits)/8)))] = ((currRead[readLengthInBytes - 1 - (static_cast<uint64_t>(floor((numOfReadsProcessed+witnessLengthInBits)/8)))] & (~ (3 << (((numOfReadsProcessed+witnessLengthInBits))%8)))) | (((~ before) & 3) << (((numOfReadsProcessed+witnessLengthInBits))%8)));
								numOfPosChangedBefore++;
								
								/* Adding the corrected base to the currReadShare array. */
								tempReadShare = currReadShare & maskForBaseInsertion;
								currReadShare = tempReadShare | (((~ before) & 3) << witnessLengthInBits);
								tempReadShare = currReadShare_RevComp & maskForBaseInsertion_RevComp;
								currReadShare_RevComp = tempReadShare | (before << (62-witnessLengthInBits));
							}
							
						}
						
						/* Get corrected value of the after and before base value. */
						if (k!=0) 
						{
							/* Set the default base to the existing base. */
							before_count = counters[8*witnessHashValue+(((~ after)&3) + 4)];
							corrected = 0; // Flag to see if the witness has been corrected.
							
							/* Check to see if the base after the witness should be corrected. */
							for (j=4; j<8; ++j) 
							{
								if (counters[8*witnessHashValue+j] > before_count && counters[8*witnessHashValue+j] > threshold) 
								{
									after = j-4;
									before_count = counters[8*witnessHashValue+j];
									corrected = 1;
								}
							}
							
							/* If corrected then check if there is more than one base above the threshold. */
							if (corrected == 1) 
							{
								numThreshold = 0;
								for (j=4; j<8; ++j) 
								{
									if (counters[8*witnessHashValue+j] > threshold) 
									{
										numThreshold++;
									}
								}
								if (numThreshold > 1) 
								{
									numOfAmbiguousCorrections++;
								}
							}
							
							/* If corrected then increment the value of the corrected base. */
							if (corrected == 1 && numThreshold == 1) 
							{
								currRead[readLengthInBytes - 1 - (static_cast<uint64_t>(floor((numOfReadsProcessed - 2)/8)))] = ((currRead[readLengthInBytes - 1 - (static_cast<uint64_t>(floor((numOfReadsProcessed - 2)/8)))] & (~ (3 << ((numOfReadsProcessed-2)%8)))) | (((~ after)&3) << ((numOfReadsProcessed-2)%8)));
								numOfPosChangedAfter++;
							}
							
						}
					}
					
					/* Add the next byte from left to the 64 bit variable after processing the right 8 bits. */
					if (numOfReadsProcessed%8==0 && numOfReadsProcessed!=0) 
					{
						if (numOfBytesRead <= readLengthInBytes) 
						{
							extraRead = currRead[counterForRead];
							extraRead_RevComp = currRead_RevComp[counterForRead_RevComp];
							
							counterForRead_RevComp++;
							counterForRead--;
							
							extraRead <<= 56;
							
							currReadShare_RevComp |= extraRead_RevComp;
							currReadShare |= extraRead;
							
							numOfBytesRead++;
						}
					}
					
				}
				
				/* Store the corrected read back into the binReads array. */
				for (d=readLengthInBytes-1, f=0; d >=0 ; --d, ++f)
					binReads[readLocation[counter]-f] = currRead[d];
			}
		}
		
		#pragma omp critical (print)
		{
			numOfPosChangedAfterTotal += numOfPosChangedAfter;
			numOfPosChangedBeforeTotal += numOfPosChangedBefore;
			numOfAmbiguousCorrectionsTotal += numOfAmbiguousCorrections;
		}
		
		/* Delete the arrays needed for the function. */
		delete [] currRead;
		delete [] currRead_RevComp;
	}
	
	numOfPosChanged = numOfPosChangedAfterTotal + numOfPosChangedBeforeTotal;
    
#ifdef VERBOSE	
	cout << "Number of changed positions\t\t" << numOfPosChanged << endl << flush;
	cout << "Number of changed positions after\t" << numOfPosChangedAfterTotal << endl << flush;
	cout << "Number of changed positions before\t" << numOfPosChangedBeforeTotal << endl << flush;
	cout << "Number of ambiguous corrections = " << numOfAmbiguousCorrectionsTotal << endl << flush;
#endif
    
}

//*********************************************************************************
//	Write the corrected reads to the output file.  The function reads in the input
//	file and replaces the original reads with the corrected reads.
//*********************************************************************************/
void writeOutput(uint8_t * &binReads, char* argv[], char &fileFormatFlag, int* goodRead, uint64_t* &readLocation)
{
	/* Create and initialize variables used in the function. */
	uint64_t readsSize=0, readIndex=0, goodReadIndex=0, chunkSize=0, readLengthInBytes=0, readLengthInBits=0, offset=0, i=0, k=0, j=0;
	int q=0, offsetCount=0;
	uint8_t baseValue=0;
	char alphabet[4]={'A','C','G','T'};

	/* Open input file provided by the user. */
	FILE *inReadfile = fopen(argv[1], "r");
	if (!inReadfile) 
	{
		cout << "ERROR: Unable to open text file." << endl << flush;
		exit( EXIT_FAILURE );
	}
	
	/* Find the size of the input and store in the readsSize variable. */
	fseek(inReadfile, 0, SEEK_END);
	readsSize = ftell(inReadfile);
	chunkSize = static_cast<uint64_t>(0.1*readsSize);
	rewind(inReadfile);
		
	/* Create and initialize an array to store the input file. */
	unsigned char * charReads;
	charReads = new unsigned char[chunkSize];
	fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile); 

	/* Open the corrected read file. */
	FILE *correctedReadfile; 
	correctedReadfile = fopen(argv[2],"w"); 
	
	/* Put the corrected reads into the charReads array. */
	for (k=0, i=0; k < readsSize; ++i, ++k) 
	{
		if (i == chunkSize && k < readsSize)
		{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
			fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
			i = 0;
		}
		
		if (goodRead[goodReadIndex] > 0)
		{	
			readLengthInBits = goodRead[goodReadIndex]*2;
			offset = (readLengthInBits % 8);
			readLengthInBytes = (offset == 0) ? readLengthInBits/8 : ((readLengthInBits/8)+1);
			
			while ((i < readsSize) && (charReads[i] != '\r') && (charReads[i] != '\n'))
			{	i++; //read header
				k++;
				if (i == chunkSize)
				{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
					fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
					i = 0;
				}
			}
		
			while ((i < readsSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
			{	i++; // blank lines
				k++;
				if (i == chunkSize)
				{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
					fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
					i = 0;
				}
			}
			while ((i < readsSize) && (charReads[i] != '\r') && (charReads[i] != '\n')) 
			{
				for (j=0; j<readLengthInBytes; ++j) 
				{
					if ((j == readLengthInBytes-1) && (offset != 0) )
						offsetCount = 7-offset;
					else 
						offsetCount = 0;
					
					for (q=6; q>=offsetCount; q-=2) 
					{
						baseValue = binReads[readIndex+j] >> q;
						baseValue &= 3;
						charReads[i]=alphabet[baseValue];
						i++;
						k++;
						if (i == chunkSize)
						{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
							fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
							i = 0;
						}
					}
				}
			}
			readIndex+=readLengthInBytes;
			goodReadIndex++;
			
			while ((i < readsSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
			{	i++; // blank lines
				k++;
				if (i == chunkSize)
				{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
					fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
					i = 0;
				}
			}			
			/* Process the quality scores in the fastq files. */
			if (fileFormatFlag == '@')	
			{
				while ((i < readsSize) && (charReads[i] != '\r') && (charReads[i] != '\n'))
				{
					i++; // next line
					k++;
					if (i == chunkSize)
					{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
						fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
						i = 0;
					}
				}
				while ((i < readsSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
				{
					i++; // blank lines
					k++;
					if (i == chunkSize)
					{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
						fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
						i = 0;
					}
				}
				while ((i < readsSize) && (charReads[i] != '\r') && (charReads[i] != '\n'))
				{
					i++; // quality scores
					k++;
					if (i == chunkSize)
					{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
						fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
						i = 0;
					}
				}
				while ((i < readsSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
				{
					i++; // blank lines
					k++;
					if (i == chunkSize)
					{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
						fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
						i = 0;
					}
				}
			}		
		}
		else
		{
			readLengthInBits = (goodRead[goodReadIndex]*(-1))*2;
			offset = (readLengthInBits % 8);
			readLengthInBytes = (offset == 0) ? readLengthInBits/8 : ((readLengthInBits/8)+1);
			
			while ((i < readsSize) && (charReads[i] != '\r') && (charReads[i] != '\n'))
			{	i++; //read header
				k++;
				if (i == chunkSize)
				{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
					fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
					i = 0;
				}
			}
		
			while ((i < readsSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
			{	i++; // blank lines
				k++;
				if (i == chunkSize)
				{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
					fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
					i = 0;
				}
			}
			while ((i < readsSize) && (charReads[i] != '\r') && (charReads[i] != '\n')) 
			{	i++; // keep uncorrected read
				k++;
				if (i == chunkSize)
				{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
					fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
					i = 0;
				}
			}
				
			readIndex+=readLengthInBytes;
			goodReadIndex++;
			
			while ((i < readsSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
			{	i++; // blank lines
				k++;
				if (i == chunkSize)
				{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
					fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
					i = 0;
				}
			}			
			/* Process the quality scores in the fastq files. */
			if (fileFormatFlag == '@')	
			{
				while ((i < readsSize) && (charReads[i] != '\r') && (charReads[i] != '\n'))
				{
					i++; // next line
					k++;
					if (i == chunkSize)
					{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
						fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
						i = 0;
					}
				}
				while ((i < readsSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
				{
					i++; // blank lines
					k++;
					if (i == chunkSize)
					{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
						fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
						i = 0;
					}
				}
				while ((i < readsSize) && (charReads[i] != '\r') && (charReads[i] != '\n'))
				{
					i++; // quality scores
					k++;
					if (i == chunkSize)
					{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
						fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
						i = 0;
					}
				}
				while ((i < readsSize) && ((charReads[i] == '\r') || (charReads[i] == '\n')))
				{
					i++; // blank lines
					k++;
					if (i == chunkSize)
					{	fwrite (charReads , 1 , chunkSize , correctedReadfile);
						fread((void *) &charReads[0], sizeof(char), chunkSize, inReadfile);
						i = 0;
					}
				}
			}		
		}
	}
	fwrite (charReads , 1 , i-1 , correctedReadfile);
	
	fclose(inReadfile);
	fclose(correctedReadfile);
	
	delete [] charReads;
}

int main(int argc, char* argv[]) 
{
	// Check Arguments
	if (argc != 4) {
		cout << "ERROR: Wrong number of arguments!\n";
		cout << "Usage: ./RACER <inputReads> <correctedReads> <genomeLength>" << endl << flush;
		cout << "Please consult the readme file for more information on how to run the program." << endl << flush;
		exit(EXIT_FAILURE);
	}

	time_t t_start=0, t_end=0, t_witnessStart=0, t_witnessEnd=0, t_correctStart=0, t_correctEnd=0;
	double witnessTime=0, correctingTime=0;
	time(&t_start);

	/* Create and initialize variables. */
	char fileFormatFlag;
	int iteration=0, coverage=0, hashTableSizeIndex=0, n=0;
	
	int64_t maxReadLength = 0, readsInEachSplit=0, numberOfReads=0, witnessLength=0, previousWitnessLength=0, witLength_small=0, witLength_large=0, threshold=0, endLocation=0, startLocation=0;
	
	uint64_t numOfPosChanged=0, witnessLengthInBits=0, totalBasesPerSplit=0, totalBases=0, numOfReadSplits=0, numOfWitAbovThresholdForAfterBase=0, numOfWitAbovThresholdForBeforeBase=0, currHashTabSize=0, genomeLength = strtoull(argv[3],NULL,10), currentMemory=0, peakMemory=0, previousPeak=0, x=0, y=0, z=0;
	
	uint8_t * binReads;	// Array to store read information in binary.
	uint64_t * witness;	// Array to store the witnesses.
	uint8_t * counters;	// Array to store the counters for each witness.
	
	/* Create array to mark reads with too many N's. */
	int* goodRead;
	assert(goodRead = new int[5000000]);
	
	currentMemory += sizeof(int)*5000000;
	if (currentMemory > peakMemory)
	{	previousPeak = peakMemory;
		peakMemory = currentMemory;
#ifdef VERBOSE
		cout << "New peak Memory = " << (peakMemory/1048576) << " MB. goodRead array added " << (peakMemory-previousPeak)/1048576 << " MB." << endl << flush;
#endif
	}
    
	/* Create array to mark the loaction of each read in binReads. */
	uint64_t* readLocation;
	assert(readLocation =  new uint64_t[5000000]);
    
	currentMemory += sizeof(uint64_t)*5000000;
	if (currentMemory > peakMemory)
	{	previousPeak = peakMemory;
		peakMemory = currentMemory;
#ifdef VERBOSE
		cout << "New peak Memory = " << (peakMemory/1048576) << " MB. readLocation array added " << (peakMemory-previousPeak)/1048576 << " MB." << endl << flush;
#endif
	}
    
	/* Create and initialize array to store hash table sizes. All values are prime numbers. */
	uint64_t hashTableSize[450]={10, 20, 30, 40, 50, 100, 200, 500, 1200, 2500, 1769627, 1835027, 1900667, 1966127, 2031839, 2228483, 2359559, 2490707, 2621447, 2752679, 2883767, 3015527, 3145739, 3277283, 3408323, 3539267, 3670259, 3801143, 3932483, 4063559, 4456643, 4718699, 4980827, 5243003, 5505239, 5767187, 6029603, 6291563, 6553979, 6816527, 7079159, 7340639, 7602359, 7864799, 8126747, 8913119, 9437399, 9962207, 10485767, 11010383, 11534819, 12059123, 12583007, 13107923, 13631819, 14156543, 14680067, 15204467, 15729647, 16253423, 17825999, 18874379, 19923227, 20971799, 22020227, 23069447, 24117683, 25166423, 26214743, 27264047, 28312007, 29360147, 30410483, 31457627, 32505983, 35651783, 37749983, 39845987, 41943347, 44040383, 46137887, 48234623, 50331707, 52429067, 54526019, 56623367, 58720307, 60817763, 62915459, 65012279, 71303567, 75497999, 79691867, 83886983, 88080527, 92275307, 96470447, 100663439, 104858387, 109052183, 113246699, 117440699, 121635467, 125829239, 130023683, 142606379, 150994979, 159383759, 167772239, 176160779, 184549559, 192938003, 201327359, 209715719, 218104427, 226493747, 234882239, 243269639, 251659139, 260047367, 285215507, 301989959, 318767927, 335544323, 352321643, 369100463, 385876703, 402654059, 419432243, 436208447, 452986103, 469762067, 486539519, 503316623, 520094747, 570425399, 603979919, 637534763, 671089283, 704643287, 738198347, 771752363, 805307963, 838861103, 872415239, 905971007, 939525143, 973079279, 1006633283, 1040187419, 1140852767, 1207960679, 1275069143, 1342177379, 1409288183, 1476395699, 1543504343, 1610613119, 1677721667, 1744830587, 1811940419, 1879049087, 1946157419, 2013265967, 2080375127, 2281701827, 2415920939, 2550137039, 2684355383, 2818572539, 2952791147, 3087008663, 3221226167, 3355444187, 3489661079, 3623878823, 3758096939, 3892314659, 4026532187, 4160749883, 4563403379, 4831838783, 5100273923, 5368709219, 5637144743, 5905580687, 6174015503, 6442452119, 6710886467, 6979322123, 7247758307, 7516193123, 7784629079, 8053065599, 8321499203, 9126806147, 9663676523, 10200548819, 10737418883, 11274289319, 11811160139, 12348031523, 12884902223, 13421772839, 13958645543, 14495515943, 15032386163, 15569257247, 16106127887, 16642998803, 18253612127, 19327353083, 20401094843, 21474837719, 22548578579, 23622320927, 24696062387, 25769803799, 26843546243, 27917287907, 28991030759, 30064772327, 31138513067, 32212254947, 33285996803, 36507222923, 38654706323, 40802189423, 42949673423, 45097157927, 47244640319, 49392124247, 51539607599, 53687092307, 55834576979, 57982058579, 60129542339, 62277026327, 64424509847, 66571993199, 73014444299, 77309412407, 81604379243, 85899346727, 90194314103, 94489281203, 98784255863, 103079215439, 107374183703, 111669150239, 115964117999, 120259085183, 124554051983, 128849019059, 133143986399, 146028888179, 154618823603, 163208757527, 171798693719, 180388628579, 188978561207, 197568495647, 206158430447, 214748365067, 223338303719, 231928234787, 240518168603, 249108103547, 257698038539, 266287975727, 292057776239, 309237645803, 326417515547, 343597385507, 360777253763, 377957124803, 395136991499, 412316861267, 429496730879, 446676599987, 463856468987, 481036337207, 498216206387, 515396078039, 532575944723, 584115552323, 618475290887, 652835029643, 687194768879, 721554506879, 755914244627, 790273985219, 824633721383, 858993459587, 893353198763, 927712936643, 962072674643, 996432414899, 1030792152539, 1065151889507, 1168231105859, 1236950582039, 1305670059983, 1374389535587, 1443109012607, 1511828491883, 1580547965639, 1649267441747, 1717986918839, 1786706397767, 1855425872459, 1924145348627, 1992864827099, 2061584304323, 2130303780503, 2336462210183, 2473901164367, 2611340118887, 2748779070239, 2886218024939, 3023656976507, 3161095931639, 3298534883999, 3435973836983, 3573412791647, 3710851743923, 3848290698467, 3985729653707, 4123168604483, 4260607557707, 4672924419707, 4947802331663, 5222680234139, 5497558138979, 5772436047947, 6047313952943, 6322191860339, 6597069767699, 6871947674003, 7146825580703, 7421703488567, 7696581395627, 7971459304163, 8246337210659, 8521215117407, 9345848837267, 9895604651243, 10445360463947, 10995116279639, 11544872100683, 12094627906847, 12644383722779, 13194139536659, 13743895350023, 14293651161443, 14843406975659, 15393162789503, 15942918604343, 16492674420863, 17042430234443, 18691697672867, 19791209300867, 20890720927823, 21990232555703, 23089744183799, 24189255814847, 25288767440099, 26388279068903, 27487790694887, 28587302323787, 29686813951463, 30786325577867, 31885837205567, 32985348833687, 34084860462083, 37383395344739, 39582418600883, 41781441856823, 43980465111383, 46179488367203, 48378511622303, 50577534878987, 52776558134423, 54975581392583, 57174604644503, 59373627900407, 61572651156383, 63771674412287, 65970697666967, 68169720924167, 74766790688867, 79164837200927, 83562883712027, 87960930223163, 92358976733483, 96757023247427, 101155069756823, 105553116266999, 109951162779203, 114349209290003, 118747255800179, 123145302311783, 127543348823027, 131941395333479, 136339441846019, 149533581378263, 158329674402959, 167125767424739, 175921860444599, 184717953466703, 193514046490343, 202310139514283, 211106232536699, 219902325558107, 228698418578879, 237494511600287, 246290604623279, 255086697645023, 263882790666959, 272678883689987, 299067162755363, 316659348799919, 334251534845303, 351843720890723, 369435906934019, 387028092977819, 404620279022447, 422212465067447, 439804651111103, 457396837157483, 474989023199423, 492581209246163, 510173395291199, 527765581341227, 545357767379483, 598134325510343, 633318697599023, 668503069688723, 703687441776707, 738871813866287, 774056185954967, 809240558043419, 844424930134187, 879609302222207, 914793674313899, 949978046398607, 985162418489267, 1020346790579903, 1055531162666507, 1090715534754863};
    
	/* Read the input file and store the reads in binary. */
	readInput(binReads, argv, fileFormatFlag, totalBases, numberOfReads, maxReadLength, goodRead, readLocation, currentMemory, peakMemory);

	/* Find the coverage and the number of read splits. */
	coverage = static_cast<int>(floor(totalBases/genomeLength));
	numOfReadSplits = static_cast<uint64_t>(floor(coverage/70));
	if (numOfReadSplits > 1)
		readsInEachSplit = static_cast<int>(floor(numberOfReads/numOfReadSplits));
	else 
	{
		numOfReadSplits = 1;
		readsInEachSplit = numberOfReads;
	}
	totalBasesPerSplit = totalBases/numOfReadSplits;
	endLocation = numberOfReads;
	
	/* Set the initial size of the hash table to 3 times the size of the genome. */
	for (n=0; n<450; ++n) 
		if (hashTableSize[n] > (genomeLength*3)) 
		{
			hashTableSizeIndex = n;
			break;
		}
	currHashTabSize = hashTableSize[hashTableSizeIndex];  //Store the size of the hash table.

	/* Create the counters array. */
	assert(counters = new uint8_t[8*currHashTabSize]);

	currentMemory += sizeof(uint8_t)*8*currHashTabSize;
	if (currentMemory > peakMemory)
	{	previousPeak = peakMemory;
		peakMemory = currentMemory;
#ifdef VERBOSE
		cout << "New peak Memory = " << (peakMemory/1048576) << " MB. counters array added " << (peakMemory-previousPeak)/1048576 << " MB." << endl << flush;
#endif
	}
    
	/* Create the witness array. */
	assert(witness = new uint64_t[currHashTabSize]);

	currentMemory += sizeof(uint64_t)*currHashTabSize;
	if (currentMemory > peakMemory)
	{	previousPeak = peakMemory;
		peakMemory = currentMemory;
#ifdef VERBOSE
		cout << "New peak Memory = " << (peakMemory/1048576) << " MB. witness array added " << (peakMemory-previousPeak)/1048576 << " MB." << endl << flush;
#endif
	}

	/* Find the witness lengths needed for optimal correction. */
    computeWitLength(maxReadLength, readsInEachSplit, genomeLength, witLength_small, witLength_large);
	
	if(witLength_small == witLength_large)
		witLength_small--;
		
	/* Calculate all the possible values for the reverse complement of an 8 bit value */
	uint8_t* allRevCompOf8Bit; 
	allRevCompOf8Bit = new uint8_t[256];
	revComplOf8Bit(allRevCompOf8Bit);

	/* Output information to the screen. */
	cout << "Genome length = " << genomeLength << endl << flush;
	cout << "Coverage = " << coverage << endl << flush;
	cout << "Maximum read length = " << maxReadLength << endl << flush;
	cout << "Number of reads = " << numberOfReads << endl << flush;

#ifdef VERBOSE
	cout << "Number of splits = " << numOfReadSplits << endl << flush;
	cout << "Number of reads in each split = " << readsInEachSplit << endl << flush;
	cout << "Initial hash table size = " << currHashTabSize << endl << flush;
//	cout << "Threshold = " << threshold << endl << flush;
	cout << "witLength_small (wm) = " << witLength_small << endl << flush;
	cout << "witLength_large (wM) = " << witLength_large << endl << flush;
#endif
    
    int tid = 0;
    #pragma omp parallel
    {
        tid = omp_get_num_threads();
    }
    cout << "Number of threads = " << tid << endl;    
    cout << "correcting reads..." << flush;
    
    /* Process the reads to compute the witness pool, counters, and correct errors. */
    for (z=0; z<numOfReadSplits; ++z) 
	{	
		/* Set the number of reads to be processed by this split. */
		startLocation = endLocation-1;
		if (z == (numOfReadSplits-1)) 
			endLocation = 0;
		else
			endLocation -= readsInEachSplit;
        
#ifdef VERBOSE		
		cout << endl << "------------------------------ Split #" << z+1 << " -------------------------------" << endl << flush;
#endif
        
		/* Reset iteration count after each split. */
		iteration=0;
		
		/* Correcting the errors in the reads for as many iterations as required. */
		do {
			/* Set the witness length for each iteration. */
			iteration++;
			if(iteration == 1)		 witnessLength = witLength_small+1;
			else if (iteration == 2) witnessLength = witLength_small+1;
			else if (iteration == 3) witnessLength = witLength_large+1;
			else if (iteration == 4) witnessLength = witLength_large+1;
			else if (iteration == 5) witnessLength = witLength_small-1;
			else if (iteration == 6) witnessLength = witLength_small-1;
			else if (iteration == 7) witnessLength = witLength_large-1;
			else if (iteration == 8) witnessLength = witLength_large-1;
			else if (iteration == 9) witnessLength = witLength_large+1;
			else witnessLength = witLength_large+1;
			
			witnessLengthInBits = witnessLength*2;
			
			/* Output information for this iteration to the screen. */
            
#ifdef VERBOSE
			cout << endl << "================== Iteration #" << iteration << "  Witness Length = " << witnessLength << " ==================" << endl << flush;
#endif
            
			/* Find the witnesses and calculate the counters for the before and after bases. */
			if (witnessLength != previousWitnessLength)
			{
                
				computeThreshold(maxReadLength, readsInEachSplit, genomeLength, witnessLength, threshold);
#ifdef VERBOSE
				cout << "threshold = " << threshold << endl;
				cout << "Starting to build the witnesses and counters." << endl << flush;
#endif
       			/* Compute threshold such that values larger are correct, smaller are wrong. */

				time(&t_witnessStart);

				buildWitnessesAndCounters(binReads, witness, counters, startLocation, endLocation, allRevCompOf8Bit, witnessLengthInBits, hashTableSize, hashTableSizeIndex, currHashTabSize, maxReadLength, goodRead, readLocation, currentMemory, peakMemory);

				time(&t_witnessEnd);
				witnessTime += difftime(t_witnessEnd,t_witnessStart);
			
				/* Find the number of witnesses with a value above the threshold for the before and after bases. */
				numOfWitAbovThresholdForAfterBase=0;
				numOfWitAbovThresholdForBeforeBase=0;
				for (x=0; x < currHashTabSize; ++x) 
				{
					for (y=0; y<4; ++y) 
						if (counters[8*x+y] > threshold) 
						{
							numOfWitAbovThresholdForAfterBase++;
						}
					for (y=4; y<8; ++y) 
						if (counters[8*x+y] > threshold) 
						{
							numOfWitAbovThresholdForBeforeBase++;
						}
				}
#ifdef VERBOSE
				cout << "Witnesses above threshold in before base = " << numOfWitAbovThresholdForBeforeBase << endl << flush;
				cout << "Witnesses above threshold in after base = " << numOfWitAbovThresholdForAfterBase << endl << flush;
				cout << "Calculated the witness counters in " << difftime(t_witnessEnd,t_witnessStart) << " seconds." << endl << flush;
#endif
                
			}
			
			previousWitnessLength = witnessLength;
			
			/* Correct the reads based on the witness pool and counters. */
			time(&t_correctStart);
            
#ifdef VERBOSE
			cout << endl << "Starting to correct reads." << endl << flush;
#endif
            
			correctErrors(binReads, witness, counters, startLocation, endLocation, witnessLengthInBits, allRevCompOf8Bit, threshold, numOfPosChanged, currHashTabSize, maxReadLength, goodRead, readLocation, hashTableSize, hashTableSizeIndex);
			
			time(&t_correctEnd);
			correctingTime += difftime(t_correctEnd,t_correctStart);
            
#ifdef VERBOSE
			cout << "Corrected the reads in " << difftime(t_correctEnd,t_correctStart) << " seconds." << endl << flush;
#endif
			
		}while (iteration < 10);
	}

    cout << "done correcting reads" << endl << flush;
    
	/* Deleting the arrays that are no longer needed.  */
	delete [] counters;                                                 
	delete [] witness;
	delete [] allRevCompOf8Bit;

    cout << "writing back the reads..." << flush;
    
	/* Write the corrected reads to a fastq or fasta file. */
	writeOutput(binReads, argv, fileFormatFlag, goodRead, readLocation);

    cout << "done" << endl << flush;

	/* Delete the rest of the arrays. */
	delete [] binReads;
	delete [] goodRead;
	delete [] readLocation;
	
	time(&t_end);
    
#ifdef VERBOSE
	cout << endl << "============================== Finished ===============================" << endl << flush;
#endif
	cout << "Total time building witnesses and counters = " << witnessTime << " seconds." << endl << flush;
	cout << "Total time correcting = " << correctingTime << " seconds." << endl << flush;
	cout << "Total time = " << difftime(t_end,t_start) << " seconds." << endl << flush;
	cout << "Peak Memory = " << (peakMemory/1048576) << " MB." << endl << flush;
    
	return 0;
}
