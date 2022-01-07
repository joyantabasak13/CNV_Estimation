#include <bits/stdc++.h>
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

#define _USE_MATH_DEFINES

#define weightFile "GCvsOverlappingReads.txt"
#define contigGcFile "gcValue.txt"
#include <math.h>
#include "cgalHeader.h"

#define MAX_REC_LEN 1024

///copy number
#define MAX_COPY_NUMBER 10

class multiMapReadData{
    public:
        int contigNo;
        int insertSize;
        double prob;
};

map< long long,vector<multiMapReadData> > multiMapReads;///index= Read index, value=contig wise mapping data of that read
map<int,vector< long long> > readsInContig;///contigno index, value= read Index

int noContigs=0;
int noReads=0;
long int contigLength=0;
int maxReadLength=0;

long int fileLineCount=0;
long int totalContigLength=0;
long int tcl = 0;

//Important for contig properties

vector<char*> contigs;
vector<char*> contigNames;
vector<long int> contigLengths;

///unique mapped read
vector<long int> uniqueMappedReadCount;
double likelihoodOriginal;

vector<int> contigCopyNumber;

vector<vector<int> >ContigVsCopyVsProb;///row= a contig column = copy no. val=likelihood
map< int, vector<double> >insertVsContigWeightedLen;///key=insert column=contig val=weighted len of a contig for the insertSize





FILE *contigFile;
FILE *mapFile;
FILE *summaryFile;
FILE *outFile;

char *contigFileName;
char *mapFileName;

long int *insertCounts;
int maxInsertSize=0;
int MAX_INSERT_SIZE=100000;
double insertSizeMean;
double insertSizeVar;


long int errorTypes[5][5];
long int baseCounts[5];
long int *errorPos;
long int *inPos;
long int *inLengths;
long int *delPos;
long int *delLengths;
long int *readLengths;

long int *effectiveLengths;
///changes
double *weightedEffectiveLengths;
///changes


double errorTypeProbs[5][5];
double baseErrorRates[5];
double *errorPosDist;
double *inPosDist;
double *inLengthDist;
double *delPosDist;
double *delLengthDist;
double *insertLengthDist;

double *noErrorProbs;
int tmpCount=0;
int toAlign;

long int erroredReads=0;
long int uniqueMappedReads=0;
long int discardedReads=0;
long int totalCount,unCount;

char tempCigar[500], tempMD[500];
char noErrorCigar[500], noErrorMD[500];


///GC Bias....
int windowSize=500;
int gcCount=0;
vector<vector<double> > gcContentArr;
int gcArrIndex=0;

ofstream gcWriter;

int cn=0;

vector<double>gcVSweight;///index=GC in integer form
                    ///value = weight or number of overlapping reads


vector < vector<int> > overlappingReadsInContig;///index = contigNo,pos value = number of overlapping reads


long int findBucket(double gcVal){
    ///returns the bucket/index of the fractional gc value
    ///using this index we can get the weight from weight array
    ///the returned index must start from 0.
    return int(gcVal*100);
}



//This function prints the length of all contigs
void print_cntg_length(){

  int max_length = 0,mark;
  for(int i=0;i<contigs.size();i++){
      cout<<"Length of "<<i<<"th"<<"is "<<contigLengths[i]<<endl;
      if(contigLengths[i]>max_length){
        max_length = contigLengths[i];
        mark=i;
      }

  }
  cout<<mark<<" th"<< " contig length is maximum. Max length is "<<max_length<<endl;

}
void updateGCArr()
{
    cout<<"The number of total contigs is "<<contigs.size()<<endl;

    ///Outer loop for accessing each contig
    for(int k=0; k<contigs.size(); k++)
    {
        int currentPosition=0,dividingWindow=0;
        gcCount=0;
        vector<double> gcSingleContig; ///GC content array for each contig

        int lowerMark=0,upperMark=windowSize/2+1; ///Two markers for GC counter
        ///This loop for counting GC in first window
        for(int i=0; i<=windowSize/2; i++)
        {
            if(i<contigLengths[k])
            {
                if(contigs[k][i]=='G' || contigs[k][i]=='C')
                {
                    gcCount++;
                }
                dividingWindow++;
            }
        }
        currentPosition++;
        gcSingleContig.push_back((gcCount*1.0)/dividingWindow); ///GC content is pushed into the array
        ///Then the two markers update the GC counter according to upper and lower content

        while(currentPosition<contigLengths[k])
        {
            if(upperMark<contigLengths[k])
            {
                if((contigs[k][upperMark]=='G' || contigs[k][upperMark]=='C') )
                {
                    gcCount++;
                }
                upperMark++;
            }

            if((currentPosition-lowerMark)>windowSize/2)
            {
                if((contigs[k][lowerMark]=='G' || contigs[k][lowerMark]=='C'))
                {
                    gcCount--;
                }
                lowerMark++;
            }

            gcSingleContig.push_back((gcCount*1.0)/(upperMark-lowerMark)); ///GC content is pushed into the array for each window

            currentPosition++;

        }

        for(int p=0; p<gcSingleContig.size(); p++)
        {
            gcWriter<<gcSingleContig[p]<<"\n";
            fileLineCount++;
        }
        gcContentArr.push_back(gcSingleContig); ///Array for each contig is pushed


    }

}

///GC Bias edit ends....

void initInsertCounts(int max)
{
    maxInsertSize=max;
    insertCounts=new long int[maxInsertSize];
    for(int i=0; i<maxInsertSize; i++)
    {
        insertCounts[i]=1;
    }
}

void updateInsertCounts(int index)
{

    if(index<=0)
        return;
    if(index<maxInsertSize)
    {
        insertCounts[index]++;
    }
    else
    {

        if(index>MAX_INSERT_SIZE)
        {
            discardedReads++;
            return;
        }
        int tempInsertSize=max(maxInsertSize*2,index);
        long int *tempCounts=new long int[maxInsertSize];
        for(int i=0; i<maxInsertSize; i++)
        {
            tempCounts[i]=insertCounts[i];
        }
        insertCounts=new long int[tempInsertSize];
        for(int i=0; i<maxInsertSize; i++)
        {
            insertCounts[i]=tempCounts[i];
        }
        for(int i=maxInsertSize; i<tempInsertSize; i++)
        {
            insertCounts[i]=1;
        }

        insertCounts[index]++;
        maxInsertSize=tempInsertSize;
        delete []tempCounts;

    }

}

void initErrorTypes(int readLength)
{
    for(int i=0; i<5; i++)
        for(int j=0; j<5; j++)
            errorTypes[i][j]=1;

    for(int i=0; i<5; i++)
        baseCounts[i]=1;

    errorPos=new long int[readLength];
    inPos=new long int[readLength];
    inLengths=new long int[readLength];
    delPos=new long int[readLength];
    delLengths=new long int[readLength];
    readLengths=new long int[readLength];

    for(int i=0; i<readLength; i++)
    {
        errorPos[i]=1;
        inPos[i]=1;
        inLengths[i]=1;
        delPos[i]=1;
        delLengths[i]=1;
        readLengths[i]=0;
    }
}


int getLength(char *read)
{

    int i=0;
    while(read[i])
    {
        if(read[i]=='A')
            baseCounts[0]++;
        else if(read[i]=='C')
            baseCounts[1]++;
        else if(read[i]=='G')
            baseCounts[2]++;
        else if(read[i]=='T')
            baseCounts[3]++;
        else
            baseCounts[4]++;

        i++;
    }

    return i;
}


void processErrorTypes(char *cigar, char *md, char *read, int strandNo)
{

    int readLength=getLength(read);
    readLengths[readLength-1]++;

    if(strcmp(md,noErrorCigar)!=0)
        erroredReads++;
    else
        return;


    int mdLength=strlen(md)-5;
    int tempLength=0;

    char *temp;
    int index=0,totalLength=0;


    int curIndex=0;
    int *inserts=new int[readLength];

    for(int i=0; i<readLength; i++)
    {
        inserts[i]=0;
    }

    int cigarLength=strlen(cigar);
//	char *tempCigar=new char[cigarLength];
    char cigarChar;

    strcpy(tempCigar,cigar);


    temp=strtok(tempCigar,"IDMS^\t\n ");

    while(temp!=NULL)
    {

        tempLength=atoi(temp);
        totalLength+=strlen(temp);
        cigarChar=cigar[totalLength];

        if(cigarChar=='M')
        {
            index+=tempLength;
            curIndex+=tempLength;
        }
        else if(cigarChar=='I' || cigarChar=='S')
        {
            if(strandNo==0)
            {
                inPos[index]++;
                inLengths[tempLength-1]++;

            }
            else
            {

                inPos[readLength-index-1]++;
                inLengths[tempLength-1]++;
            }

            inserts[curIndex]=tempLength;

            index+=tempLength;
        }
        else if(cigarChar=='D' )
        {
            if(strandNo==0)
            {
                delPos[index]++;
                delLengths[tempLength-1]++;

            }
            else
            {

                delPos[readLength-index-1]++;
                delLengths[tempLength-1]++;
            }
        }
        totalLength++;
        temp=strtok(NULL,"IDMS^\t\n ");
    }


    strcpy(tempMD,md);

    strtok(tempMD,":");
    strtok(NULL,":");


    index=0,totalLength=0,tempLength=0;

    int f, t;

    while((temp=strtok(NULL,"ACGTN^\t\n "))!=NULL)
    {
        tempLength=strlen(temp);


        totalLength+=tempLength;


        if(totalLength<mdLength)
        {
            char from=md[5+totalLength];


            if(from=='^')
            {
                totalLength++;
                index+=atoi(temp);
                for(int i=totalLength; i<mdLength; i++)
                {
                    from=md[5+totalLength];
                    if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
                        totalLength++;
                    else
                        break;

                }
            }
            else if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
            {
                totalLength++;
                index+=atoi(temp)+1;



                curIndex=0;
                for(int i=0; i<index; i++)
                {
                    curIndex+=inserts[i];
                }
                char to=read[index-1+curIndex];

                if(strandNo==0)
                    errorPos[index-1+curIndex]++;
                else
                    errorPos[readLength-index-curIndex]++;


                switch(from)
                {
                case 'A':
                    f=0;
                    break;
                case 'C':
                    f=1;
                    break;
                case 'G':
                    f=2;
                    break;
                case 'T':
                    f=3;
                    break;
                default:
                    f=4;
                }

                switch(to)
                {
                case 'A':
                    t=0;
                    break;
                case 'C':
                    t=1;
                    break;
                case 'G':
                    t=2;
                    break;
                case 'T':
                    t=3;
                    break;
                default:
                    t=4;
                }

                if(f==t)
                {

                }
                else

                    errorTypes[f][t]++;

            }
            else
                break;
        }

    }
    delete []inserts;

}

void computeProbabilitesForGC()
{

	int errorCount=0;

	for(int i=0;i<5;i++)
	{
		errorCount=0;
		for(int j=0;j<5;j++)
		{
			errorCount+=errorTypes[i][j];
		}
		for(int j=0;j<5;j++)
		{
			errorTypeProbs[i][j]=(double)errorTypes[i][j]/errorCount;
		}

		baseErrorRates[i]=errorCount/(double)baseCounts[i];
	}

	double sum=0;
	for(int i=0;i<4;i++)
		sum+=baseErrorRates[i];

	for(int i=0;i<4;i++)
	{
		baseErrorRates[i]=4*baseErrorRates[i]/sum;
	}

	baseErrorRates[4]=1;

	for(int i=maxReadLength-1;i>0;i--)
	{
		readLengths[i-1]=readLengths[i]+readLengths[i-1];
	}

	errorPosDist=new double[maxReadLength];

	for(int i=0;i<maxReadLength;i++)
	{
		errorPosDist[i]=(double)errorPos[i]/readLengths[i];
	}

	inPosDist=new double[maxReadLength];

	for(int i=0;i<maxReadLength;i++)
	{
		inPosDist[i]=(double)inPos[i]/readLengths[i];
	}

	inLengthDist=new double[maxReadLength];


	int inCount=0;

	for(int i=0;i<maxReadLength;i++)
	{
		inCount+=inLengths[i];
	}

	for(int i=0;i<maxReadLength;i++)
	{
		inLengthDist[i]=(double)inLengths[i]/inCount;
	}

	delPosDist=new double[maxReadLength];


	for(int i=0;i<maxReadLength;i++)
	{
		delPosDist[i]=(double)delPos[i]/readLengths[i];
	}

	delLengthDist=new double[maxReadLength];

	int delCount=0;

	for(int i=0;i<maxReadLength;i++)
	{
		delCount+=delLengths[i];
	}

	for(int i=0;i<maxReadLength;i++)
	{
		delLengthDist[i]=(double)delLengths[i]/delCount;
	}


	insertLengthDist=new double[maxInsertSize];

	long int insCount=discardedReads;

	sum=0;

	for(int i=0;i<maxInsertSize;i++)
	{
		insCount+=insertCounts[i];
		sum+=i*insertCounts[i];
	}
	insertSizeMean=sum/insCount;

	sum=0;

	for(int i=0;i<maxInsertSize;i++)
	{
		insertLengthDist[i]=(double)insertCounts[i]/insCount;

		sum+=insertCounts[i]*(insertSizeMean-i)*(insertSizeMean-i);
	}

	insertSizeVar=sum/insCount;


	noErrorProbs=new double[maxReadLength];

	double noErrorProb=1.0;

	for(int i=0;i<maxReadLength;i++)
	{
		noErrorProb*=(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
		noErrorProbs[i]=noErrorProb;
	}

	effectiveLengths=new long int[maxInsertSize];

	for(int i=0;i<maxInsertSize;i++)
	{
		effectiveLengths[i]=-1;
	}

	long int totalContigLength=0;
	for(int i=0;i<contigLengths.size();i++)
	{
		totalContigLength+=contigLengths[i];
	}
	effectiveLengths[0]=totalContigLength;

}

void computeProbabilites()
{


    ///changes
    weightedEffectiveLengths=new double[maxInsertSize];
    ///changes done


    for(int i=0; i<maxInsertSize; i++)
    {
        ///changes
        weightedEffectiveLengths[i]=-1;
        ///changes done
    }

    double totalContigSiteWeight=0;///no confusion

    /*ifstream inFile;

    inFile.open(weightFile);
    if (!inFile) {
        cout << "Unable to open file "<<weightFile<<endl;
        exit(1); // terminate with error
    }
    ///take these file values in gcVSweight array
    double z;
    double bucket;
    while(1){
        if(inFile >> bucket){///1st value of a line.The bucket
            inFile >> z;///2nd value of the line. The main weight(overlapping reads)
            gcVSweight.push_back(z);
        }
        else break;
    }
    inFile.close();
    ///done taking values into array*/
    long int index;
    for(int i=0; i<contigLengths.size(); i++)
    {
        ///changing
        long int currentContigLength = contigLengths[i];
        for(int j=0;j<currentContigLength;j++){
            ///find the position in the gcContentArr where the contig has started with its gc contents.
            ///From there get all the weights of that contig.

            index = findBucket(gcContentArr[i][j]);
            totalContigSiteWeight+=gcVSweight[index];/**each weight value found from file for this contig**/
        }
        /// done changing
    }
    ///changes
    weightedEffectiveLengths[0]=totalContigSiteWeight;
    ///changes done
}

double dnorm(double x,double mean, double variance)
{
    double val=1/sqrt(M_PI*2*variance);
    val*=exp(-((x-mean)*(x-mean))/(2*variance));
    return val;
}

void processMapping(char *line)
{

    char * temp;
    char *qname, *rname, *mapq;
    int	pos,flag,contigNo;
    char * cigar, * readString; // * md, *nhstring;

    char md[500];
    char nhstring[500];

    int nh;

    int strandNo=0;


    qname=strtok(line,"\t");

    temp=strtok(NULL,"\t");
    flag=atoi(temp);


    strandNo=(flag&16)>>4;
///contig no
    temp=strtok(NULL,"\t");
    contigNo=atoi(temp);

    temp=strtok(NULL,"\t");
    pos=atoi(temp);


    cigar=strtok(NULL,"\t");


    temp=strtok(NULL,"\t");

    readString=strtok(NULL,"\t");

    int insertSize=atoi(temp);

    while((temp=strtok(NULL,"\t\n"))!=NULL)
    {
        if(temp[0]=='M' && temp[1]=='D')
        {
            strcpy(md,temp);
        }
        else if(temp[0]=='I' && temp[1]=='H')
        {
            strcpy(nhstring,(temp+5));
            nh=atoi(nhstring) ;
        }

    }


    if(nh==1 && md[5]!='^')
    {
        updateInsertCounts(insertSize);
        processErrorTypes(cigar,md,readString,strandNo);
        uniqueMappedReads++;
    }


}

long int getEffectiveLength(int insertSize)
{
    if(insertSize<0)
        return effectiveLengths[0];

    if(insertSize>=maxInsertSize)
    {
        long int effectiveLength=0;
        for(int i=0; i<contigLengths.size(); i++)
        {
            if(contigLengths[i]>=insertSize)
                effectiveLength+=(contigLengths[i]-insertSize+1);
        }
        return effectiveLength;

    }
    if(effectiveLengths[insertSize]==-1)
    {
        long int effectiveLength=0;
        for(int i=0; i<contigLengths.size(); i++)
        {
            if(contigLengths[i]>=insertSize)
                effectiveLength+=(contigLengths[i]-insertSize+1);
        }
        effectiveLengths[insertSize]=effectiveLength;
    }
    return effectiveLengths[insertSize];
}

///Effective length correction using GC bias Error Model.
///This function calculates effective length for a given insertion (fragment)...
///considering the weight factor of each position.
///Here Weight of a position = overlapping reads of that position.
double getWeightedEffectiveLength(int insertSize)
{
    if(insertSize<0)
        return weightedEffectiveLengths[0];

    if(insertSize>=maxInsertSize)
    {
        double weightedEffectiveLength=0;
        long int index;
        for(int i=0; i<contigLengths.size(); i++)
        {
            if(contigLengths[i]>=insertSize)
            for(int j=0;j<=contigLengths[i]-insertSize;j++){
                index = findBucket(gcContentArr[i][j]);
                weightedEffectiveLength+= gcVSweight[index];
            }
        }
        return weightedEffectiveLength;

    }
    if(weightedEffectiveLengths[insertSize]==-1)
    {


        double weightedEffectiveLength=0;
        long int index;
        for(int i=0; i<contigLengths.size(); i++)
        {
            double contigWeight=0.0;
            if(contigLengths[i]>=insertSize){
                for(int j=0;j<=contigLengths[i]-insertSize;j++){
                    index = findBucket(gcContentArr[i][j]);
                    weightedEffectiveLength+= gcVSweight[index];
                    contigWeight+=gcVSweight[index];
                }
            }
            insertVsContigWeightedLen[insertSize].push_back(contigWeight);
        }
        weightedEffectiveLengths[insertSize]=weightedEffectiveLength;
    }
    return weightedEffectiveLengths[insertSize];
}


double computeErrorProb(char *cigar, char *md, char *read, int strandNo)
{

    int readLength=strlen(read);



    double errorProb=noErrorProbs[readLength-1];

    if(md[5]=='^')
        return errorProb;


    char tempMD[1000], tempCigar[1000];

    int mdLength=strlen(md)-5;
    int tempLength=0;

    char *temp;
    int index=0,totalLength=0;


    int curIndex=0;
    int *inserts=new int[readLength];

    for(int i=0; i<readLength; i++)
    {
        inserts[i]=0;
    }

    int cigarLength=strlen(cigar);
    char cigarChar;

    strcpy(tempCigar,cigar);

    temp=strtok(tempCigar,"IDM^\t\n ");

    while(temp!=NULL)
    {

        tempLength=atoi(temp);
        totalLength+=strlen(temp);
        cigarChar=cigar[totalLength];

        if(cigarChar=='M')
        {
            index+=tempLength;
            curIndex+=tempLength;
        }
        else if(cigarChar=='I')
        {
            int i;
            if(strandNo==0)
            {
                //look up insert probs
                i=index;

            }
            else
            {
                i=readLength-index-1;
            }

            errorProb=errorProb*inPosDist[i]*inLengthDist[tempLength-1]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);

            inserts[curIndex]=tempLength;

            index+=tempLength;
        }
        else if(cigarChar=='D')
        {
            int i;
            if(strandNo==0)
            {
                i=index;
                //	look up delete probs
            }
            else
            {
                i=readLength-index-1;
            }

            errorProb=errorProb*delPosDist[i]*delLengthDist[tempLength-1]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
        }
        totalLength++;
        temp=strtok(NULL,"IDM^\t\n ");
    }


    strcpy(tempMD,md);

    strtok(tempMD,":");
    strtok(NULL,":");

    index=0,totalLength=0,tempLength=0;

    int f, t;

    while((temp=strtok(NULL,"ACGTN^\t\n "))!=NULL)
    {
        tempLength=strlen(temp);

        totalLength+=tempLength;

        if(totalLength<mdLength)
        {
            char from=md[5+totalLength];

            if(from=='^')
            {
                totalLength++;
                index+=atoi(temp);
                for(int i=totalLength; i<mdLength; i++)
                {
                    from=md[5+totalLength];
                    if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
                        totalLength++;
                    else
                        break;
                }
            }
            else if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
            {
                totalLength++;
                index+=atoi(temp)+1;


                curIndex=0;
                for(int i=0; i<index; i++)
                {
                    curIndex+=inserts[i];
                }
                char to=read[index-1+curIndex];

                int i;
                if(strandNo==0)
                    i=index-1+curIndex;
                else
                    i=readLength-index-curIndex;


                errorProb=errorProb*errorPosDist[i]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);


                switch(from)
                {
                case 'A':
                    f=0;
                    break;
                case 'C':
                    f=1;
                    break;
                case 'G':
                    f=2;
                    break;
                case 'T':
                    f=3;
                    break;
                default:
                    f=4;
                }

                switch(to)
                {
                case 'A':
                    t=0;
                    break;
                case 'C':
                    t=1;
                    break;
                case 'G':
                    t=2;
                    break;
                case 'T':
                    t=3;
                    break;
                default:
                    t=4;
                }

                if(f==t)
                {


                }
                else
                {
                    //errorTypeProb
                    errorProb*=baseErrorRates[f]*errorTypeProbs[f][t];
                }
            }
            else
                break;

        }

    }

    delete []inserts;
    return errorProb;
}


double computeLikelihoodForGC(char *file)
{

    mapFile=fopen(file, "r");
    char *line1= new char[MAX_REC_LEN];
    char *line2= new char[MAX_REC_LEN];

    char *qname1,*qname2,preqname1[500],preqname2[500];

    int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);

    double sum=0.0;
    double logsum=0.0;

    char * temp;
    char *rname1, *rname2;
    int	pos1,pos2,flag,strandNo1, strandNo2, insertSize1, insertSize2,contigNo1,contigNo2;
    char *cigar1, *cigar2, *readString1, *readString2, md1[1000], md2[1000];


    double insertSizeProb;

    double errorProb1, errorProb2;

    int noUnmappedReads=0;
    int noUniqueMappedReads=0;
    int noMultMappedReads=0;

    preqname1[0]=preqname2[0]='*';
    preqname1[1]=preqname2[1]=0;

    ///edited from raihan vai
    char *contigField1, *contigField2;
    vector<info> multiMapProb1, multiMapProb2;
	vector<bs> bsCollection1,bsCollection2;
	map<long int, string> cigarMap1,cigarMap2;
	map<long int, string> readMap1,readMap2;
	int fileCount=0;


    ///done


    int it=0;

    while(fgets(line1, MAX_FILE_READ, mapFile)!=NULL)
    {

        if(line1[0]=='@')
            continue;
//????
        if(fgets(line2, MAX_FILE_READ, mapFile)==NULL)
            break;



        qname1=strtok(line1,"\t");
        temp=strtok(NULL,"\t");
        flag=atoi(temp);


        strandNo1=(flag&16)>>4;

        ///changing for contigNo
        temp=strtok(NULL,"\t");
        contigField1=temp;
        contigNo1=atoi(temp);


        temp=strtok(NULL,"\t");
        pos1=atoi(temp);


        cigar1=strtok(NULL,"\t");


        temp=strtok(NULL,"\t");



        insertSize1=atoi(temp);


        readString1=strtok(NULL,"\t");


        while((temp=strtok(NULL,"\t\n"))!=NULL)
        {
            if(temp[0]=='M' && temp[1]=='D')
            {
                strcpy(md1,temp);
            }
        }

//		cout<<insertSize1<<" "<<cigar1<<" "<<md1<<endl;


//second of the pair

        qname2=strtok(line2,"\t");
        temp=strtok(NULL,"\t");
        flag=atoi(temp);


        strandNo2=(flag&16)>>4;

        ///changing for contigNo
        temp=strtok(NULL,"\t");
        contigField2=temp;
        contigNo2=atoi(temp);

        temp=strtok(NULL,"\t");
        pos2=atoi(temp);



        cigar2=strtok(NULL,"\t");


        temp=strtok(NULL,"\t");

        insertSize2=atoi(temp);

        readString2=strtok(NULL,"\t");


        while((temp=strtok(NULL,"\t\n"))!=NULL)
        {
            if(temp[0]=='M' && temp[1]=='D')
            {
                strcpy(md2,temp);
            }
        }

//		cout<<insertSize2<<" "<<cigar2<<" "<<md2<<endl;



        int insertSize=max(insertSize1, insertSize2);


        insertSizeProb=0;

        if(insertSize>=0 && insertSize<maxInsertSize)
        {
            insertSizeProb=insertLengthDist[insertSize];
        }

        if(insertSizeProb==0)
        {
            insertSizeProb=1/(double)uniqueMappedReads;
        }


        errorProb1=computeErrorProb(cigar1,md1,readString1,strandNo1);


        errorProb2=computeErrorProb(cigar2,md2,readString2,strandNo2);


        long int totalEffectiveLength=getEffectiveLength(insertSize);

        double prob=(1/(double)(totalEffectiveLength))*insertSizeProb*errorProb1*errorProb2;



///changes
       /* int index = findBucket(gcContentArr[contigNo1][pos1]);
        double currentWeight = gcVSweight[index];

        double totalWeightedEffectiveLength=getWeightedEffectiveLength(insertSize);
        double prob=(currentWeight/(double)(totalWeightedEffectiveLength))*insertSizeProb*errorProb1*errorProb2;*/
///changes done

//		cout<<errorProb1<<" "<<errorProb2<<" "<<insertSizeProb<<" "<<prob<<endl;

        if(strcmp(qname1,preqname1)==0 && strcmp(qname2,preqname2)==0)
        {
            sum+=prob;

            ///added from raihan vai
            createInfo(insertSize,qname1, errorProb1, pos1, contigField1, readString1, multiMapProb1);
			createInfo(insertSize,qname2, errorProb2, pos2, contigField2, readString2, multiMapProb2);
			///done

        }
        else if(strcmp("*",preqname1)!=0 && strcmp("*",preqname2)!=0)
        {
            ///edited from raihan vai
            writeInfoToFile(multiMapProb1, cigarMap1, readMap1);
			writeInfoToFile(multiMapProb2, cigarMap2, readMap2);

			multiMapProb1.clear();
			multiMapProb2.clear();


            createInfo(insertSize,qname1, errorProb1, pos1, contigField1, readString1, multiMapProb1);
			createInfo(insertSize,qname2, errorProb2, pos2, contigField2, readString2, multiMapProb2);

			///done
        }
        else
        {
            ///edited from raihan vai
            createInfo(insertSize,qname1, errorProb1, pos1, contigField1, readString1, multiMapProb1);
			createInfo(insertSize,qname2, errorProb2, pos2, contigField2, readString2, multiMapProb2);
			///done

        }

        strcpy(preqname1,qname1);
        strcpy(preqname2,qname2);
        it++;


    }

    fclose(mapFile);

    delete[] line1;
    delete[] line2;

    return 0;
}

double computeLikelihood(char *file)
{
    //int multiCounter=0,totalMultiCounter=0,prevIh=1;
    int readsPerContig[1000000]={0};
    int cnt=0;
    int pq=50;
    mapFile=fopen(file, "r");
    long long qIndex=0;
    char *line1= new char[MAX_REC_LEN];
    char *line2= new char[MAX_REC_LEN];

    char *qname1,*qname2,preqname1[500],preqname2[500];

    int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);

    double sum=0.0;
    double logsum=0.0;

    char * temp;
    char *rname1, *rname2;
    int	pos1,pos2,flag,strandNo1, strandNo2, insertSize1, insertSize2,contigNo1,contigNo2;
    char *cigar1, *cigar2, *readString1, *readString2, md1[1000], md2[1000],ih1[100],ih2[100];
    bool multimapped1,multimapped2;

    double insertSizeProb;

    double errorProb1, errorProb2;

    int noUnmappedReads=0;
    int noUniqueMappedReads=0;
    int noMultMappedReads=0;

    preqname1[0]=preqname2[0]='*';
    preqname1[1]=preqname2[1]=0;

    ///edited from raihan vai
    char *contigField1, *contigField2;
    vector<info> multiMapProb1, multiMapProb2;
	vector<bs> bsCollection1,bsCollection2;
	map<long int, string> cigarMap1,cigarMap2;
	map<long int, string> readMap1,readMap2;
	int fileCount=0;


    ///done

    int it=0;

    ///clearing previous mappings
    readsInContig.clear();
    multiMapReads.clear();
    //uniqueMappedReadCount.clear();

    while(fgets(line1, MAX_FILE_READ, mapFile)!=NULL)
    {

        if(line1[0]=='@')
            continue;
//????
        if(fgets(line2, MAX_FILE_READ, mapFile)==NULL)
            break;



        qname1=strtok(line1,"\t");
        temp=strtok(NULL,"\t");
        flag=atoi(temp);


        strandNo1=(flag&16)>>4;

        ///changing for contigNo
        temp=strtok(NULL,"\t");
        contigField1=temp;
        contigNo1=atoi(temp);


        temp=strtok(NULL,"\t");
        pos1=atoi(temp);


        cigar1=strtok(NULL,"\t");


        temp=strtok(NULL,"\t");



        insertSize1=atoi(temp);


        readString1=strtok(NULL,"\t");


        while((temp=strtok(NULL,"\t\n"))!=NULL)
        {
            if(temp[0]=='M' && temp[1]=='D')
            {
                strcpy(md1,temp);
            }
            ///changed
            else if(temp[0]=='I' && temp[1]=='H')
            {
                temp=temp+5;
                strcpy(ih1,temp);
            }
        }

        readsPerContig[contigNo1]++;

        int t=atoi(ih1);
        if(t==1){
            uniqueMappedReadCount[contigNo1]++;
        }
        ///done

//second of the pair

        qname2=strtok(line2,"\t");
        temp=strtok(NULL,"\t");
        flag=atoi(temp);


        strandNo2=(flag&16)>>4;

        ///changing for contigNo
        temp=strtok(NULL,"\t");
        contigField2=temp;
        contigNo2=atoi(temp);

        temp=strtok(NULL,"\t");
        pos2=atoi(temp);



        cigar2=strtok(NULL,"\t");


        temp=strtok(NULL,"\t");

        insertSize2=atoi(temp);

        readString2=strtok(NULL,"\t");


        while((temp=strtok(NULL,"\t\n"))!=NULL)
        {
            if(temp[0]=='M' && temp[1]=='D')
            {
                strcpy(md2,temp);
            }
        }
        ///done


        int insertSize=max(insertSize1, insertSize2);


        insertSizeProb=0;

        if(insertSize>=0 && insertSize<maxInsertSize)
        {
            insertSizeProb=insertLengthDist[insertSize];
        }

        if(insertSizeProb==0)
        {
            insertSizeProb=1.0/(double)uniqueMappedReads;
        }


        errorProb1=computeErrorProb(cigar1,md1,readString1,strandNo1);


        errorProb2=computeErrorProb(cigar2,md2,readString2,strandNo2);


 //       long int totalEffectiveLength=getEffectiveLength(insertSize);

 //       double prob=(1/(double)(totalEffectiveLength))*insertSizeProb*errorProb1*errorProb2;

///changes
        int index = findBucket(gcContentArr[contigNo1][pos1]);
        double currentWeight = gcVSweight[index];

        double totalWeightedEffectiveLength=getWeightedEffectiveLength(insertSize);
        double prob=(currentWeight/(double)(totalWeightedEffectiveLength))*insertSizeProb*errorProb1*errorProb2;
///changes done

        if(strcmp(qname1,preqname1)==0 && strcmp(qname2,preqname2)==0)
        {
            sum+=prob;
            if(readsInContig.find(contigNo1)==readsInContig.end())
            {
                readsInContig[contigNo1].push_back(qIndex);
            }
            else
            {
                if(readsInContig[contigNo1].back()!=(qIndex))
                    readsInContig[contigNo1].push_back(qIndex);
            }

        }
        else if(strcmp("*",preqname1)!=0 && strcmp("*",preqname2)!=0)
        {
            if(sum<1e-320 || isnan(sum))
            {
                sum=1e-320;
            }
            logsum+=log(sum);

            sum=prob;
            if(t>1){
                qIndex++;
                readsInContig[contigNo1].push_back(qIndex);
            }

        }
        else
        {
            sum=prob;
            if(t>1){
                readsInContig[contigNo1].push_back(qIndex);
            }
        }

        ///changing
        if(t>1){

            multiMapReadData ob;
            ob.contigNo=contigNo1;
            ob.insertSize=insertSize;
            if(prob<1e-320 || isnan(prob))
            {
                prob=1e-320;
            }
            ob.prob=prob;/// correct
            multiMapReads[qIndex].push_back(ob);
        }
        ///done

        strcpy(preqname1,qname1);
        strcpy(preqname2,qname2);
        it++;


        if(isinf( logsum ))
        {
            cout<<it<<endl;
            exit(1);
        }


    }
    if(sum!=0)
        logsum+=log(sum);
    //cout<<"total qindex : "<<qIndex<<endl;

    fclose(mapFile);
    for(int i=0; i<contigs.size(); i++)
    {
        cout<<"Reads in contig "<<i<<" = "<<readsPerContig[i]<<endl;
    }

    map<int, vector<long long> >::iterator it1;
    for(it1=readsInContig.begin();it1!=readsInContig.end();it1++){///iterates over every contig and its set of reads
        cout<<"Multimapped reads in contig "<<it1->first<<" is "<<it1->second.size()<<endl;
    }

    return logsum;
}

/*double computeLikelihoodForCNV(int copyCount, double prevLikelihood, int contigNo)
{
    int cnt=0;
    double likelihoodEffective=prevLikelihood,likelihoodMulti=prevLikelihood;
    //double sumInsert=0;
    long int n=uniqueMappedReadCount[contigNo];
    cout<<"uniqe mapped reads "<<n<<endl;
    cout<<"Multi mapped reads "<<readsInContig[contigNo].size()<<endl;
    //cout<<"likeli before unique "<<prevLikelihood<<"contig "<<contigNo<<endl;
    double likeli=prevLikelihood - n*log(1) + n*log(copyCount);///uniquely mapped handled
    cout<<"likeli after unique "<<likeli<<endl;
    double extraLen;///extra length is the amount of
                ///increase of effective length due to increase of a copy number
    double prevEffective,newEffective;
    for(int i=0;i<MAX_INSERTSIZE_CNV;i++){///i=insert length
        if(insertLenCount[i]!=0){///not equal to 0 means this insert length has some count
            extraLen=insertVsContigWeightedLen[i][contigNo]*(copyCount-1);
            prevEffective=weightedEffectiveLengths[i];
            newEffective=weightedEffectiveLengths[i]+extraLen;
            //cout<<"likeli before insertSize incorporation "<<likeli<<endl;
            likeli-=((insertLenCount[i]/2.0)*(log(1.0/prevEffective)-log(1.0/newEffective)));
            likelihoodEffective-=((insertLenCount[i]/2.0)*(log(1.0/prevEffective)-log(1.0/newEffective)));
            //sumInsert+=insertLenCount[i]/2.0;
            //cout<<"likeli after insertSize incorporation "<<likeli<<endl;
            //cout<<"insertSize hocche "<<insertLenCount[i]<<endl;
        }
    }
    //cout<<"total inserts "<<sumInsert<<endl;
    cout<<"likeli after insertSize incorporation "<<likelihoodEffective<<endl;
    map<int, vector<long long> >::iterator it1;
    for(it1=readsInContig.begin();it1!=readsInContig.end();it1++){///iterates over every contig and its set of reads

        if(it1->first==contigNo){///basay change(mim)
            vector< long long > :: iterator it2;
            for(it2=it1->second.begin();it2!=it1->second.end();it2++){
                int sz=multiMapReads[*it2].size();
                double sum=0.0,prevProb=0.0,newProb=0.0;
                for(int j=0;j<sz;j++){
                    sum+=multiMapReads[*it2][j].prob;
                }
                prevProb=log(sum);
                double sum1=0.0;
                for(int j=0;j<sz;j++){///this iteration is for new prob calculation
                    if(multiMapReads[*it2][j].contigNo==contigNo){
                        sum1+=multiMapReads[*it2][j].prob;
                    }
                }
                //cout<<"multi sum "<<sum1<<endl;
                newProb=log(sum+sum1*(copyCount-1));///one copy's prob already has been calculated in prevProb.So copycount-1 is used
                //cout<<"likeli sum change, old"<<likeli<<"contig "<<contigNo<<endl;
                likeli=likeli-prevProb+newProb;
                likelihoodMulti=likelihoodMulti-prevProb+newProb;
                //cout<<"likeli sum change, new"<<likeli<<"contig "<<contigNo<<endl;
                //cout<<"prevsum newsum "<<prevProb<<" "<<newProb<<endl;
            }
        }
    }
    cout<<"likeli after multimap incorporation "<<likelihoodMulti<<endl;
    return likeli;
}*/


double computeLikelihoodForCNVdependent(int copyCount, double prevLikelihoodFinal, int contigNo)
{
    int cnt=0;
    double likelihoodEffective=prevLikelihoodFinal,likelihoodMulti=prevLikelihoodFinal;
    //double sumInsert=0;
    long int n=uniqueMappedReadCount[contigNo];
    //cout<<"likeli before unique "<<prevLikelihoodFinal<<"contig "<<contigNo<<endl;
    double likeli=prevLikelihoodFinal - n*log(1) + n*log(copyCount);///uniquely mapped handled
    //cout<<"likeli after unique "<<likeli<<endl;
    double extraLen;///extra length is the amount of
                ///increase of effective length due to increase of a copy number

    double prevEffective,newEffective;
    for(int i=0;i<MAX_INSERTSIZE_CNV;i++){///i=insert length
        double extraLenforCNV=0;
        if(insertLenCount[i]!=0){///not equal to 0 means this insert length has some count

            extraLen=insertVsContigWeightedLen[i][contigNo]*(copyCount-1);
            for(int j=0;j<contigCopyNumber.size();j++){
                extraLenforCNV+=insertVsContigWeightedLen[i][j]*(contigCopyNumber[j]-1);
            }

            prevEffective=weightedEffectiveLengths[i]+extraLenforCNV;
            newEffective=weightedEffectiveLengths[i]+extraLenforCNV+extraLen;
            //cout<<"likeli before insertSize incorporation "<<likeli<<endl;
            likeli-=((insertLenCount[i]/2.0)*(log(1.0/prevEffective)-log(1.0/newEffective)));
            likelihoodEffective-=((insertLenCount[i]/2.0)*(log(1.0/prevEffective)-log(1.0/newEffective)));
            //sumInsert+=insertLenCount[i]/2.0;
            //cout<<"likeli after insertSize incorporation "<<likeli<<endl;
            //cout<<"insertSize hocche "<<insertLenCount[i]<<endl;
        }
    }
    //cout<<"total inserts "<<sumInsert<<endl;
    //cout<<"likeli after insertSize incorporation "<<likelihoodEffective<<endl;
    map<int, vector<long long> >::iterator it1;
    for(it1=readsInContig.begin();it1!=readsInContig.end();it1++){///iterates over every contig and its set of reads

        if(it1->first==contigNo){///basay change(mim)
            vector< long long > :: iterator it2;
            for(it2=it1->second.begin();it2!=it1->second.end();it2++){
                int sz=multiMapReads[*it2].size();
                double sum=0.0,prevProb=0.0,newProb=0.0;
                int copyFact=1;
                for(int j=0;j<sz;j++){
                    if(contigCopyNumber.size()>multiMapReads[*it2][j].contigNo){
                        copyFact=contigCopyNumber[multiMapReads[*it2][j].contigNo];
                    }
                    else copyFact=1;
                    sum+=(multiMapReads[*it2][j].prob)*copyFact;
                }
                prevProb=log(sum);
                double sum1=0.0;
                for(int j=0;j<sz;j++){///this iteration is for new prob calculation
                    if(multiMapReads[*it2][j].contigNo==contigNo){
                        sum1+=multiMapReads[*it2][j].prob;
                    }
                }
                //cout<<"multi sum "<<sum1<<endl;
                newProb=log(sum+sum1*(copyCount-1));///one copy's prob already has been calculated in prevProb.So copycount-1 is used
                //cout<<"likeli sum change, old"<<likeli<<"contig "<<contigNo<<endl;
                likeli=likeli-prevProb+newProb;
                likelihoodMulti=likelihoodMulti-prevProb+newProb;
                //cout<<"likeli sum change, new"<<likeli<<"contig "<<contigNo<<endl;
                //cout<<"prevsum newsum "<<prevProb<<" "<<newProb<<endl;
            }
        }
    }
    //cout<<"likeli after multimap incorporation "<<likelihoodMulti<<endl;
    cout<<"likelihood final for contig "<<contigNo<<" copy number "<<copyCount<<" is : "<<likeli<<endl;
    return likeli;
}

double computeCNVdependent()
{
    cout<<"entering Compute CNV"<<endl;
    ContigVsCopyVsProb.clear();
    contigCopyNumber.clear();

    int limit=MAX_COPY_NUMBER;

    double prevLikelihood=likelihoodOriginal, currLikelihood=likelihoodOriginal;
    double prevLikelihoodFinal=likelihoodOriginal;

    int totalContigs=contigLengths.size();
    for(int i=0; i<totalContigs; i++)
    {
        vector <int>copyVsLike;///index=copy number val= likelihood
        copyVsLike.push_back(-1);///0 copy count .so bad value
        copyVsLike.push_back(currLikelihood);///1st position of the vector.copy number 1.
                                             ///contains main likelihood
        int copyCount=1;
        prevLikelihoodFinal=prevLikelihood;
        cout<<"uniqe mapped reads "<<uniqueMappedReadCount[i]<<endl;
        while(limit>copyCount)
        {
            copyCount++;
            currLikelihood=computeLikelihoodForCNVdependent(copyCount,prevLikelihoodFinal, i);
            cout<<"Difference "<<currLikelihood-prevLikelihood<<endl;

            if(currLikelihood<prevLikelihood ||limit==copyCount)
            {
                //cout<<"contig "<<i<<"copyCount with low likelihood "<<currLikelihood<<endl;
                ContigVsCopyVsProb.push_back(copyVsLike);
                contigCopyNumber.push_back(copyCount-1);
                currLikelihood=prevLikelihood;
                //prevLikelihood=likelihoodOriginal;

                break;
            }
            else
            {
                copyVsLike.push_back(currLikelihood);
                prevLikelihood=currLikelihood;
            }
        }
    }

    return prevLikelihood;
}

/*double computeCNV()
{
    cout<<"entering Compute CNV"<<endl;
    ContigVsCopyVsProb.clear();
    contigCopyNumber.clear();

    int limit=MAX_COPY_NUMBER;

    double prevLikelihood=likelihoodOriginal, currLikelihood=likelihoodOriginal;

    int totalContigs=contigLengths.size();
    for(int i=0; i<totalContigs; i++)
    {
        vector <int>copyVsLike;///index=copy number val= likelihood
        copyVsLike.push_back(-1);///0 copy count .so bad value
        copyVsLike.push_back(currLikelihood);///1st position of the vector.copy number 1.
                                             ///contains main likelihood
        int copyCount=1;
        prevLikelihood=likelihoodOriginal;

        while(limit>copyCount)
        {
            copyCount++;
            currLikelihood=computeLikelihoodForCNV(copyCount,likelihoodOriginal, i);

            if(currLikelihood<prevLikelihood ||limit==copyCount)
            {
                //cout<<"contig "<<i<<"copyCount with low likelihood "<<currLikelihood<<endl;
                ContigVsCopyVsProb.push_back(copyVsLike);
                contigCopyNumber.push_back(copyCount-1);
                currLikelihood=likelihoodOriginal;
                //prevLikelihood=likelihoodOriginal;

                break;
            }
            else
            {
                copyVsLike.push_back(currLikelihood);
                prevLikelihood=currLikelihood;
            }
        }
    }

    return prevLikelihood;
}
*/


void calculateGcVsWeight()
{
    ifstream inFile("infoOutput1.txt");

    long int count =0, curr_pos=0,pos, curr_contig=0, prev_contig, prev_pos=-1;
    long int readLen,insertSize;///pos is the starting position of the read.
    double prob;
    map <long int, vector<long int> > mp;
    inFile >> curr_contig >> pos >> prob >> readLen >> insertSize;
    prev_contig = curr_contig;
    priority_queue<long int, vector<long int>, greater<long int> > pq;
    //joyanta
    cout<<"Entering gcVsweight calculation While 1 loop"<<endl;
    while(1){
        long int pushVal=min(pos+readLen-1,contigLengths[curr_contig]-1);
        pq.push(pushVal);

        while(curr_pos<pos){///we haven't reached the file given read starting pos
            if(pq.top()>=curr_pos){
                if(pq.empty()==false){
                    mp[curr_contig].push_back(count);
                }
            }
            else{
                while(pq.top()<curr_pos && pq.empty()==false){
                    pq.pop();
                    count--;
                }
                mp[curr_contig].push_back(count);
            }
            curr_pos++;
        }
        //joyanta
        //cout<< "reached the file given pos"<<endl;
        ///now reached the file given pos
        while(pq.top()<curr_pos && pq.empty()==false){
                pq.pop();
                count--;
        }
        count++;
        if(prev_pos==pos) ///checking if the previous read of the file was also from the same pos
        {
            mp[curr_contig].pop_back();
        }
        mp[curr_contig].push_back(count);

        if(inFile >> curr_contig >> pos >> prob >> readLen >> insertSize){
        	//joyanta
            //cout<< "In file line 1841"<<endl;
            prev_pos=curr_pos;
            if(pos>prev_pos && curr_contig==prev_contig)///eitar mane ki???
            {
                curr_pos++;
            }
            if(curr_contig!=prev_contig){

                ///filling the prev contig
                curr_pos++;
                int contigSize=contigLengths[prev_contig];
                while(curr_pos<contigSize){
                    if(pq.top()>=curr_pos){
                        if(pq.empty()==false){
                            mp[prev_contig].push_back(count);
                        }
                    }
                    else{
                        while(pq.empty()==false && pq.top()<curr_pos){
                            //cout<<pq.top()<<endl;
                            pq.pop();
                            count--;
                        }
                        mp[prev_contig].push_back(count);
                    }
                    curr_pos++;
                }

                count=0;
                curr_pos=0;
                prev_contig=curr_contig;
                prev_pos=-1;
                while(pq.empty()==false){
                    pq.pop();
                }
            }
        }
        else
        {
        	//joyanta
        	//cout<< "In filling the prev contig conditional"<<endl;
            ///filling the prev contig
                curr_pos++;
                int contigSize=contigLengths[prev_contig];
                while(curr_pos<contigSize){
                    if(pq.top()>=curr_pos){
                            if(pq.empty()==false){
                            mp[prev_contig].push_back(count);
                        }
                    }
                    else{
                        while(pq.top()<curr_pos && pq.empty()==false){
                            pq.pop();
                            count--;
                        }
                        mp[prev_contig].push_back(count);
                    }
                    curr_pos++;
                }
                break;
        }
    }

    inFile.close();
    /*for(int i=0; i<3;i++){
        for(int j=0;j<5;j++){
            cout<<mp[i][j]<<" ";
        }
        cout<<endl;
    }*/

    gcVSweight.clear();
    vector<int> counter;
    //joyanta
    cout<< "Cleared gcVSWeight" << endl;
    for(int i=0; i<100; i++)
    {
        gcVSweight.push_back(0.0);
        counter.push_back(0);
    }
    
    //joyanta
    cout<< "contigLengths.size() "<< contigLengths.size() << endl;
    int gc_row = sizeof(gcContentArr) / sizeof(gcContentArr[0]);
    int gc_column = sizeof(gcContentArr[0])/sizeof(gcContentArr[0][0]);
    cout<< "gcContentArr size: Row "<< gc_row << " Column " << gc_column << endl;
    int mp_row = sizeof(mp) / sizeof(mp[0]);
    int mp_column = sizeof(mp[0])/sizeof(mp[0][0]);
    cout<< "mp Size: Row "<< mp_row << " Column "<< mp_column << endl;
    cout<< "gcVSweight Size : "<< sizeof(gcVSweight)<<endl; 

    for(int i=0; i<contigLengths.size(); i++)
    {
        int contigSize=contigLengths[i];
        //joyanta
        cout<< "ContigSIze of "<< i << " th contig is "<< contigSize << endl;

        for(int j=0; j<contigSize; j++)
        {
            int gcValue=int(100.0*gcContentArr[i][j]); 
            gcVSweight[gcValue]+=mp[i][j];
            if(mp[i][j]!=0)
                counter[gcValue]++;
        }
        cout<< "After for loop" << endl;
    }

    for(int i=0;i<100; i++)
    {
        if(counter[i]!=0)
        {
            //gcVSweight[i]=gcVSweight[i]/counter[i];
        }
        gcVSweight[i]=1.0;
        //gcVSweight[i]+=1.0;
        //cout<<gcVSweight[i]<<endl;
    }


    return;
}

void printHelp()
{

    cout<<"cgal v0.9.5-beta"<<endl;
    cout<<"----------------"<<endl;
    cout<<endl;
    cout<<"cgal - computes likelihood"<<endl;
    cout<<"Usage:"<<endl;
    cout<<"cgal [options] <contigfile.sam>"<<endl;
    cout<<endl;
    cout<<"Required arguments:"<<endl;
    cout<<"<contigfile.sam>\t Assembly file in FASTA format"<<endl;
    cout<<endl;
    cout<<"Options:"<<endl;
    cout<<"-h [--help]\t\t Prints this message"<<endl;
    cout<<endl;
    cout<<"Output: "<<endl;
    cout<<"(In file out.txt) <numberContigs> <totalLikelihood> <mappedLikelihood> <unmappedLikelihood> <noReads> <noReadsUnmapped>"<<endl;
    cout<<"<numberContigs>\t\t Number of contigs"<<endl;
    cout<<"<totalLikelihood>\t Total log likelihood value"<<endl;
    cout<<"<mappedLikelihood>\t Likelihood value of reads mapped by the mapping tool"<<endl;
    cout<<"<unmappedLikelihood>\t Likelihood value corresponding to reads not mapped by alignment tool"<<endl;
    cout<<"<noReads>\t\t Total number of paired-end reads"<<endl;
    cout<<"<noReadsUnmapped>\t Number of reads not mapped by the alignment tool"<<endl;
    cout<<endl;
    exit(1);

}

///Change
void freeMemory(){
    delete[] errorPosDist;
	delete[] inPosDist;
	delete[] inLengthDist;
	delete[] delPosDist;
	delete[] delLengthDist;
	delete[] insertLengthDist;
	delete[] noErrorProbs;
	delete[] effectiveLengths;

}

///change done


int main(int argc, char *argv[])
{
    if(argc<2)
    	printHelp();

    if(strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-h")==0)
    	printHelp();

    contigFileName=argv[1];

    ///edited from raihan vai
    infoFile.open("info1.txt");
    ///done

    contigFile=fopen(contigFileName, "r");
    outFile=fopen("out.txt", "w");


    if (contigFile == NULL)
    {
        printf("Can't open contig file\n");
        exit(1);
    }
    char *line= new char[MAX_REC_LEN];

    char *line1= new char[MAX_REC_LEN];
    char *line2= new char[MAX_REC_LEN];

    int read;
    int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);


    long int bufferLength=1024;

    char *contig=new char[bufferLength];
    contig[0]='\0';
    char *newcontig;
    char *contigName;
    contigLength=0;


    long int tempContigLength=0;


    while(fgets(line, MAX_FILE_READ, contigFile)!=NULL)
    {
        if(line[0]==';')
        {
            continue;
        }
        else if(line[0]=='>')
        {
            contigName=new char[strlen(line)];
            strcpy(contigName,line+1);
            contigName[strlen(contigName)-1]='\0';
            contigNames.push_back(contigName);
            if(contigLength>0)
            {
                noContigs++;
                contigs.push_back(contig);
                contigLengths.push_back(contigLength);

                totalContigLength+=contigLength;
                contigLength=0;
                bufferLength=1024;
                contig=new char[bufferLength];
                contig[0]='\0';
            }
        }
        else
        {
            read=strlen(line);
            tempContigLength=contigLength;
            if(read<MAX_FILE_READ-1)
            {
                contigLength+=(read-1);
            }
            else
            {
                contigLength+=MAX_FILE_READ-1;
                read++;

            }
            if(contigLength>bufferLength)
            {
                bufferLength=max(bufferLength*2,contigLength+1);
                newcontig=new char[bufferLength];
                strcpy(newcontig,contig);
                line[read-1]='\0';
                strcat(newcontig, line);
                delete []contig;
                contig=newcontig;
            }
            else
            {
                line[read-1]='\0';
                strcpy(contig+tempContigLength, line);
            }

        }

    }
    noContigs++;
    contigs.push_back(contig);
    contigLengths.push_back(contigLength);

    totalContigLength+=contigLength;

    fclose(contigFile);

    ///changed
    for(int i=0; i<contigs.size();i++)
        uniqueMappedReadCount.push_back(0);
    ///done

    gcWriter.open("gcValue.txt");
    updateGCArr(); ///GC content updating
    /// gc content read from file and input those values in array
    ///initializeGcArrays();
    gcWriter.close();

    print_cntg_length();

    cout<<"Total line in file is :"<<fileLineCount<<endl;
    cout<<"Total contig length is : "<<tcl<<" "<<totalContigLength<<endl;


//	cout<<"after reading contig file"<<endl;

    /*
    	use bfast or some other tool to map reads and save mapping
    */

    mapFileName="myout.sam";

    mapFile=fopen(mapFileName, "r");

    summaryFile=fopen("stat.txt","r");

    if (mapFile == NULL)
    {
        printf("Can't open map file\n");
        exit(1);
    }

    int count=0;

    fscanf(summaryFile,"%ld %ld %d %d %d",&totalCount, &unCount, &toAlign, &maxReadLength, &MAX_INSERT_SIZE);


    initInsertCounts(5000);


    initErrorTypes(maxReadLength);


    itoa(maxReadLength, noErrorCigar, 10);
    strcpy(noErrorMD,"MD:Z:");
    strcat(noErrorMD,noErrorCigar);
    strcat(noErrorCigar,"M");

//	cout<<"after allocation"<<endl;

    int noMatches=0;

    while(fgets(line1, MAX_FILE_READ, mapFile)!=NULL)
    {

        if(line1[0]=='@')
            continue;

        processMapping(line1);

        count+=1;

    }

    fclose(mapFile);
    fclose(summaryFile);

    ///changes

    computeProbabilitesForGC();


    double abcd=computeLikelihoodForGC(mapFileName);


    ///edited from raihan vai
    infoFile.close();

	unixSort();

	cout<<"before gcVsweight calculation"<<endl;

	calculateGcVsWeight();


	///done


    ///changes
    cout<<"After gcVsWeight caculation"<<endl;

    computeProbabilites();
    cout<<"after prob calculation"<<endl;

    double val1=computeLikelihood(mapFileName);

    cout<<"After likelihood calculation"<<endl;

    freeMemory();

    double value=val1,val2=0.0;

//	cout<<"after val 1"<<endl;

    /*val2=computeLikelihood("unmappedOut.sam");

//	cout<<"after val 2"<<endl;


    if(unCount==0)
    {
        value=val1;
    }
    else if(unCount>toAlign)
    {
        val2=val2/toAlign*unCount;
        value=val1+val2;
    }
    else
    {
        value=val1+val2;
    }*/

    ///changed
    likelihoodOriginal = value;

    ///done

    cout<<value<<endl;

    fprintf(outFile,"%d\t%f\t%f\t%f\t%ld\t%ld\n",noContigs,value,val1,val2,totalCount,unCount);

    fclose(outFile);
    ///copy number detection
    computeCNVdependent();
    ofstream cnvOut;
    cnvOut.open("ContigVsCopyNo.txt");
    for(int i=0;i<ContigVsCopyVsProb.size();i++){

        for(int j=0;j<ContigVsCopyVsProb[i].size();j++){
            cnvOut<<"Contig Number: "<<i<<" Copy Number "<<j<<" Likelihood: "<<ContigVsCopyVsProb[i][j]<<endl;
            cout<<"Contig Number: "<<i<<" Copy Number "<<j<<" Likelihood: "<<ContigVsCopyVsProb[i][j]<<endl;
        }
        cnvOut<<endl<<endl;
    }
    cnvOut.close();


    return 0;
}
