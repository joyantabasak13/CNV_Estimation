#include <iostream>
#include <string>
#include <limits>
#include <vector>
using namespace std;

#include "cgal.h"
#include <math.h>

#define MAX_REC_LEN 1024
#define MAX_READLENGTH 500
#define MAX_NAMELENGTH 500

struct SAM
{
	char qname[MAX_NAMELENGTH];
	int flag;
	char rname[MAX_NAMELENGTH];
	int pos;
	int mapq;
	char cigar[MAX_READLENGTH];
	char rnext[MAX_NAMELENGTH];
	int pnext;
	int tlen;
	char seq[MAX_NAMELENGTH];
	char qual[MAX_NAMELENGTH];
	char md[MAX_NAMELENGTH];
	int ih;
};


FILE *outFile;
FILE *mapFile;
FILE *unFile;
FILE *statFile;

vector <SAM *> reads1;
vector <SAM *> reads2;

char line1[MAX_REC_LEN];
char line2[MAX_REC_LEN];

long int unCount=0;
long int totalCount=0;
int maxReadLength=0;
double insertSizeMean=0;
double insertSizeVar=0;
double squaredError=0;

long int MAX_FRAGMENT_SIZE=5000;

void reverse(char *reverse, char *read)
{

	char ch='A';	
	int readLength=strlen(read);
	for(int i=1; i<=readLength;i++)
	{
		ch=read[i-1];
		if(ch=='A')
			reverse[readLength-i]='T';
		else if(ch=='C')
			reverse[readLength-i]='G';
		else if(ch=='G')
			reverse[readLength-i]='C';
		else if(ch=='T')
			reverse[readLength-i]='A';
		else		
			reverse[readLength-i]='N';						
	}
	reverse[readLength]='\0';

}




void writeSam(SAM* read, FILE *out)
{
/*	fprintf(out,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->rname,read->pos,
		read->mapq,read->cigar,read->rnext,read->pnext,read->tlen,read->seq,read->qual,read->md,read->ih);
*/

	fprintf(out,"%s\t%d\t%d\t%s\t%d\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->pos,
		read->cigar,read->tlen,read->seq,read->md,read->ih);

}


void printVectors(FILE * out, FILE * u)
{
		
	SAM *read1, *read2;
	int readLength1,readLength2;
	int pos1, pos2;
	int insertSize;

	char temp[MAX_READLENGTH];

	int ih=0;

	for(int i=0;i<reads1.size();i++)
	{
		read1=reads1[i];
			
		for(int j=0;j<reads2.size();j++)
		{
			read2=reads2[j];

			if(strcmp(read1->rname,"*")==0 || strcmp(read2->rname,"*")==0 || strcmp(read1->rname,read2->rname)!=0)
			{

			}
			else
			{
				ih++;
			}

		}	
	}

	if(ih==0)
	{

		int i=0;
		int j=0;
		int ncount1=0;
		int ncount2=0;
		
		readLength1=strlen(read1->seq);
		if(readLength1>maxReadLength)
			maxReadLength=readLength1;
	
		readLength2=strlen(read2->seq);
		if(readLength2>maxReadLength)
			maxReadLength=readLength2;
	

		for(int i=0;i<readLength1;i++)
		{
			if(read1->seq[i]=='N'||read1->seq[i]=='n')
			{
				ncount1++;
			}
		}

		for(int i=0;i<readLength2;i++)
		{

			if(read2->seq[i]=='N'||read2->seq[i]=='n')
			{
				ncount2++;
			}
		}


		
		
		if(ncount1/(double)readLength1 < 0.8 && ncount2/(double)readLength2 < 0.8)

		{
			unCount++;
			totalCount++;

			


			fputs("@",u);
			fputs(read1->qname,u);
			fputs("\n",u);


			int strandNo=(read1->flag&16)>>4;
			if(strandNo==1)
			{
				reverse(temp,read1->seq);
				fputs(temp,u);
			}
			else
			{				
				fputs(read1->seq,u);
				
			}
			fputs("\n",u);

			fputs("+",u);
			fputs(read1->qname,u);
			fputs("\n",u);

			fputs(read1->qual,u);
			fputs("\n",u);

			fputs("@",u);
			fputs(read2->qname,u);
			fputs("\n",u);
		
			strandNo=(read2->flag&16)>>4;
			if(strandNo==1)
			{
				reverse(temp,read2->seq);
				fputs(temp,u);
			}
			else
			{				
				fputs(read2->seq,u);
				
			}
			fputs("\n",u);

			fputs("+",u);
			fputs(read2->qname,u);
			fputs("\n",u);

			fputs(read2->qual,u);
			fputs("\n",u);

		}

		for(int k=0;k<reads1.size();k++)
			delete reads1[k];

		for(int k=0;k<reads2.size();k++)
			delete reads2[k];
	

		reads1.clear();
		reads2.clear();

		return;
	


	}

	int i=0;
	int j=0;
	int ncount1=0;
	int ncount2=0;
	
	readLength1=strlen(read1->seq);
	if(readLength1>maxReadLength)
		maxReadLength=readLength1;
	
	readLength2=strlen(read2->seq);
	if(readLength2>maxReadLength)
		maxReadLength=readLength2;
	

	for(int i=0;i<readLength1;i++)
	{
		if(read1->seq[i]=='N'||read1->seq[i]=='n')
		{
			ncount1++;
		}
	}

	for(int i=0;i<readLength2;i++)
	{
		if(read2->seq[i]=='N'||read2->seq[i]=='n')
		{
			ncount2++;
		}
	}


		
		
	if(ncount1/(double)readLength1 >= 0.8 || ncount2/(double)readLength2 >= 0.8)
	{
		for(int k=0;k<reads1.size();k++)
			delete reads1[k];

		for(int k=0;k<reads2.size();k++)
			delete reads2[k];
	

		reads1.clear();
		reads2.clear();

		return;
	}



	totalCount++;

	for(int i=0;i<reads1.size();i++)
	{
		read1=reads1[i];
		pos1=read1->pos;
		readLength1=strlen(read1->seq);
		if(readLength1>maxReadLength)
			maxReadLength=readLength1;
	
		for(int j=0;j<reads2.size();j++)
		{
			read2=reads2[j];
			pos2=read2->pos;
			readLength2=strlen(read2->seq);
			if(readLength2>maxReadLength)
				maxReadLength=readLength2;

			if(strcmp(read1->rname,"*")==0 || strcmp(read2->rname,"*")==0 || strcmp(read1->rname,read2->rname)!=0 || abs(read2->pos-read1->pos)>MAX_FRAGMENT_SIZE)
			{
			}
			else
			{
				if(pos2>=pos1)
					insertSize=pos2-pos1+readLength2;
				else
					insertSize=pos2-pos1-readLength1;

				read1->tlen=insertSize;
				read2->tlen=-insertSize;

				writeSam(read1,out);

		//second one

				writeSam(read2,out);

				
				
			}

		}	
	}
	if(insertSize<0)
		insertSize=-insertSize;

	if(insertSize<MAX_FRAGMENT_SIZE)
	{
		insertSizeMean=(insertSizeMean*(totalCount-unCount-1)+insertSize)/(totalCount-unCount);
		squaredError+=(insertSize-insertSizeMean)*(insertSize-insertSizeMean);
	}
	
	for(int i=0;i<reads1.size();i++)
		delete reads1[i];

	for(int i=0;i<reads2.size();i++)
		delete reads2[i];
	

	reads1.clear();
	reads2.clear();

}


SAM *getSAM(char *line)
{
	SAM *sam=new SAM;
	char *temp;

	temp=strtok(line,"\t\n ");
	strcpy(sam->qname,temp);

	temp=strtok(NULL,"\t\n ");
	sam->flag=atoi(temp);

	temp=strtok(NULL,"\t\n ");
	strcpy(sam->rname,temp);

	temp=strtok(NULL,"\t\n ");
	sam->pos=atoi(temp);

	temp=strtok(NULL,"\t\n ");
	sam->mapq=atoi(temp);

	temp=strtok(NULL,"\t\n ");
	strcpy(sam->cigar,temp);
	
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->rnext,temp);

	temp=strtok(NULL,"\t\n ");
	sam->pnext=atoi(temp);

	temp=strtok(NULL,"\t\n ");
	sam->tlen=atoi(temp);

	temp=strtok(NULL,"\t\n ");
	strcpy(sam->seq,temp);
	
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->qual,temp);

	while((temp=strtok(NULL,"\t\n "))!=NULL)
	{
		if(temp[0]=='I' && temp[1]=='H')
		{
			sam->ih=atoi(temp+5);
		}
		else if(temp[0]=='M' && temp[1]=='D')
		{
			strcpy(sam->md,temp);
		}
	}


	return sam;
}

void printHelp()
{

	cout<<"cgal v0.9.5-beta"<<endl;
	cout<<"----------------"<<endl;
	cout<<endl;
	cout<<"bfastconvert - converts the map file in sam format outputted by BFAST into an internal format"<<endl;
	cout<<"Usage:"<<endl;
	cout<<"bfastconvert [options] <mapfile.sam> [maxFragmentLength]"<<endl; 
	cout<<endl;
	cout<<"Required arguments:"<<endl;
	cout<<"<mapfile.sam>\t\t Map file outputted by BFAST in sam format"<<endl;
	cout<<endl;
	cout<<"Optional arguments:"<<endl;
	cout<<"[maxFragmentLength]\t Maximum insert size (fragment length). Default 5000"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
	cout<<"-h [--help]\t\t Prints this message"<<endl;
	cout<<endl;
	exit(1);

	exit(1);

}

int main(int argc, char *argv[])
{
	/*	input contig file name, read file name
		contig file - fasta format
		read file - fastq format		
	*/


	if(argc<2)
		printHelp();

	if(strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-h")==0)
		printHelp();


	if(argc==3)
		MAX_FRAGMENT_SIZE=atoi(argv[2]);
	else
		MAX_FRAGMENT_SIZE=5000;


	char *line= new char[MAX_REC_LEN];
	char *templine= new char[MAX_REC_LEN];
	
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
	
	int nh=1;

	char * mapFileName=argv[1];


	mapFile=fopen(mapFileName, "r");

	if (mapFile == NULL) 
	{
		printf("Can't open map file\n");
		exit(1);
	}


	outFile=fopen("myout.sam","w");

	unFile=fopen("unmapped.txt","w");

	statFile=fopen("stat.txt","w");

	char *temp,nhstring[500];



	
	int it=0; 
		
	while(fgets(line, MAX_FILE_READ, mapFile)!=NULL)
	{
		
		it++;

		
		if(line[0]=='@')
			continue;

		SAM *read1=getSAM(line);

		reads1.push_back(read1);

//		cout<<read1.ih<<endl;

		nh=read1->ih;
		
		for(int i=0;i<nh-1;i++)
		{
			fgets(line, MAX_FILE_READ, mapFile);
			reads1.push_back(getSAM(line));
			
		}


		fgets(line, MAX_FILE_READ, mapFile);

		SAM *read2=getSAM(line);

		reads2.push_back(read2);

		nh=read2->ih;
		
		for(int i=0;i<nh-1;i++)
		{
			fgets(line, MAX_FILE_READ, mapFile);
			reads2.push_back(getSAM(line));
		}



	
		printVectors(outFile, unFile);
		
		if(totalCount>=10550000)
			break;

	}

	
	insertSizeVar=squaredError/(totalCount-unCount);

	int max=(int)(insertSizeMean+5*sqrt(insertSizeVar));

	fprintf(statFile,"%ld %ld %d %d",totalCount, unCount, maxReadLength, max);


	fclose(mapFile);
	fclose(outFile);
	fclose(unFile);
	fclose(statFile);

	return 0;
}
