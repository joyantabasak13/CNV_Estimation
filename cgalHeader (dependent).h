#include<bits/stdc++.h>
#include <iostream>
#include<fstream>
#include <string>
#include <limits>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
using namespace std;

#define _USE_MATH_DEFINES
#include <math.h>
#include "cgal.h"

#define MAX_INSERTSIZE_CNV 1000
#define MAX_REC_LEN 1024
  ///  edited  ///
//#define readSize 101

int insertLenCount[MAX_INSERTSIZE_CNV];///index=insert length...value=count


long int infoLineCount=0;

typedef struct information
{
    double prob;
    long int pos;
    char contigField[15];
    long int readLen;
    char* conName;
    int insertSize;
} info;

typedef struct buildString
{
    long int pos;
    string cigar;
    string read;
} bs;

typedef struct readCigar
{
    string cigar;
    string read;
    int contigNo;
} CR;

static int countInfoFile=0;

 /// edited  ///

ofstream infoFile; /// edited

///edited
void writeInfoToFile(vector<info> &multiMap, map<long int, string> &, map<long int, string> &);
void createInfo(double errorProb1, long int pos1, char* contigField1, char* readString1, vector<info> &M);
void unixSort();




void createInfo(int insertSize,char *qName, double errorProb1, long int pos1, char* contigField1, char* readString1, vector<info> &M)
{
    info temp1;
    temp1.prob = errorProb1;
    temp1.pos = pos1;
    strcpy(temp1.contigField, contigField1);
    temp1.readLen = strlen(readString1);
    temp1.conName = qName;
    temp1.insertSize=insertSize;
    M.push_back(temp1);
}
void unixSort()
{
    cout << "Sorting info1" << endl;
    system("sort -n -k1,1 -k2,2 info1.txt > infoOutput1.txt");

}

void writeInfoToFile(vector<info> &multiMap, map<long int, string> &cigar, map<long int, string> &read)
{
    vector<double> cdf;
    vector<info> data(multiMap);// = getClonedMultiMap(multiMap);
    //cout<<data.size()<<"  "<<multiMap.size()<<endl;
    double s=0,s1=0;
    for(int i=0; i<multiMap.size(); i++)
    {
        info temp = multiMap[i];
        //fprintf(infoFile, "(%ld, %ld) --> %lf\n",temp.pos1,temp.pos2,temp.prob);
        s += temp.prob;


    }
    srand(time(NULL));
    for(int i=0; i<multiMap.size(); i++)
    {
        info temp = multiMap[i];
        //cout<<temp.prob<<endl;
        multiMap[i].prob = temp.prob/s;
    }
    double r = ((double) rand() / (RAND_MAX));
    ///basay change(mim)
    for(int i=0; i<multiMap.size(); i++)
    {
        info temp = multiMap[i];
        s1 += temp.prob;

        info temp1 = data[i];
        /*for(int k=0;k<multiMap.size();k++){
            info ob=data[k];
            if(ob.insertSize>=MAX_INSERTSIZE_CNV-1){
                insertLenCount[MAX_INSERTSIZE_CNV-1]++;
            }
            else{///approximation...insert length greater than 1000 are approximated to be 1000
                insertLenCount[ob.insertSize]++;///increasing the count of this insertsize
            }
        }*/
        ///basay change done(mim)
        if(r<=s1)
        {
            if(temp1.insertSize>=MAX_INSERTSIZE_CNV-1){
                insertLenCount[MAX_INSERTSIZE_CNV-1]++;
            }
            else{///approximation...insert length greater than 1000 are approximated to be 1000
                insertLenCount[temp1.insertSize]++;///increasing the count of this insertsize
            }
            infoFile << temp1.contigField<< " " << temp1.pos << " " << temp1.prob << " " << temp1.readLen << " " << temp1.insertSize << endl;
            infoLineCount++;
            int contigNo=atoi(temp1.contigField);
            //overlappingReadsInContig[contigNo];

            //cout<<iteration<<endl;
            break;
        }
    }

    countInfoFile++;
//    cout<<"Write info to file"<<endl;
}

