/*
== Author: Yu-Jung Chang 2017/06~
*/
//=============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include "./kmc_api/kmc_file.h"
#include <math.h>
#include <stdlib.h> // for qsort
#include <set>
#include <omp.h>
using namespace std;

//=============================================================================
#define LINE_BUF_SIZE 10000
#define QOffset 33
#define QTH 15
#define NumOfBatch 800
#define NumOfThreads 16
#define TileStep 5
 
class ReadRecord
{
public:
	string read, qstr;
};

class LabelCopyRecord
{
public:
	char label;
	size_t copy;
};

template <class T>
int cmp(const void * a, const void * b)
{
 	if ( *(T*)a <  *(T*)b )
		return -1;
	else if ( *(T*)a >  *(T*)b )
		return 1;
	else
 		return 0;
}

//=============================================================================
void QaryCountingSort(vector<char>& Qary)
{
	size_t Qcnts[50] = {0};
	// Get counting result Qcnts[] first
	for (size_t i=0; i < Qary.size(); ++i) // i: 0-41
		Qcnts[ Qary[i] ]++; 
	
	for (size_t i=0, k=0; k < Qary.size(); ++i) // i: 0-41
		while (Qcnts[i]--)
			Qary[k++] = (char)i;
}


/*
numpy.percentile
interpolation : {．linear・, ．lower・, ．higher・, ．midpoint・, ．nearest・}
This optional parameter specifies the interpolation method to use when the desired quantile lies between two data points i < j:
linear: i + (j - i) * fraction, where fraction is the fractional part of the index surrounded by i and j.
lower: i.
higher: j.
nearest: i or j, whichever is nearest.
midpoint: (i + j) / 2.
*/

/*
// nearest
size_t PercentTile2Rank(double percent, size_t N) // return the 0-based rank of percent tile
{
	// use The Nearest Rank method
	if (percent == 0)
		return 0;
	else
		return (size_t)ceil(percent * (double)N / 100.0) - 1;
}
*/

// linear interpolation of closest ranks
template <class T>
double GetPercentile(double percent, vector<T>& A)
{
	if (percent <= 0.0)
		return (double)A.front();
	if (percent >= 1.0)
		return (double)A.back();
		
	double portion = percent * (double)(A.size()-1);
	size_t idx1 = (size_t)floor(portion);
	size_t idx2 = (size_t)ceil(portion);
//	printf(">%lu,%lu %.2f,%.2f\n", idx1, idx2, portion, A[idx1] + (A[idx2]-A[idx1])*(portion-(double)idx1));

	if (idx1 == idx2)
		return (double)A[idx1];
	else
		return (double)A[idx1] + (double)(A[idx2]-A[idx1])*(portion-(double)idx1);
}

string GetDNAKey(string seq)
{
	string ret = "";
	
	// RevComp
	for (size_t i=seq.size()-1; i>=1; i--)
	{
		switch (seq[i])
		{
			case 'A': ret += 'T'; break;
			case 'C': ret += 'G'; break;
			case 'G': ret += 'C'; break;
			case 'T': ret += 'A'; break;
			default: ret += 'N'; break;
		}
	}
	
	if (ret > seq)
		return ret = seq;

	return ret;
}

double GetHiComplex(string seq)
{
	set<string> myset;

	for (size_t i=0; i<seq.size()-2; i++)
	{
		string Lmer = seq.substr(i, 3);
		myset.insert(GetDNAKey(Lmer));
	}
	
	return (double)myset.size() / (double)(seq.size()-2); // L-3+1=L-2
}

//=============================================================================
// File operations
void OpenFile(char* filename, const char* mode, FILE*& fp)
{
	fp = fopen(filename, mode);
	if (fp == NULL)
	{
		printf("Error: cannot open file (%s)!\n", filename);
		printf("%s\n", strerror(errno));
		exit(EXIT_FAILURE);;
	}
}

// readID,label,rcopy,Quality score,K-mer CN-76
string GetFeaturesA(ReadRecord& rrec, CKMCFile& kc, int ReadLenMax)
{
	string read = rrec.read;
	string qstr = rrec.qstr;
	string ret = "";
	
	// Get KCN
	vector<uint32> kcnts = {0};
	kc.GetCountersForRead(read, kcnts);

	char tmps[1000];
	for (size_t i=0; i<qstr.size(); i++)
	{
		sprintf(tmps, ",%d", qstr[i]-33);
		ret += tmps;
	}

	if (qstr.size() < ReadLenMax)
	{
		for (size_t i=qstr.size(); i<ReadLenMax; i++)
			ret += ",0";
	}

	for (size_t i=0; i<kcnts.size(); i++)
	{
		sprintf(tmps, ",%u", kcnts[i]);
		ret += tmps;
	}

	if (kcnts.size() < ReadLenMax-24) // ReadLenMax-25+1 = ReadLenMax-24
	{
		for (size_t i=kcnts.size(); i<ReadLenMax-24; i++)
			ret += ",0";
	}
	ret += "\n";
	
	return ret;
}

//=============================================================================
bool kcnLandscapeNewPrawA(char *fq, char *db, char *class2, char * joincsv, int NumOfRec, int ReadLenMax)
{
	int grpBy = 1; 
/*	if (grpBy < 0 || NumOfRec < 0)
	{
		printf("Error: NumOfRecords and GroupBy both need to be unsigned integers\n");
		return false;
	}*/
	
	printf("Input:\t%s\nKMC-DB:\t%s\n2-class CSV:\t%s\nOutput:\t%s\nNumOfRecords:\t%d\nReadLenMax:\t%d\n\n", fq, db, class2, joincsv, NumOfRec, ReadLenMax);

	// Open files
	FILE *ffq, *fc2, *fls;
	char *p, buf[LINE_BUF_SIZE];

	// a. Load class2
	OpenFile(class2, "rt", fc2);
	vector<LabelCopyRecord> readClass2;
	while ((p=fgets(buf, LINE_BUF_SIZE, fc2)) != NULL)
	{
		LabelCopyRecord rec;
		strtok(buf, ",\n"); // skip readID
		p = strtok(NULL, ",\n");
		rec.label = p[0];
		rec.copy = atoi( strtok(NULL, ",\n") );

		readClass2.push_back(rec);
	}
	fclose(fc2);
	printf("%lu read-records have been loaded from %s.\n", readClass2.size(), class2);

	// b. Open kcDB file for random acess
	CKMCFile kc;
	if (!kc.OpenForRA(db))
	{
		printf("Error: Cannot open the KMC-DB (%s)\n", db);
		return false;
	}

	// c. Load fq and open joincsv
	OpenFile(fq, "rt", ffq);

	OpenFile(joincsv, "wt", fls);
	fprintf(fls, "readID,label,rcopy,Q[1..%d],KCN[1..%d]\n", ReadLenMax, ReadLenMax-25+1);

	omp_set_num_threads(NumOfThreads);

	size_t ReadCount = 0;
	double TotalLen = 0.0;
	size_t line = 0; // # of lines
	size_t MinLen=0xFFFFFFFF, MaxLen=0;
	size_t idx = 0;
	vector<ReadRecord> ReadAry;
	vector<string> OutAry1, OutAry2;
	// read fastq file
	while ((p=fgets(buf, LINE_BUF_SIZE, ffq)) != NULL) // get the 1st line per 4 lines
	{
		//string read, qstr, readID;
		ReadRecord rrec;
		
		// line 1
		line++;
		if ( p == NULL || buf[0] != '@') // get the 1st line: the title. If fails
		{
			printf("FASTQ file format error at line#%lu: \"%s\"\n", line, p);
			return false;
		}
		//strtok(buf, "\n");
		//readID = &buf[1]; // skip @ 

		// line 2
		line++;
		if (fgets(buf, LINE_BUF_SIZE, ffq) == NULL) // get the 2nd line: the seq. If fails
		{
			printf("FASTQ file format error at line#%lu: \"%s\"\n", line, buf);
			return false;
		}
		strtok(buf, "\n");
		rrec.read = buf;

		// line 3
		line++;
		if (fgets(buf, LINE_BUF_SIZE, ffq) == NULL || buf[0] != '+') // get the 3rd line: the q title. If fails
		{
			printf("FASTQ file format error at line#%lu: \"%s\"\n", line, buf);
			return false;
		}

		// line 4
		line++;
		if (fgets(buf, LINE_BUF_SIZE, ffq) == NULL) // get the 4th line: the qscore. If fails
		{
			printf("FASTQ file format error at line#%lu: \"%s\"\n", line, buf);
			return false;
		}
		strtok(buf, "\n");
		rrec.qstr = buf;
		if (rrec.read.size() != rrec.qstr.size())
		{
			printf("FASTQ file format error at line#%lu: incorrect length of Q-string\n", line);
			return false;
		}

		char tmps[100];
		sprintf(tmps, "%lu,%c,%lu", ReadCount, readClass2[ReadCount].label, readClass2[ReadCount].copy);
		OutAry1.push_back(tmps);

		ReadAry.push_back(rrec);
		if (ReadAry.size() == NumOfBatch)
		{
			OutAry2.resize(ReadAry.size());
			
			#pragma omp parallel for
			for (int i=0; i<ReadAry.size(); i++)
			{
				OutAry2[i] = GetFeaturesA(ReadAry[i], kc, ReadLenMax);
			}

			for (int i=0; i<ReadAry.size(); i++)
			{
				fprintf(fls, "%s%s", OutAry1[i].c_str(), OutAry2[i].c_str());
			}

			ReadAry.clear();
			OutAry1.clear();
			OutAry2.clear();
		}
		
		// statistics
		ReadCount++;
		TotalLen += (double)rrec.read.size();
		MinLen = (rrec.read.size() < MinLen) ? rrec.read.size() : MinLen;
		MaxLen = (rrec.read.size() > MaxLen) ? rrec.read.size() : MaxLen;
		
		if (NumOfRec <= 0)
			continue;
		if (ReadCount >= NumOfRec)
			break;
	}
	fclose(ffq);

	// Compute & output the last batch
	if (ReadAry.size() > 0)
	{
//			printf("#read of LastBatch: %lu\n", ReadAry.size());
			OutAry2.resize(ReadAry.size());
			
			#pragma omp parallel for
			for (int i=0; i<ReadAry.size(); i++)
			{
				OutAry2[i] = GetFeaturesA(ReadAry[i], kc, ReadLenMax);
			}

			for (int i=0; i<ReadAry.size(); i++)
			{
				fprintf(fls, "%s%s", OutAry1[i].c_str(), OutAry2[i].c_str());
			}

			ReadAry.clear();
			OutAry1.clear();
			OutAry2.clear();
	}
	fclose(fls);
	
	printf("#read generated:\t%lu\n", ReadCount);
	printf("#base generated:\t%.0f\n", TotalLen);
	printf("Mean ReadLen:\t%.1f bp\n", TotalLen/(double)ReadCount);
	printf("Min ReadLen:\t%lu bp\n", MinLen);
	printf("Max ReadLen:\t%lu bp\n", MaxLen);
			
	return true;
}



//=============================================================================
int main(int argc, char **argv)
{
//=============================================================================
	printf("=== kcnLandscapeNewPrawA: Read fastq, its KMC-db & class2.csv to gen join read records ===\n");
	if (argc == 7)
	{
		if ( !kcnLandscapeNewPrawA( argv[1], argv[2], argv[3], argv[4], atoi(argv[5]), atoi(argv[6]) ) )
		{
			printf("Error in kcnLandscapeNew()\n");
			return -1;
		}
	}
	else
	{
		printf("Usage: kcnLandscapeNewPrawA Input.fastq KMC-db all-reads_class2.csv all-reads_join.csv NumOfRecords ReadLenMax\n");
		printf("option: Output all records: set NumOfRecords = 0\n");
		printf("=== By Yu-Jung Chang, Last modified: 2017.10 ===\n\n");

		return 1;
	}
	

//=============================================================================
	return 0;
}

