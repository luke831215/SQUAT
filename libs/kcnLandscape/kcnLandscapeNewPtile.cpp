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
#define NumOfBatch 400
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

string GetFeatures(ReadRecord& rrec, CKMCFile& kc)
{
	string read = rrec.read;
	string qstr = rrec.qstr;
	string ret = "";
	char tmps[1000];
	
	// Get KCN
	vector<uint32> kcnts = {0};
	kc.GetCountersForRead(read, kcnts); // if rlen<K, kcnts.size()==0
//printf("readlen=%lu, #kmer=%lu\n", read.size(), kcnts.size());
/*
		// Compute features: HCKRatio, HCK%, MCK%, LCK%,ECK%, UCK%
		vector<char> KCNtypes;
		KCNtypes.resize(kcnts.size());

		// Initialize KCNtypes
		char types[] = {'H', 'M', 'L', 'U', 'E'};
		for (size_t i=0; i < kcnts.size(); i++) 
		{
			if (kcnts[i] >= t2)
			{
				KCNtypes[i] = 0; //'H'
			}
			else if (kcnts[i] <= t1)
			{
				if (kcnts[i] > 0)
					KCNtypes[i] = 2; //'L'
				else
					KCNtypes[i] = 3; //'U'
			}
			else
			{
				KCNtypes[i] = 1; //'M'
			}
		}
		
		// mark E
		// kcn 0 --> read 0..24
		// kcn i --> read i..i+25-1
		size_t cntTypes[5] = {0};
		int Kmer = (int)(read.size() - kcnts.size() + 1);
		for (int i=0; i<KCNtypes.size(); i++)
		{
			// update the type if necessary
			if (KCNtypes[i] == 3) // is L
			{
				for (size_t j=i; j < i+Kmer ; j++)
				{
					if(qstr[j] - QOffset <= QTH) // Is E
					{
						KCNtypes[i] == 4; // is E 
						break;
					}
				}
			}
			
			cntTypes[KCNtypes[i]]++;
		}

		double kmernum = (double)KCNtypes.size();
		double TypeCKP[5] = {0}; // type-copy-kmer-percent
		for (int i=0; i<5; i++)
		{
			TypeCKP[i] = (double)cntTypes[i] / kmernum;
		}

		double Hratio;
		if (TypeCKP[0] == 0)
			Hratio = 0;
		else // cntH > 0 
			Hratio = (double)cntTypes[0] / (double)(cntTypes[0] + cntTypes[1]);
		

		// Get footcopy, meanKCN
		vector<size_t> LMH;
		double meanKCN = 0.0;
		for (int i=0; i<kcnts.size(); i++)
		{
			meanKCN += (double)kcnts[i];
			if (KCNtypes[i] < 3) // 0,1,2 => H,M,L
				LMH.push_back(kcnts[i]);
		}
		meanKCN = meanKCN / kmernum;
		
		size_t footCopy = 0;
		if (LMH.size() > 0)
		{
			qsort(&LMH[0], LMH.size(), sizeof(size_t), cmp<size_t>);
			footCopy = LMH[0];
		}

		// Get meanQ
		size_t sumQ = 0;
		for (size_t i=0; i<qstr.size(); i++)
		{
			sumQ += (size_t)(qstr[i] - QOffset);
		}
		double meanQ = (double)sumQ / (double)qstr.size();
		
		
		// Output
		// readID, label, rcopy, HCK, LCK, UCK+ECK, Hratio, footCopy, meanQ, hiCompl
		sprintf(tmps, "%.4f,%.4f,%.4f,%.4f,%.1f,%lu,%.1f,%.4f",
			Hratio, TypeCKP[0], TypeCKP[2], TypeCKP[3]+TypeCKP[4], meanKCN, footCopy, meanQ, GetHiComplex(read));

	ret = tmps;
*/	

	// Get meanQ
	size_t sumQ = 0;
	for (size_t i=0; i<qstr.size(); i++)
	{
		sumQ += (size_t)(qstr[i] - QOffset);
	}
	double meanQ = (double)sumQ / (double)qstr.size();

	// Get meanKCN
	double meanKCN = 0.0;
	if (kcnts.size() > 0) // // if rlen<K, kcnts.size()==0
	{
		for (int i=0; i<kcnts.size(); i++)
		{
			meanKCN += (double)kcnts[i];
		}
		meanKCN = meanKCN / (double)kcnts.size();
		sprintf(tmps, "%.1f,%.1f", meanQ, meanKCN);
	}
	else
	{
		sprintf(tmps, "%.1f,-0.0", meanQ);
	}
	ret = tmps;

	// Get Q percentile
	vector<char> Qary;
	for (size_t i=0; i<qstr.size(); i++)
		Qary.push_back(qstr[i] - QOffset);
	
	QaryCountingSort(Qary);
	for (int i=0; i<=100; i+=TileStep)
	{
		sprintf(tmps, ",%.2f", GetPercentile<char>((double)i/100.0, Qary) );
		ret += tmps;
	}
	
	// Get K percentile
	if (kcnts.size() == 0) // // if rlen<K, kcnts.size()==0
	{
		for (int i=0; i<=100; i+=TileStep)
		{
			sprintf(tmps, ",%.2f", 0.0);
			ret += tmps;
		}
	}
	else
	{
		qsort(&kcnts[0], kcnts.size(), sizeof(uint32), cmp<uint32>);
		for (int i=0; i<=100; i+=TileStep)
		{
			sprintf(tmps, ",%.2f", GetPercentile<uint32>((double)i/100.0, kcnts) );
			ret += tmps;
		}
	}
	
	ret += "\n";
	
	return ret;
}

//=============================================================================
bool kcnLandscapeNewPtile(char *fq, char *db, char *class2, char * joincsv, int NumOfRec)
{
	int grpBy = 1; 
/*	if (grpBy < 0 || NumOfRec < 0)
	{
		printf("Error: NumOfRecords and GroupBy both need to be unsigned integers\n");
		return false;
	}*/
	
	printf("Input:\t%s\nKMC-DB:\t%s\n2-class CSV:\t%s\nOutput:\t%s\nNumOfRecords:\t%d\n\n", fq, db, class2, joincsv, NumOfRec);

	// Open files
	FILE *ffq, *fc2, *fls;
	char *p, buf[LINE_BUF_SIZE];

	// a. Load class2
//	printf("a. Load class2 file (%s).\n", class2); fflush(stdout);
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
//	printf("b. Open KMC file (%s).\n", db); fflush(stdout);
	CKMCFile kc;
	if (!kc.OpenForRA(db))
	{
		printf("Error: Cannot open the KMC-DB (%s)\n", db);
		return false;
	}

	// c. Load fq and open joincsv
//	printf("c. Load fastq file (%s) & the output join.csv file (%s).\n", fq, joincsv); fflush(stdout);
	OpenFile(fq, "rt", ffq);

	OpenFile(joincsv, "wt", fls);
	fprintf(fls, "readID,label,rcopy,meanQ,meanKCN");
//	fprintf(fls, "readID,label,rcopy,Hratio,HCK%%,LCK%%,UECK%%,meanKCN,footCopy,meanQ,hiCompl%%");

	for (int i=0; i<=100; i+=TileStep)
		fprintf(fls, ",Q%d%%", i);

	for (int i=0; i<=100; i+=TileStep)
		fprintf(fls, ",K%d%%", i);

	fprintf(fls, "\n");
	
	omp_set_num_threads(NumOfThreads);

	size_t ReadCount = 0;
	double TotalLen = 0.0;
	size_t line = 0; // # of lines
	size_t MinLen=0xFFFFFFFF, MaxLen=0;
	size_t idx = 0;
	vector<ReadRecord> ReadAry;
	vector<string> OutAry1, OutAry2;
	// read fastq file
//	printf("c-1. Read the fastq file & generate the join.csv file.\n"); fflush(stdout);
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

//	printf("c-1-1. [%lu] read seq=%s\n", ReadCount, rrec.read.c_str()); getchar();

		char tmps[1000];
		sprintf(tmps, "%lu,%c,%lu,", ReadCount, readClass2[ReadCount].label, readClass2[ReadCount].copy);
		OutAry1.push_back(tmps);

		ReadAry.push_back(rrec);
		if (ReadAry.size() == NumOfBatch)
		{
			OutAry2.resize(ReadAry.size());
			
			#pragma omp parallel for
			for (int i=0; i<ReadAry.size(); i++)
			{
				OutAry2[i] = GetFeatures(ReadAry[i], kc);
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
				OutAry2[i] = GetFeatures(ReadAry[i], kc);
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
	printf("=== kcnLandscapeNewPtile: Gen join read records (per 5%% percentiles for Q & K) ===\n");
	if (argc == 6)
	{
		if ( !kcnLandscapeNewPtile( argv[1], argv[2], argv[3], argv[4], atoi(argv[5]) ) )
		{
			printf("Error in kcnLandscapeNewPtile()\n");
			return -1;
		}
	}
	else
	{
		printf("Usage: kcnLandscapeNewPtile Input.fastq KMC-db all-reads_class2.csv all-reads_join.csv NumOfRecords\n");
		printf("option: Output all records: set NumOfRecords = 0\n");
		printf("=== By Yu-Jung Chang, Last modified: 2017.10 ===\n\n");

		return 1;
	}
	

//=============================================================================
	return 0;
}

