/*
== Author: Yu-Jung Chang; update: 2017/10
Give %HiQ(x) and read PE FASTQ files to select PEs 
*/
//=============================================================================
#include <stdio.h>
#include <string.h>
#include <string>
#include <math.h>
#include <stdlib.h> // for atoi

using namespace std;

#define MIN(x,y) ((x < y) ? x : y)
#define MAX(x,y) ((x > y) ? x : y)

//=============================================================================
#define LINE_BUF_SIZE 10000
#define QSCORE_SIZE 42
#define QCharOffset 33
#define SEQ_SIZE 400
#define HiQCellSize 200
#define DefalutHiQTh 20

//=============================================================================
// Select PEs By LowQPercent
bool peSelect_HiQ(char *r1, char *r2, char *outPrjName, double HiQPercentTh, char HiQTh)
{
	FILE *fp1, *fp2, *fpo1, *fpo2, *fpcsv;
	char line_buf[LINE_BUF_SIZE], line_buf2[LINE_BUF_SIZE];
	double TotalLen = 0.0, sTotalLen = 0.0; // Total length of all output sequences
	size_t PECount = 0, sPECount = 0; // # of PE sequences
	size_t line = 0; // # of fq lines
	char tmps[512];
	double Qsum = 0.0, sQsum = 0.0;
	size_t MinSeqLen, sMinSeqLen;
	size_t MaxSeqLen, sMaxSeqLen;
	MinSeqLen = sMinSeqLen = LINE_BUF_SIZE;
	MaxSeqLen = sMaxSeqLen = 0;
	char outFiles[LINE_BUF_SIZE];

	// Open files
	fp1 = fopen(r1, "rt");
	if (fp1 == NULL)
	{
		printf("Read FQ File 1 Error!\n");
		return false;
	}

	fp2 = fopen(r2, "rt");
	if (fp2 == NULL)
	{
		printf("Read FQ File 2 Error!\n");
		return false;
	}

	sprintf(tmps, "%s.csv", outPrjName);
	fpcsv = fopen(tmps, "wt");
	if (fpcsv == NULL)
	{
		printf("Open csv OUT File (%s) Error!\n", tmps);
		return false;
	}
	sprintf(outFiles, "%s", tmps);

	sprintf(tmps, "%s-r1.fq", outPrjName);
	fpo1 = fopen(tmps, "wt");
	if (fpo1 == NULL)
	{
		printf("Open r1 OUT File (%s) Error!\n", tmps);
		return false;
	}
	sprintf(outFiles, "%s,%s", outFiles, tmps);

	sprintf(tmps, "%s-r2.fq", outPrjName);
	fpo2 = fopen(tmps, "wt");
	if (fpo2 == NULL)
	{
		printf("Open r2 OUT File (%s) Error!\n", tmps);
		return false;
	}
	sprintf(outFiles, "%s,%s", outFiles, tmps);

	// Read the file and select
	string fqrec1[4], fqrec2[4]; // line 1-4
	while (fgets(line_buf, LINE_BUF_SIZE, fp1) != NULL && fgets(line_buf2, LINE_BUF_SIZE, fp2) != NULL) // get the 1st line per 4 lines
	{
		// line 1
		line++;
		if (line_buf[0] != '@' || line_buf2[0] != '@') // Is not the 1st title line
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}
		fqrec1[0] = line_buf;
		fqrec2[0] = line_buf2;

		// line 2
		line++;
		if (fgets(line_buf, LINE_BUF_SIZE, fp1) == NULL || fgets(line_buf2, LINE_BUF_SIZE, fp2) == NULL) // get the 2nd line: the seq. If fails
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}
		fqrec1[1] = line_buf;
		fqrec2[1] = line_buf2;
		
		strtok(line_buf, "\n");
		strtok(line_buf2, "\n");
		size_t SeqLen1 = strlen(line_buf);
		size_t SeqLen2 = strlen(line_buf2);

		size_t tmpi = MIN(SeqLen1, SeqLen2);
		MinSeqLen = MIN(MinSeqLen, tmpi);
		tmpi = MAX(SeqLen1, SeqLen2);
		MaxSeqLen = MAX(MaxSeqLen, tmpi);

		// line 3
		line++;
		if (fgets(line_buf, LINE_BUF_SIZE, fp1) == NULL || line_buf[0] != '+' || fgets(line_buf2, LINE_BUF_SIZE, fp2) == NULL || line_buf2[0] != '+') // get the 3rd line: the q title. If fails
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}
		fqrec1[2] = line_buf;
		fqrec2[2] = line_buf2;

		// line 4
		line++;
		if (fgets(line_buf, LINE_BUF_SIZE, fp1) == NULL || fgets(line_buf2, LINE_BUF_SIZE, fp2) == NULL) // get the 4th line: the qscore. If fails
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}
		fqrec1[3] = line_buf;
		fqrec2[3] = line_buf2;

		strtok(line_buf, "\n");
		strtok(line_buf2, "\n");
		size_t QLen1 = strlen(line_buf);
		size_t QLen2 = strlen(line_buf2);
		if (QLen1 != SeqLen1 || QLen2 != SeqLen2)
		{
			printf("FASTQ file format error at line#%lu: incorrect length of Q-string\n", line);
			return false;
		}

		double tmpPEQsum = 0.0;
		int HiQCnt1 = 0, HiQCnt2 = 0;
		for (size_t i=0; i<QLen1; i++)
		{
			char tmpq = line_buf[i]-QCharOffset;
			if (tmpq >= HiQTh)
				HiQCnt1 ++;
			tmpPEQsum += (double)tmpq;
		}
		for (size_t i=0; i<QLen2; i++)
		{
			char tmpq = line_buf2[i]-QCharOffset;
			if (tmpq >= HiQTh)
				HiQCnt2 ++;
			tmpPEQsum += (double)tmpq;
		}
		Qsum += tmpPEQsum;
		PECount++;
		TotalLen += (double)(SeqLen1+SeqLen2);
		double hiqp = (int)floor(HiQCellSize*MIN((double)HiQCnt1/(double)QLen1, (double)HiQCnt2/(double)QLen2)) * 100.0 / HiQCellSize;
//		printf("lowqp=%f\n", lowqp); getchar();
		if (hiqp >= HiQPercentTh)
		{
			// output 
			for (int i=0; i < 4; i++)
			{
				fprintf(fpo1, "%s", fqrec1[i].c_str());
				fprintf(fpo2, "%s", fqrec2[i].c_str());
			}
			sPECount++;
			sTotalLen += (double)(SeqLen1+SeqLen2);
			sQsum += tmpPEQsum;

			tmpi = MIN(SeqLen1, SeqLen2);
			sMinSeqLen = MIN(sMinSeqLen, tmpi);
			tmpi = MAX(SeqLen1, SeqLen2);
			sMaxSeqLen = MAX(sMaxSeqLen, tmpi);
		}

	}
	fclose(fp1);
	fclose(fp2);
	fclose(fpo1);
	fclose(fpo2);

	// Start outputing
	sprintf(tmps, "--- Summary ---\n");
	fputs(tmps, stdout);
	fputs(tmps, fpcsv);

	sprintf(tmps, "PE's '%%HighQ(%d) >= %.1f\n", HiQTh, HiQPercentTh);
	fputs(tmps, stdout);
	fputs(tmps, fpcsv);

	sprintf(tmps, "Name,Original,Subset,Sub/Ori\n");
	fputs(tmps, stdout);
	fputs(tmps, fpcsv);

	sprintf(tmps, "#PE,%lu,%lu,%.1f%%\n", PECount, sPECount, 100.0*(double)sPECount/(double)PECount);
	fputs(tmps, stdout);
	fputs(tmps, fpcsv);

	sprintf(tmps, "#Base,%.0f,%.0f,%.1f%%\n", TotalLen, sTotalLen, 100.0*(double)sTotalLen/(double)TotalLen);
	fputs(tmps, stdout);
	fputs(tmps, fpcsv);

	double tmpf, tmpf2;
	tmpf = TotalLen/(double)(PECount*2);
	tmpf2 = sTotalLen/(double)(sPECount*2);
	sprintf(tmps, "AvgReadLen,%.1f,%.1f,%.1f%%\n", tmpf, tmpf2, 100.0*tmpf2/tmpf);
	fputs(tmps, stdout);
	fputs(tmps, fpcsv);

	sprintf(tmps, "MaxReadLen,%lu,%lu,%.1f%%\n", MaxSeqLen, sMaxSeqLen, 100.0*(double)sMaxSeqLen/(double)MaxSeqLen);
	fputs(tmps, stdout);
	fputs(tmps, fpcsv);

	sprintf(tmps, "MinReadLen,%lu,%lu,%.1f%%\n", MinSeqLen, sMinSeqLen, 100.0*(double)sMinSeqLen/(double)MinSeqLen);
	fputs(tmps, stdout);
	fputs(tmps, fpcsv);

	tmpf = Qsum/TotalLen;
	tmpf2 = sQsum/sTotalLen;
	sprintf(tmps, "BaseQmean,%.1f,%.1f,%.1f%%\n", tmpf, tmpf2, 100.0*tmpf2/tmpf);
	fputs(tmps, stdout);
	fputs(tmps, fpcsv);

	printf("OutFile: %s\n", outFiles);
	return true;
}



//=============================================================================
int main(int argc, char **argv)
{
//=============================================================================
	if(argc != 5 && argc != 6)
	{
		printf("=== peQsubset: Select the PE subset from the input PE FASTQ files by %%HighQ(QTh) ===\n\n");
		printf("Usage: peQsubset r1.fq r2.fq outPrjName %%HighQ [QTh]\n");
		printf("Input: r1.fq (read1 fastq of PE), r2.fq (read2 fastq of PE)\n");
		printf(" -QTh: Optional. The Q-value lowerbound (Range: 0-41; default 20)\n");
		printf(" -%%HighQ: The lowerbound percent of read1's/read2's' bases with Q-values >= QTh (Range: 100.0-0.0)\n");
		printf("Output: outPrjName-r1.fq, outPrjName-r2.fq, outPrjName.csv\n");
		printf("Verson: 0.9 (2017/10) \n");
		printf("Author: Yu-Jung Chang\n\n");

		return 1;
	}

//=============================================================================
	if (argc == 5)
		peSelect_HiQ(argv[1], argv[2], argv[3], atof(argv[4]), DefalutHiQTh);
	else if (argc == 6)
		peSelect_HiQ(argv[1], argv[2], argv[3], atof(argv[4]), (char)atoi(argv[5]));

//=============================================================================
	return 0;
}

