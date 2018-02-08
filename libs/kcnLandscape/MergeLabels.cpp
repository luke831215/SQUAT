/*
== Author: Yu-Jung Chang 2017/02~
Evaluate repeat by KCN
*/
//=============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include <errno.h>

using namespace std;

#define MIN(x,y) ((x < y) ? x : y)
#define MAX(x,y) ((x > y) ? x : y)

//=============================================================================
#define LINE_BUF_SIZE 10000

class LabelCopyRecord
{
public:
	char label;
	size_t copy;
};

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

//=============================================================================
bool MergeLabels(char *class_txt, char *repeat_stock, char *out_csv)
{
/*
// 1st filter
0 N: reads with N

// after BWA
1 F: reads failed to map
2 M: multiply mapped
3 U0: uniquely mapped (0 error)
4 Us: uniquely mapped (substution error)
5 Ux: uniquely mapped (other error)
new labels for sub-U: { u,s,x } 
*/
	int NumOfClass = 7;
	size_t tagCounts[NumOfClass] = {0};
	char tagChar[NumOfClass] = {'N', 'F', 'M', 'u', 's', 'x', 'V'};

	FILE *fp1, *fp2, *fpout;
	char line_buf[LINE_BUF_SIZE];

	// Open files
	OpenFile(class_txt, "rt", fp1);
	OpenFile(repeat_stock, "rt", fp2);
	OpenFile(out_csv, "wt", fpout);
	vector<LabelCopyRecord> readtag;
	
	// 1. Load all-read_class: readID\tClassName 
	while (fgets(line_buf, LINE_BUF_SIZE, fp1) != NULL)
	{
		strtok(line_buf, "\t\n");
		string tagstr = strtok(NULL, "\t\n");
//		printf("<%s><%s>\n", line_buf, tagstr.c_str());
//		getchar();
		LabelCopyRecord rec;
		char tagno;
		if (tagstr[0] == 'U') 
		{
			rec.copy = 1;
			if (tagstr[1] == '0')
				tagno = 3;
			else if (tagstr[1] == 's')
				tagno = 4;
			else // Ux
				tagno = 5;
		}
		else // N, M, F, R
		{
			rec.copy = 0;
			if (tagstr[0] == 'N')
				tagno = 0;
			else if (tagstr[0] == 'F')
				tagno = 1;
			else if (tagstr[0] == 'M')
				tagno = 2;
			else if (tagstr[0] == 'V')
				tagno = 6;
		}
		
		tagCounts[tagno]++;
		rec.label = tagChar[tagno]; 
		readtag.push_back( rec );
	}
	fclose(fp1);

	// 2. Load repeat.stock
	size_t NumOfRepStock = 0;
	while (fgets(line_buf, LINE_BUF_SIZE, fp2) != NULL)
	{
		size_t rid = (size_t)atoi(strtok(line_buf, "\t\n"));
		size_t copy = (size_t)atoi(strtok(NULL, "\t\n"));
		readtag[rid].copy = copy;
		NumOfRepStock++;
	}
	fclose(fp2);

	printf("\n#read = %lu\n", readtag.size());
	for (char i=0; i<NumOfClass; i++)
	{
		printf("  #%c = %lu\t(%.1f%%)\n", tagChar[i], tagCounts[i], (double)tagCounts[i]*100.0/(double)readtag.size());
	}
	printf("#reads-in-RepStock = %lu\n", NumOfRepStock);

	// 3. output csv
//	fprintf(fpout, "readID,label,copy\n");
	for (size_t i=0; i<readtag.size(); i++)
	{
		fprintf(fpout, "%lu,%c,%lu\n", i, readtag[i].label, readtag[i].copy);
	}
	fclose(fpout);

	return true;
}



//=============================================================================
int main(int argc, char **argv)
{
//=============================================================================
printf("=== MergeLabels: merge all-reads_class.txt and *repeats.stock into csv ===\n\n");
	if(argc != 4)
	{
		printf("Usage: MergeLabels all-reads_class.txt repeats.stock all-reads_class2.csv\n");
		printf("=== By Yu-Jung Chang, Last modified: 2017.05 ===\n\n");

		return 1;
	}

//=============================================================================
	MergeLabels(argv[1], argv[2], argv[3]);

//=============================================================================
	return 0;
}

