/*
== Author: Yu-Jung Chang; update: 2018/03
Read PE-FASTQ files and generate quality distribution of the data
*/
//=============================================================================
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <string>
//#include <algorithm> // for heap
//#include <stdlib.h> // for qsort

using namespace std;

#define MIN(x,y) ((x < y) ? x : y)
#define MAX(x,y) ((x > y) ? x : y)

//=============================================================================
#define LINE_BUF_SIZE 10000
#define ALPHABET_SIZE 256
#define QSCORE_SIZE 42
#define QCharOffset 33
#define SEQ_SIZE 400
#define HiQCellSize 200

#define HiQParamSize 2
int HiQTh[HiQParamSize] = {15, 20};

//=============================================================================
int cmpchar(const void *arg1, const void *arg2)
{
  return  (int)(*(char *)arg1 - *(char *)arg2);
}

void counting_sort(char Qstr[], int Qlen, int Qcnt[])
{
	// Get counting result Qcnt[] first
	
	for (int i=0, k=0; k < Qlen; ++i) // i: 0-41
		while (Qcnt[i]--)
			Qstr[k++] = (char)i;
}

//=============================================================================
bool ProbeFASTQPE(char *r1, char *r2, char *outPrjName)
{
	FILE *fp1, *fp2, *fpout, *fphtm;
	char line_buf[LINE_BUF_SIZE], line_buf2[LINE_BUF_SIZE];
	char OutCSV[1024], OutHTM[1024];
	size_t AlphabetCount[ALPHABET_SIZE] = {0}; // counting alphabet occurrence
	double TotalLen = 0.0; // Total length of all sequences
	size_t PECount = 0; // # of PE sequences
	size_t line = 0; // # of lines
	size_t MinSeqLen = LINE_BUF_SIZE;
	size_t MaxSeqLen = 0;

	size_t QCount[QSCORE_SIZE] = {0}; // counting Q scores
	size_t MinQCount[QSCORE_SIZE] = {0}; // counting MinQ scores

	// for HiQ%
	size_t HiQPercentCount[HiQParamSize][HiQCellSize+1] = {0};
	
	// for GC%
	size_t CntGCRead[101] = {0}; // Counts of GC% for read

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

	sprintf(OutCSV, "%s.csv", outPrjName);
	fpout = fopen(OutCSV, "wt");
	if (fpout == NULL)
	{
		printf("Open OUT CSV File Error!\n");
		return false;
	}

	sprintf(OutHTM, "%s.htm", outPrjName);
	fphtm = fopen(OutHTM, "wt");
	if (fphtm == NULL)
	{
		printf("Open OUT HTM File Error!\n");
		return false;
	}

	// Read the file and check
	while (fgets(line_buf, LINE_BUF_SIZE, fp1) != NULL && fgets(line_buf2, LINE_BUF_SIZE, fp2) != NULL) // get the 1st line per 4 lines
	{
		// line 1
		line++;
		if (line_buf[0] != '@' || line_buf2[0] != '@') // Is not the 1st title line
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}

		// line 2
		line++;
		if (fgets(line_buf, LINE_BUF_SIZE, fp1) == NULL || fgets(line_buf2, LINE_BUF_SIZE, fp2) == NULL) // get the 2nd line: the seq. If fails
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}
		
		// read 1 & 2
		strtok(line_buf, "\n");
		strtok(line_buf2, "\n");
		size_t SeqLen1 = strlen(line_buf);
		size_t SeqLen2 = strlen(line_buf2);
		TotalLen += (double)(SeqLen1+SeqLen2);
		PECount++;

		size_t tmpi = MIN(SeqLen1, SeqLen2);
		MinSeqLen = MIN(MinSeqLen, tmpi);
		tmpi = MAX(SeqLen1, SeqLen2);
		MaxSeqLen = MAX(MaxSeqLen, tmpi);

		// GC% of read1
		char GCvalue;
		size_t GCcnt=0, ATcnt=0;
		for (size_t i=0; i<SeqLen1; i++)
		{
			AlphabetCount[line_buf[i]]++;
			if (line_buf[i] == 'G' || line_buf[i] == 'C')
				GCcnt++;
			else if (line_buf[i] == 'A' || line_buf[i] == 'T')
				ATcnt++;
		}
		GCvalue = (char)round(100.0*(double)GCcnt/(double)(GCcnt+ATcnt));
		CntGCRead[GCvalue]++;
//		printf("[%lu] NumGC=%lu, NumAT=%lu, GCvalue=%d\n", PECount, GCcnt[0], ATcnt[0], GCvalue); getchar();

		// GC% of read2
		GCcnt=ATcnt=0;
		for (size_t i=0; i<SeqLen2; i++)
		{
			AlphabetCount[line_buf2[i]]++;
			if (line_buf2[i] == 'G' || line_buf2[i] == 'C')
				GCcnt++;
			else if (line_buf2[i] == 'A' || line_buf2[i] == 'T')
				ATcnt++;
		}
		GCvalue = (char)round(100.0*(double)GCcnt/(double)(GCcnt+ATcnt));
		CntGCRead[GCvalue]++;

		// line 3
		line++;
		if (fgets(line_buf, LINE_BUF_SIZE, fp1) == NULL || line_buf[0] != '+' || fgets(line_buf2, LINE_BUF_SIZE, fp2) == NULL || line_buf2[0] != '+') // get the 3rd line: the q title. If fails
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}

		// line 4
		line++;
		if (fgets(line_buf, LINE_BUF_SIZE, fp1) == NULL || fgets(line_buf2, LINE_BUF_SIZE, fp2) == NULL) // get the 4th line: the qscore. If fails
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}
		strtok(line_buf, "\n");
		strtok(line_buf2, "\n");
		size_t QLen1 = strlen(line_buf);
		size_t QLen2 = strlen(line_buf2);
		if (QLen1 != SeqLen1 || QLen2 != SeqLen2)
		{
			printf("FASTQ file format error at line#%lu: incorrect length of Q-string\n", line);
			return false;
		}

		char minq1 = QSCORE_SIZE-1, minq2 = QSCORE_SIZE-1;
		int HiQCnt[HiQParamSize][2] = {0};

		for (size_t i=0; i<QLen1; i++)
		{
			char tmpq = line_buf[i]-QCharOffset;
			QCount[tmpq]++;
			minq1 = MIN(minq1, tmpq);
			for (int j=0; j < HiQParamSize; j++)
				if (tmpq >= HiQTh[j])
					HiQCnt[j][0] ++;
		}
		for (size_t i=0; i<QLen2; i++)
		{
			char tmpq = line_buf2[i]-QCharOffset;
			QCount[tmpq]++;
			minq2 = MIN(minq2, tmpq);
			for (int j=0; j < HiQParamSize; j++)
				if (tmpq >= HiQTh[j])
					HiQCnt[j][1] ++;
		}

		// for MinQ
		MinQCount[MIN(minq1, minq2)]++;

		// for HiQ%
		for (int j=0; j < HiQParamSize; j++)
		{
			// For each PE, decide its HiQ% tile number
			// LowQ by ceil; HiQ by floor
			int HiqpTileNo = (int)floor(HiQCellSize*MIN((double)HiQCnt[j][0]/(double)QLen1, (double)HiQCnt[j][1]/(double)QLen2));
			HiQPercentCount[j][HiqpTileNo]++;

		}

		// Progress
		if (PECount % 5000000 == 0) // print . per 5M PE
		{
			printf(".");
			fflush(stdout);
		}

		// for Test
//		if (PECount == 20000000)
//			break;;
	}
	fclose(fp1);
	fclose(fp2);
	printf("done\n");

// --- summary
	// to screen
	fprintf(stdout, "--- Summary of PE FASTQ ---\n");
	fprintf(stdout, "InputFile: %s,%s\n", r1, r2);
//	fprintf(stdout, "#Line,%lu\n", line);
	fprintf(stdout, "#PE: %lu\n", PECount);
	fprintf(stdout, "#Read: %lu\n", PECount*2);
	fprintf(stdout, "#Base: %.0f\n", TotalLen);
	fprintf(stdout, "AvgReadLen: %.2f\n", TotalLen/(double)(PECount*2));
	fprintf(stdout, "MinReadLen: %lu\n", MinSeqLen);
	fprintf(stdout, "MaxReadLen: %lu\n", MaxSeqLen);
	fprintf(stdout, "OutFile: %s,%s\n", OutCSV, OutHTM);

	// Start outputing to files
	fprintf(fpout, "--- Summary of PE FASTQ ---\n");
	fprintf(fpout, "FileName,%s,%s\n", r1, r2);
	fprintf(fpout, "#PE,%lu\n", PECount);
	fprintf(fpout, "#Read,%lu\n", PECount*2);
	fprintf(fpout, "#Base,%.0f\n", TotalLen);
	fprintf(fpout, "AvgReadLen,%.2f\n", TotalLen/(double)(PECount*2));
	fprintf(fpout, "MaxReadLen,%lu\n", MaxSeqLen);
	fprintf(fpout, "MinReadLen,%lu\n", MinSeqLen);

	// Output html
	fprintf(fphtm, "<html>\n<head>\n<script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script></head>\n");
	fprintf(fphtm, "<body>\n");
	fprintf(fphtm, "--- Summary of PE FASTQ ---\n");
	fprintf(fphtm, "<li>InputFile: %s,%s\n", r1, r2);
	fprintf(fphtm, "<li>#PE: %lu\n", PECount);
	fprintf(fphtm, "<li>#Read: %lu\n", PECount*2);
	fprintf(fphtm, "<li>#Base: %.0f\n", TotalLen);
	fprintf(fphtm, "<li>AvgReadLen: %.2f\n", TotalLen/(double)(PECount*2));
	fprintf(fphtm, "<li>MinReadLen: %lu\n", MinSeqLen);
	fprintf(fphtm, "<li>MaxReadLen: %lu\n", MaxSeqLen);

// --- Alphabet freq
	fprintf(fpout, "\n--- Alphabet Occurrence Count/Frequency ---\n");
	fprintf(fpout, "Name,Count,Freq%%\n");
	for (size_t i=0; i<ALPHABET_SIZE; i++)
	{
		if (AlphabetCount[i] > 0)
		{
			if (i == 'N')
				continue;
			
			double tmpf = 100.0 * (double)AlphabetCount[i] / TotalLen;
			fprintf(fpout, """%c"",%lu,%.2f%%\n", (char)i, AlphabetCount[i], tmpf);
		}
	}
	if (AlphabetCount['N'] > 0)
	{
		double tmpf = 100.0 * (double)AlphabetCount['N'] / TotalLen;
		fprintf(fpout, """%c"",%lu,%.2f%%\n", 'N', AlphabetCount['N'], tmpf);
	}

	// htm
	fprintf(fphtm, "\n<br><br>--- Alphabet Occurrence Count/Frequency ---\n");
	fprintf(fphtm, "<table border=1><tr><td>Name</td><td>Count</td><td>Freq%%</td></tr>\n");
	for (size_t i=0; i<ALPHABET_SIZE; i++)
	{
		if (AlphabetCount[i] > 0)
		{
			if (i == 'N')
				continue;
			
			double tmpf = 100.0 * (double)AlphabetCount[i] / TotalLen;
			fprintf(fphtm, "<tr><td>""%c""</td><td>%lu</td><td>%.2f%%</td></tr>\n", (char)i, AlphabetCount[i], tmpf);
		}
	}
	if (AlphabetCount['N'] > 0)
	{
		double tmpf = 100.0 * (double)AlphabetCount['N'] / TotalLen;
		fprintf(fphtm, "<tr><td>""%c""</td><td>%lu</td><td>%.2f%%</td></tr>\n", 'N', AlphabetCount['N'], tmpf);
	}
	fprintf(fphtm, "</table>\n\n");

// --- GC% dist of reads
	fprintf(fpout, "\n--- GC%% of Reads ---\n");
	fprintf(fpout, "GC%%,Count,Freq%%\n");
	for (char i=0; i<=100; i++)
	{
		if (CntGCRead[i] == 0)
			continue;
		
		double tmpf = 100.0 * (double)CntGCRead[i] / (double)(PECount*2);
		fprintf(fpout, "%d%%,%lu,%.2f%%\n", i, CntGCRead[i], tmpf);
	}

	// htm
	fprintf(fphtm, "<br><table border=0 width=900px>");
	fprintf(fphtm, "  <tr><td><br>--- GC%% distribution of Reads ---</td></tr>\n");
	fprintf(fphtm, "  <tr><td id=gc style=\"height: 300px\"></td></tr>\n");

	fprintf(fphtm, "  <tr><td><br>--- Base Q-value distribution ---</td></tr>\n");
	fprintf(fphtm, "  <tr><td id=bq style=\"height: 300px\"></td></tr>\n");

	fprintf(fphtm, "  <tr><td><br>--- MinimalQ distribution ---</td></tr>\n");
	fprintf(fphtm, "  <tr><td id=mq style=\"height: 300px\"></td></tr>\n");

	fprintf(fphtm, "  <tr><td><br>--- %%HighQ(15) distribution ---</td></tr>\n");
	fprintf(fphtm, "  <tr><td id=hq15 style=\"height: 300px\"></td></tr>\n");

	fprintf(fphtm, "  <tr><td><br>--- %%HighQ(18) distribution ---</td></tr>\n");
	fprintf(fphtm, "  <tr><td id=hq18 style=\"height: 300px\"></td></tr>\n");

	fprintf(fphtm, "  <tr><td><br>--- %%HighQ(20) distribution ---</td></tr>\n");
	fprintf(fphtm, "  <tr><td id=hq20 style=\"height: 300px\"></td></tr>\n");

	fprintf(fphtm, "  <tr><td><br>--- %%HighQ(25) distribution ---</td></tr>\n");
	fprintf(fphtm, "  <tr><td id=hq25 style=\"height: 300px\"></td></tr>\n");
	fprintf(fphtm, "</table>\n\n");

	// script
	fprintf(fphtm, "<script type=\"text/javascript\">\n");
	fprintf(fphtm, "function DrawDist() {\n");

	// option
	fprintf(fphtm, "  var opt1 = {\n");
	fprintf(fphtm, "    title: 'GC%% distribution of reads', hAxis: { title: 'GC%%' }, vAxis: { title: 'Freq%%' }, colors: ['#a52714', '#097138']\n");
	fprintf(fphtm, "  };\n");

	fprintf(fphtm, "  var d1 = new google.visualization.DataTable();\n");
	fprintf(fphtm, "  d1.addColumn('number', 'GC%%');\n");
	fprintf(fphtm, "  d1.addColumn('number', 'Freq');\n");
	fprintf(fphtm, "  d1.addRows( [ "); 

	for (char i=0; i<=100; i++)
	{
		double tmpf = 100.0 * (double)CntGCRead[i] / (double)(PECount*2);
		fprintf(fphtm, "[%d,%.2f],", i, tmpf);
		if (i % 10 == 9)
			fprintf(fphtm, "\n");
	}

	fprintf(fphtm, "  ] );\n");
	fprintf(fphtm, "  var chart1 = new google.visualization.LineChart(document.getElementById('gc'));\n");
	fprintf(fphtm, "  chart1.draw(d1, opt1);\n");
	fprintf(fphtm, "\n");


// --- Base Q dist
	fprintf(fpout, "\n--- Base Q-value Count ---\n");
	fprintf(fpout, "QValue,Count,Freq%%\n");
	for (size_t i=0; i<QSCORE_SIZE; i++)
	{
		if (QCount[i] > 0)
		{
			double tmpf = 100.0 * (double)QCount[i] / TotalLen;
			fprintf(fpout, "%lu,%lu,%.2f%%\n", i, QCount[i], tmpf);
		}
	}

	// option
	fprintf(fphtm, "  var opt2 = {\n");
	fprintf(fphtm, "    title: 'Base Q-value distribution', hAxis: { title: 'Q-value' }, vAxis: { title: 'Freq%%', format: 'percent' }, colors: ['#a52714', '#097138']\n");
	fprintf(fphtm, "  };\n");

	fprintf(fphtm, "  var d2 = new google.visualization.DataTable();\n");
	fprintf(fphtm, "  d2.addColumn('number', 'Q-value');\n");
	fprintf(fphtm, "  d2.addColumn('number', 'Freq');\n");
	fprintf(fphtm, "  d2.addRows( [ "); 

	for (size_t i=0; i<QSCORE_SIZE; i++)
	{
		double tmpf = 100.0 * (double)QCount[i] / TotalLen;
		fprintf(fphtm, "[%lu,%.2f], ", i, tmpf);
		if (i % 10 == 9)
			fprintf(fphtm, "\n");
	}

	fprintf(fphtm, "  ] );\n");
	fprintf(fphtm, "  var chart2 = new google.visualization.LineChart(document.getElementById('bq'));\n");
	fprintf(fphtm, "  chart2.draw(d2, opt2);\n");
	fprintf(fphtm, "\n");

// --- MinQ dist
	fprintf(fpout, "\n--- PE MinQ-value Count ---\n");
	fprintf(fpout, "MinQ,Count,Freq%%,cumuFreq%%\n");
	for (size_t i=0, cumuCnt=PECount; i<QSCORE_SIZE; i++)
	{
		if (MinQCount[i] > 0)
		{
			double tmpf = 100.0*(double)MinQCount[i] / (double)PECount;
			double tmpf2 = 100.0*(double)cumuCnt / (double)PECount;
			fprintf(fpout, "%lu,%lu,%.2f%%,%.2f%%\n", i, MinQCount[i], tmpf, tmpf2);
			cumuCnt -= MinQCount[i];
		}
	}

	// option
	fprintf(fphtm, "  var opt3 = {\n");
	fprintf(fphtm, "    title: 'MinQ distribution', hAxis: { title: 'MinQ-value (in reverse dir.)', direction: -1, viewWindow: { max: 41 } }, vAxis: { title: '%%', format: 'percent' }, colors: ['#a52714', '#097138']\n");
	fprintf(fphtm, "  };\n");

	fprintf(fphtm, "  var d3 = new google.visualization.DataTable();\n");
	fprintf(fphtm, "  d3.addColumn('number', 'MinQ');\n");
	fprintf(fphtm, "  d3.addColumn('number', 'Freq');\n");
	fprintf(fphtm, "  d3.addColumn('number', 'SubsetSize');\n");
	fprintf(fphtm, "  d3.addRows( [ "); 

	for (size_t i=0, cumuCnt=PECount; i<QSCORE_SIZE; i++)
	{
		double tmpf = (double)MinQCount[i] / (double)PECount;
		double tmpf2 = (double)cumuCnt / (double)PECount;
		fprintf(fphtm, "[%lu,%.4f,%.4f], ", i, tmpf, tmpf2);
		if (i % 10 == 9)
			fprintf(fphtm, "\n");
		cumuCnt -= MinQCount[i];
	}

	fprintf(fphtm, "  ] );\n");
	fprintf(fphtm, "  var chart3 = new google.visualization.LineChart(document.getElementById('mq'));\n");
	fprintf(fphtm, "  chart3.draw(d3, opt3);\n");
	fprintf(fphtm, "\n");

// --- Dist of %HighQ(x) 
	for (int j=0; j < HiQParamSize; j++)
	{
		// option
		fprintf(fphtm, "  var optH%d = {\n", HiQTh[j]);
		fprintf(fphtm, "    title: '%%HighQ(%d) distribution', hAxis: { title: '%%HighQ(%d) (in reverse dir.)', direction: -1, viewWindow: { max: 100, min: 50 } }, vAxis: { title: '%%', format: 'percent' }, colors: ['#a52714', '#097138']\n", HiQTh[j], HiQTh[j]);
		fprintf(fphtm, "  };\n");
	
		fprintf(fphtm, "  var dH%d = new google.visualization.DataTable();\n", HiQTh[j]);
		fprintf(fphtm, "  dH%d.addColumn('number', '%%HighQ');\n", HiQTh[j]);
		fprintf(fphtm, "  dH%d.addColumn('number', 'Freq');\n", HiQTh[j]);
		fprintf(fphtm, "  dH%d.addColumn('number', 'SubsetSize');\n", HiQTh[j]);
		fprintf(fphtm, "  dH%d.addRows( [ ", HiQTh[j]); 

		fprintf(fpout, "\n--- PE %%HighQ(%d)-value Count ---\n", HiQTh[j]);
		fprintf(fpout, "%%HighQ(%d),Count,Freq%%,CumuFreq%%\n", HiQTh[j]);
		for (long i=HiQCellSize, cumuCnt=0; i >= 0; i--)
		{
			if (HiQPercentCount[j][i] == 0)
				continue;

			cumuCnt += HiQPercentCount[j][i];
			double tmpf = (double)HiQPercentCount[j][i] / (double)PECount;
			double tmpf2 = (double)cumuCnt / (double)PECount;
			fprintf(fpout, "%.1f,%lu,%.2f%%,%.2f%%\n", (double)i*(100.0/HiQCellSize), HiQPercentCount[j][i], 100.0*tmpf, 100.0*tmpf2);

			fprintf(fphtm, "[%.1f,%.4f,%.4f], ", (double)i*(100.0/HiQCellSize), tmpf, tmpf2);
			if (i % 10 == 9)
				fprintf(fphtm, "\n");
		}
	
		fprintf(fphtm, "  [0,0,1] ] );\n");
		fprintf(fphtm, "  var chartH%d = new google.visualization.LineChart(document.getElementById('hq%d'));\n", HiQTh[j], HiQTh[j]);
		fprintf(fphtm, "  chartH%d.draw(dH%d, optH%d);\n", HiQTh[j], HiQTh[j], HiQTh[j]);
		fprintf(fphtm, "\n");
	}
	fprintf(fpout, "\n");


	// end of htm
	fprintf(fphtm, "}\n");
	fprintf(fphtm, "google.charts.load('current', {packages: ['corechart', 'line']});\n");
	fprintf(fphtm, "google.charts.setOnLoadCallback(DrawDist);"); 
	fprintf(fphtm, "</script>\n");
	fprintf(fphtm, "</body>\n</html>\n");

	fclose(fpout);
	fclose(fphtm);

	return true;
}



//=============================================================================
int main(int argc, char **argv)
{
//=============================================================================
	if(argc != 4)
	{
		printf("=== peQdist: Read PE-FASTQ files and generate quality distribution and GC%% of paired-end FASTQ files ===\n\n");
		printf("Usage: peQdist r1.fq r2.fq outPrjName\n");
		printf("Input: r1.fq (read1 fastq of PE), r2.fq (read2 fastq of PE) outPrjName(name of output)\n");
		printf("Output: outPrjName.htm, outPrjName.csv\n");
		printf("Verson: 0.91 (2017/11) \n");
		printf("Author: Yu-Jung Chang\n\n");

		return 1;
	}

//=============================================================================
	ProbeFASTQPE(argv[1], argv[2], argv[3]);

//=============================================================================
	return 0;
}

