/*
== Author: Yu-Jung Chang; update: 2018/03
Read single FASTQ file and generate quality distribution of the data
*/
//=============================================================================
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <string>
#include <time.h>
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

string GetCurrentTime(void)
{
  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  string ret = asctime (timeinfo);

  return ret;
}

string AddCommas(size_t value)
{
	char buf[1000];
	sprintf(buf, "%lu", value);
	string ret = buf;
	for (int i=(int)ret.size()-3; i > 0; i-=3)
	{
		// The new contents are inserted before the character at position i
		ret.insert(i, ",");
	}
	
	return ret;
}

string GetFileName(const string& s)
{
	 char sep = '/';
	#ifdef _WIN32
	   sep = '\\';
	#endif
	size_t i = s.rfind(sep, s.length());
	if (i != string::npos) {
		string str = s.substr(i+1, s.length() - i);
		return(str);
	}
	
	return("");
}

//=============================================================================
bool ProbeFASTQ(char *r1, char *outPrjName)
{
	FILE *fp1, *fphtm;
//	FILE *fpout;
	char line_buf[LINE_BUF_SIZE];
	char OutHTM[1024];
//	char OutCSV[1024];
	size_t AlphabetCount[ALPHABET_SIZE] = {0}; // counting alphabet occurrence
	size_t TotalLen = 0; // Total length of all sequences
	size_t ReadCount = 0; // # of read sequences
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
/*
	sprintf(OutCSV, "%s.csv", outPrjName);
	fpout = fopen(OutCSV, "wt");
	if (fpout == NULL)
	{
		printf("Open OUT CSV File Error!\n");
		return false;
	}
*/
	sprintf(OutHTM, "%s.htm", outPrjName);
	fphtm = fopen(OutHTM, "wt");
	if (fphtm == NULL)
	{
		printf("Open OUT HTM File Error!\n");
		return false;
	}

	// Read the file and check
	while (fgets(line_buf, LINE_BUF_SIZE, fp1) != NULL) // get the 1st line per 4 lines
	{
		// line 1
		line++;
		if (line_buf[0] != '@') // Is not the 1st title line
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}

		// line 2
		line++;
		if (fgets(line_buf, LINE_BUF_SIZE, fp1) == NULL) // get the 2nd line: the seq. If fails
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}
		
		// read 1 & 2
		strtok(line_buf, "\n");
		size_t SeqLen1 = strlen(line_buf);
		TotalLen += SeqLen1;
		ReadCount++;

		MinSeqLen = MIN(MinSeqLen, SeqLen1);
		MaxSeqLen = MAX(MaxSeqLen, SeqLen1);

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

		// line 3
		line++;
		if (fgets(line_buf, LINE_BUF_SIZE, fp1) == NULL || line_buf[0] != '+') // get the 3rd line: the q title. If fails
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}

		// line 4
		line++;
		if (fgets(line_buf, LINE_BUF_SIZE, fp1) == NULL) // get the 4th line: the qscore. If fails
		{
			printf("FASTQ file format error at line#%lu\n", line);
			return false;
		}
		strtok(line_buf, "\n");
		size_t QLen1 = strlen(line_buf);
		if (QLen1 != SeqLen1)
		{
			printf("FASTQ file format error at line#%lu: incorrect length of Q-string\n", line);
			return false;
		}

		char minq1 = QSCORE_SIZE-1;
		int HiQCnt[HiQParamSize] = {0};

		for (size_t i=0; i<QLen1; i++)
		{
			char tmpq = line_buf[i]-QCharOffset;
			QCount[tmpq]++;
			minq1 = MIN(minq1, tmpq);
			for (int j=0; j < HiQParamSize; j++)
				if (tmpq >= HiQTh[j])
					HiQCnt[j] ++;
		}

		// for MinQ
		MinQCount[minq1]++;

		// for HiQ%
		for (int j=0; j < HiQParamSize; j++)
		{
			// For each PE, decide its HiQ% tile number
			// LowQ by ceil; HiQ by floor
			int HiqpTileNo = (int)floor(HiQCellSize*(double)HiQCnt[j]/(double)QLen1);
			HiQPercentCount[j][HiqpTileNo]++;

		}

		// Progress
		if (ReadCount % 10000000 == 0) // print . per 5M PE
		{
			printf(".");
			fflush(stdout);
		}

	}
	fclose(fp1);
	printf("done\n");

// --- summary
	// to screen
	fprintf(stdout, "--- Summary of FASTQ ---\n");
	fprintf(stdout, "InputFile: %s\n", GetFileName(r1).c_str());
	fprintf(stdout, "#Read: %lu\n", ReadCount);
	fprintf(stdout, "#Base: %lu\n", TotalLen);
	fprintf(stdout, "AvgReadLen: %.2f\n", (double)TotalLen/(double)ReadCount);
	fprintf(stdout, "MinReadLen: %lu\n", MinSeqLen);
	fprintf(stdout, "MaxReadLen: %lu\n", MaxSeqLen);
//	fprintf(stdout, "OutFile: %s,%s\n", OutCSV, OutHTM);
	fprintf(stdout, "OutFile: %s\n", OutHTM);
/*
	// Start outputing to files
	fprintf(fpout, "--- Summary of FASTQ ---\n");
	fprintf(fpout, "FileName,%s\n", r1);
	fprintf(fpout, "#Read,%lu\n", ReadCount);
	fprintf(fpout, "#Base,%lu\n", TotalLen);
	fprintf(fpout, "AvgReadLen,%.2f\n", (double)TotalLen/(double)ReadCount);
	fprintf(fpout, "MaxReadLen,%lu\n", MaxSeqLen);
	fprintf(fpout, "MinReadLen,%lu\n", MinSeqLen);
*/
	// Output html
	fprintf(fphtm, "<html>\n<head>\n<title>Pre-assembly SQUAT report</title><script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script>\n");
	fprintf(fphtm, "<link rel=\"stylesheet\" type=\"text/css\" href=\"link/template.css\">\n<script src=\"link/template.js\"></script>\n");
	fprintf(fphtm, "</head>\n<body>\n");
	
	fprintf(fphtm, "<div class=\"header\">\n  <div id=\"header_title\">Pre-Assembly SQUAT Report</div>\n  <div id=\"header_filename\">%s</div>\n</div>", GetCurrentTime().c_str());

	fprintf(fphtm, "<div class=\"summary\">\n  <h2 style=\"text-align: center\">Summary</h2>\n");
	
	fprintf(fphtm, "<div class=\"ac\"\n>");
	fprintf(fphtm, "<input class=\"ac-input\" id=\"ac-1\" name=\"ac-1\" type=\"checkbox\" /\n>");
	fprintf(fphtm, "<label class=\"ac-label\" for=\"ac-1\">Basic Statistics</label>");
	fprintf(fphtm, "<article class=\"ac-text\">");
	fprintf(fphtm, "<div class=\"ac-sub\">");
	fprintf(fphtm, "<span class=\"ac-row\" onclick=\"link('#Fs')\">");
	fprintf(fphtm, "Attributes of FASTQ");
	fprintf(fphtm, "</span>");
	fprintf(fphtm, "</div>");
	fprintf(fphtm, "<div class=\"ac-sub\">");
	fprintf(fphtm, "<span class=\"ac-row\" onclick=\"link('#Fa')\">");
	fprintf(fphtm, "Alphabet Frequency & GC content");
	fprintf(fphtm, "</span>");
	fprintf(fphtm, "</div>");
	fprintf(fphtm, "</article>");
	fprintf(fphtm, "</div>");

	fprintf(fphtm, "<div class=\"ac\"\n>");
	fprintf(fphtm, "<input class=\"ac-input\" id=\"ac-2\" name=\"ac-2\" type=\"checkbox\" /\n>");
	fprintf(fphtm, "<label class=\"ac-label\" for=\"ac-2\">Quality Statistics</label>");
	fprintf(fphtm, "<article class=\"ac-text\">");
	fprintf(fphtm, "<div class=\"ac-sub\">");
	fprintf(fphtm, "<span class=\"ac-row\" onclick=\"link('#Fb')\">");
	fprintf(fphtm, "Distribution of Bases' Quality Values");
	fprintf(fphtm, "</span>");
	fprintf(fphtm, "</div>");
	fprintf(fphtm, "<div class=\"ac-sub\">");
	fprintf(fphtm, "<span class=\"ac-row\" onclick=\"link('#Fm')\">");
	fprintf(fphtm, "Distribution of Reads' MinimaQ Values");
	fprintf(fphtm, "</span>");
	fprintf(fphtm, "</div>");
	fprintf(fphtm, "<div class=\"ac-sub\">");
	fprintf(fphtm, "<span class=\"ac-row\" onclick=\"link('#Fd')\">");
	fprintf(fphtm, "Covergae of Reads with Sufficient High-Quality Bases");
	fprintf(fphtm, "</span>");
	fprintf(fphtm, "</div>");
	fprintf(fphtm, "</article>");
	fprintf(fphtm, "</div>");
	fprintf(fphtm, "</div>");
		
		
		  
		    

	/*<button class=\"accordion\">Basic Statistics</button>\n");
	fprintf(fphtm, "    <ul><li><a href='#Fs'>Attributes of FASTQ</a></li>\n");
	fprintf(fphtm, "    <li><a href='#Fa'>Alphabet Frequency & GC content</a></li></ul>\n");
	fprintf(fphtm, "  <button class=\"accordion\">Quality Statistics</a></button>\n");
	fprintf(fphtm, "    <ul><li><a href='#Fb'>Distribution of Bases' Quality Values</a></li>\n");
	fprintf(fphtm, "    <li><a href='#Fm'>Distribution of Reads' MinimaQ Values</a></li>\n");
	fprintf(fphtm, "    <li><a href='#Fd'>Covergae of Reads with Sufficient High-Quality Bases</a></li></ul>\n");
	fprintf(fphtm, "</div>\n<div class=\"main\">\n"); */

	fprintf(fphtm, "<div class=\"main\" id=\"main\" onscroll=scrollFunction()>");
	fprintf(fphtm, "<button onclick=\"topFunction()\" id=\"btpBtn\" title=\"Go to top\"><i class=\"up\"></i>Top</button>");
	
	fprintf(fphtm, "  <div id=Fs>\n");
	fprintf(fphtm, "  <h3 style='color: darkblue;'>Attributes of FASTQ</h3>\n");
	fprintf(fphtm, "  <table border=1 width=60%% style='margin-left:120px;'><tr align=center><th width=50%%>Name</th><th>Value</th></tr>\n");
	fprintf(fphtm, "    <tr align=center><td>InputFile</td><td>%s</td></tr>\n", GetFileName(r1).c_str());
	fprintf(fphtm, "    <tr align=center><td>#Read</td><td>%s</td></tr>\n", AddCommas(ReadCount).c_str());
	fprintf(fphtm, "    <tr align=center><td>#Base</td><td>%s</td></tr>\n", AddCommas(TotalLen).c_str());
	fprintf(fphtm, "    <tr align=center><td>AvgReadLen</td><td>%.2f</td></tr>\n", TotalLen/(double)ReadCount);
	fprintf(fphtm, "    <tr align=center><td>MinReadLen</td><td>%lu</td></tr>\n", MinSeqLen);
	fprintf(fphtm, "    <tr align=center><td>MaxReadLen</td><td>%lu</td></tr>\n", MaxSeqLen);
	fprintf(fphtm, "  </table>\n  </div><br><br><br>\n");

// --- Alphabet freq
/*
	fprintf(fpout, "\n--- Alphabet Occurrence Count/Frequency ---\n");
	fprintf(fpout, "Name,Count,Freq%%\n");
	for (size_t i=0; i<ALPHABET_SIZE; i++)
	{
		if (AlphabetCount[i] > 0)
		{
			if (i == 'N')
				continue;
			
			double tmpf = 100.0 * (double)AlphabetCount[i] / TotalLen;
			fprintf(fpout, """%c"",%s,%.2f%%\n", (char)i, AddCommas(AlphabetCount[i]).c_str(), tmpf);
		}
	}
	if (AlphabetCount['N'] > 0)
	{
		double tmpf = 100.0 * (double)AlphabetCount['N'] / TotalLen;
		fprintf(fpout, """%c"",%lu,%.2f%%\n", 'N', AlphabetCount['N'], tmpf);
	}
*/

	// htm
	fprintf(fphtm, "  <div id=Fa>\n");
	fprintf(fphtm, "  <h3 style='color: darkblue;'>Alphabet Frequency & GC content</h3>\n");
	fprintf(fphtm, "  <table border=1 width=60%% style='margin-left:120px;'><tr align=center><th>Name</th><th>Count</th><th>Freq%%</th></tr>\n");
	for (size_t i=0; i<ALPHABET_SIZE; i++)
	{
		if (AlphabetCount[i] > 0)
		{
			if (i == 'N')
				continue;
			
			double tmpf = 100.0 * (double)AlphabetCount[i] / (double)TotalLen;
			fprintf(fphtm, "    <tr align=center><td>""%c""</td><td>%s</td><td>%.2f%%</td></tr>\n", (char)i, AddCommas(AlphabetCount[i]).c_str(), tmpf);
		}
	}
	if (AlphabetCount['N'] > 0)
	{
		double tmpf = 100.0 * (double)AlphabetCount['N'] / (double)TotalLen;
		fprintf(fphtm, "    <tr align=center><td>""%c""</td><td>%s</td><td>%.2f%%</td></tr>\n", 'N', AddCommas(AlphabetCount['N']).c_str(), tmpf);
	}

	// GC%
	{
		double tmpf = 100.0 * (double)(AlphabetCount['C']+AlphabetCount['G']) / (double)TotalLen;
		fprintf(fphtm, "    <tr align=center><td>GC%%</td><td>-</td><td>%.2f%%</td></tr>\n", tmpf);
	}
	fprintf(fphtm, "  </table>\n");

// --- GC% dist of reads
/*	fprintf(fpout, "\n--- GC%% of Reads ---\n");
	fprintf(fpout, "GC%%,Count,Freq%%\n");
	for (char i=0; i<=100; i++)
	{
		if (CntGCRead[i] == 0)
			continue;
		
		double tmpf = 100.0 * (double)CntGCRead[i] / (double)ReadCount;
		fprintf(fpout, "%d%%,%lu,%.2f%%\n", i, CntGCRead[i], tmpf);
	}
*/
	// htm gc plot
	fprintf(fphtm, "  <table border=0 width=80%%>");
	fprintf(fphtm, "    <tr><td id=gc style='height: 300px'></td></tr>\n");
	fprintf(fphtm, "  </table>\n");
	fprintf(fphtm, "  </div><br><br><br>\n");

	// BaseQ
	fprintf(fphtm, "  <div id=Fb>\n");
	fprintf(fphtm, "  <h3 style='color: darkblue;'>Distribution of Bases' Quality Values</h3>\n");

	size_t BQsum[4] = {0}; // <15, 15-19, 20-29, 30+
	for (size_t i=0; i<QSCORE_SIZE; i++)
	{
		if (i < 15)
			BQsum[0] += QCount[i];
		else if (i >= 15 && i < 20)
			BQsum[1] += QCount[i];
		else if (i >= 20 && i < 30)
			BQsum[2] += QCount[i];
		else
			BQsum[3] += QCount[i];
	}
	fprintf(fphtm, "  <table border=1 width=60%% style='margin-left:120px;'><tr align=center><th width=50%%>Name</th><th>AreaFreq</th></tr>\n");
	fprintf(fphtm, "    <tr align=center><td>Q30 & above</td><td>%.1f%%</td></tr>\n", 100.0*(double)BQsum[3]/(double)TotalLen);
	fprintf(fphtm, "    <tr align=center><td>Q20-Q29</td><td>%.1f%%</td></tr>\n", 100.0*(double)BQsum[2]/(double)TotalLen);
	fprintf(fphtm, "    <tr align=center><td>Q15-Q19</td><td>%.1f%%</td></tr>\n", 100.0*(double)BQsum[1]/(double)TotalLen);
	fprintf(fphtm, "    <tr align=center><td>< Q15</td><td>%.1f%%</td></tr>\n", 100.0*(double)BQsum[0]/(double)TotalLen);
	fprintf(fphtm, "  </table>\n\n");

	fprintf(fphtm, "  <table border=0 width=80%%>");
	fprintf(fphtm, "  <tr><td id=bq style=\"height: 300px\"></td></tr>\n");
	fprintf(fphtm, "  </table>\n  </div><br><br><br>\n");

	// MinQ
	fprintf(fphtm, "  <div id=Fm>\n");
	fprintf(fphtm, "  <h3 style='color: darkblue;'>Distribution of Reads' MinimalQ Values</h3>\n");

	size_t MQsum[3] = {0}; // >= 10, >=15, >=20
	for (size_t i=0; i<QSCORE_SIZE; i++)
	{
		if (i >= 10)
		{
			MQsum[0] += MinQCount[i];
			if (i >= 15)
			{
				MQsum[1] += MinQCount[i];
				if (i >= 20)
				{
					MQsum[2] += MinQCount[i];
				}
			}
		}
	}
	fprintf(fphtm, "  <table border=1 width=60%% style='margin-left:120px;'><tr align=center><th width=50%%>Name</th><th>AreaFreq</th></tr>\n");
	fprintf(fphtm, "    <tr align=center><td>%% of reads whose bases are all Q20 & above</td><td>%.1f%%</td></tr>\n", 100.0*(double)MQsum[2]/(double)ReadCount);
	fprintf(fphtm, "    <tr align=center><td>%% of reads whose bases are all Q15 & above</td><td>%.1f%%</td></tr>\n", 100.0*(double)MQsum[1]/(double)ReadCount);
	fprintf(fphtm, "    <tr align=center><td>%% of reads whose bases are all Q10 & above</td><td>%.1f%%</td></tr>\n", 100.0*(double)MQsum[0]/(double)ReadCount);
	fprintf(fphtm, "  </table>\n\n");

	fprintf(fphtm, "  <table border=0 width=80%%>");
	fprintf(fphtm, "  <tr><td id=mq style=\"height: 300px\"></td></tr>\n");
	fprintf(fphtm, "  </table>\n  </div><br><br><br>\n");

	fprintf(fphtm, "  <div id=Fd>\n");
	fprintf(fphtm, "  <h3 style='color: darkblue;'>Coverage of Reads with Sufficient High-Quality Bases</h3>\n");

	double HQcov[4] = {0}; // 95-Q20, 95-Q15, 90-Q20, 90-Q15
	for (long i=HiQCellSize, cumuCnt[HiQParamSize]={0}; i >= 0; i--)
	{
		double cumuRatio[HiQParamSize];
		double cov = (double)i*100.0 / (double)HiQCellSize;
		for (int k=0; k<HiQParamSize; k++)
		{
			cumuCnt[k] += HiQPercentCount[k][i];
//			printf("cov=%.1f, k=%d, cumuCnt[k]=%lu\n", cov, k, cumuCnt[k]);
//			getchar();
			if (cov == 95.0)
			{
				if (k == 1) // Q20
					HQcov[0] = (double)cumuCnt[k] / (double)ReadCount;
				else if (k == 0) // Q15
					HQcov[1] = (double)cumuCnt[k] / (double)ReadCount;
			}
			else if (cov == 90.0)
			{
				if (k == 1) // Q20
					HQcov[2] = (double)cumuCnt[k] / (double)ReadCount;
				else if (k == 0) // Q15
					HQcov[3] = (double)cumuCnt[k] / (double)ReadCount;
			}
		}
	}
	fprintf(fphtm, "  <table border=1 width=60%% style='margin-left:120px;'><tr align=center><th width=50%%>Name</th><th>Freq</th></tr>\n");
	fprintf(fphtm, "    <tr align=center><td>Coverage of reads that >=95%% of bases with Q20 & above</td><td>%.1f%%</td></tr>\n", 100.0*HQcov[0]);
	fprintf(fphtm, "    <tr align=center><td>Coverage of reads that >=95%% of bases with Q15 & above</td><td>%.1f%%</td></tr>\n", 100.0*HQcov[1]);
	fprintf(fphtm, "    <tr align=center><td>Coverage of reads that >=90%% of bases with Q20 & above</td><td>%.1f%%</td></tr>\n", 100.0*HQcov[2]);
	fprintf(fphtm, "    <tr align=center><td>Coverage of reads that >=90%% of bases with Q15 & above</td><td>%.1f%%</td></tr>\n", 100.0*HQcov[3]);
	fprintf(fphtm, "  </table>\n\n");

	fprintf(fphtm, "  <table border=0 width=80%%>");
	fprintf(fphtm, "    <tr><td align=center><img src=link/HighQ.gif alt='%%HighQ(q)' style='width: 400'></img></td></tr>\n");
	fprintf(fphtm, "    <tr><td id=hq style=\"height: 300px\"></td></tr>\n");
	fprintf(fphtm, "  </table>\n  </div><br><br><br>\n");
	fprintf(fphtm, "</div>");
	
	// script
	fprintf(fphtm, "<script type=\"text/javascript\">\n");
	fprintf(fphtm, "function DrawDist() {\n");

	// option
	fprintf(fphtm, "  var opt1 = {\n");
	fprintf(fphtm, "    title: 'GC%% of reads', hAxis: { title: 'GC%%' }, vAxis: { title: 'Freq%%', format: 'percent' }, colors: ['#a52714', '#097138']\n");
	fprintf(fphtm, "  };\n");

	fprintf(fphtm, "  var d1 = new google.visualization.DataTable();\n");
	fprintf(fphtm, "  d1.addColumn('number', 'GC%%');\n");
	fprintf(fphtm, "  d1.addColumn('number', 'Freq');\n");
	fprintf(fphtm, "  d1.addRows( [ "); 

	for (char i=0; i<=100; i++)
	{
		double tmpf = (double)CntGCRead[i] / (double)ReadCount;
		fprintf(fphtm, "[%d,%.4f],", i, tmpf);
		if (i % 10 == 9)
			fprintf(fphtm, "\n");
	}

	fprintf(fphtm, "  ] );\n");
	fprintf(fphtm, "  var chart1 = new google.visualization.LineChart(document.getElementById('gc'));\n");
	fprintf(fphtm, "  chart1.draw(d1, opt1);\n");
	fprintf(fphtm, "\n");

// --- Base Q dist
/*
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
*/
	// option
	fprintf(fphtm, "  var opt2 = {\n");
//	fprintf(fphtm, "    title: 'Frequency of Bases Quality Values', hAxis: { title: 'Quality value' }, vAxis: { title: 'Freq%%', format: 'percent' }, colors: ['#a52714', '#097138']\n");
	fprintf(fphtm, "    title: 'Frequency of base quality values', hAxis: { title: 'Quality value' }, vAxis: { title: 'Freq%%', format: 'percent' }, colors: ['#a52714', '#097138']\n");
	fprintf(fphtm, "  };\n");

	fprintf(fphtm, "  var d2 = new google.visualization.DataTable();\n");
	fprintf(fphtm, "  d2.addColumn('number', 'Quality value');\n");
	fprintf(fphtm, "  d2.addColumn('number', 'Freq');\n");
	fprintf(fphtm, "  d2.addRows( [ "); 

	for (size_t i=0; i<QSCORE_SIZE; i++)
	{
		double tmpf = (double)QCount[i] / (double)TotalLen;
		fprintf(fphtm, "[%lu,%.4f], ", i, tmpf);
		if (i % 10 == 9)
			fprintf(fphtm, "\n");
	}

	fprintf(fphtm, "  ] );\n");
	fprintf(fphtm, "  var chart2 = new google.visualization.AreaChart(document.getElementById('bq'));\n");
	fprintf(fphtm, "  chart2.draw(d2, opt2);\n");
	fprintf(fphtm, "\n");

// --- MinQ dist
/*	fprintf(fpout, "\n--- Read MinQ-value Count ---\n");
	fprintf(fpout, "MinQ,Count,Freq%%,cumuFreq%%\n");
	for (size_t i=0, cumuCnt=ReadCount; i<QSCORE_SIZE; i++)
	{
		if (MinQCount[i] > 0)
		{
			double tmpf = 100.0*(double)MinQCount[i] / (double)ReadCount;
			double tmpf2 = 100.0*(double)cumuCnt / (double)ReadCount;
			fprintf(fpout, "%lu,%lu,%.2f%%,%.2f%%\n", i, MinQCount[i], tmpf, tmpf2);
			cumuCnt -= MinQCount[i];
		}
	}
*/ 

	// option
	fprintf(fphtm, "  var opt3 = {\n");
	fprintf(fphtm, "    title: 'MinimalQ distribution', hAxis: { title: 'MinmalQ value', viewWindow: { max: 41 } }, vAxis: { title: 'Freq%%', format: 'percent' }, colors: ['#a52714', '#097138']\n");
	fprintf(fphtm, "  };\n");

	fprintf(fphtm, "  var d3 = new google.visualization.DataTable();\n");
	fprintf(fphtm, "  d3.addColumn('number', 'MinQ');\n");
	fprintf(fphtm, "  d3.addColumn('number', 'Freq');\n");
	fprintf(fphtm, "  d3.addRows( [ "); 

	for (size_t i=0, cumuCnt=ReadCount; i<QSCORE_SIZE; i++)
	{
		double tmpf = (double)MinQCount[i] / (double)ReadCount;
//		double tmpf2 = (double)cumuCnt / (double)ReadCount;
		fprintf(fphtm, "[%lu,%.4f], ", i, tmpf);
		if (i % 10 == 9)
			fprintf(fphtm, "\n");
		cumuCnt -= MinQCount[i];
	}

	fprintf(fphtm, "  ] );\n");
	fprintf(fphtm, "  var chart3 = new google.visualization.AreaChart(document.getElementById('mq'));\n");
	fprintf(fphtm, "  chart3.draw(d3, opt3);\n");
	fprintf(fphtm, "\n");

// --- Dist of %HighQ(x) 
	// option
	fprintf(fphtm, "  var optH = {\n");
	fprintf(fphtm, "    title: 'Coverage of reads with %%HighQ(q) >= X%%', hAxis: { title: 'X%% (X%% from 100 downto 50)', direction: -1, viewWindow: { max: 100, min: 50 } }, vAxis: { title: 'Coverage%%', format: 'percent' }, colors: ['#a52714', '#097138']\n");
	fprintf(fphtm, "  };\n");
	
	fprintf(fphtm, "  var dH = new google.visualization.DataTable();\n");
	fprintf(fphtm, "  dH.addColumn('number', 'X%%');\n");
	for (int k=0; k<HiQParamSize; k++)
		fprintf(fphtm, "  dH.addColumn('number', 'q=%d');\n", HiQTh[k]);
	fprintf(fphtm, "  dH.addRows( [ "); 

/*		fprintf(fpout, "\n--- PE %%HighQ(%d)-value Count ---\n", HiQTh[j]);
		fprintf(fpout, "%%HighQ(%d),Count,Freq%%,CumuFreq%%\n", HiQTh[j]);
*/
	for (long i=HiQCellSize, cumuCnt[HiQParamSize]={0}; i >= 0; i--)
	{
		// Output each record		
		fprintf(fphtm, "[%.1f", (double)i*(100.0/HiQCellSize));
		double cumuRatio[HiQParamSize];
		for (int k=0; k<HiQParamSize; k++)
		{
			cumuCnt[k] += HiQPercentCount[k][i];
			cumuRatio[k] = (double)cumuCnt[k] / (double)ReadCount;
			fprintf(fphtm, ",%.4ff", cumuRatio[k]);
		}
		fprintf(fphtm, "],");
//		double tmpf = (double)HiQPercentCount[j][i] / (double)ReadCount;
//			fprintf(fpout, "%.1f,%lu,%.2f%%,%.2f%%\n", (double)i*(100.0/HiQCellSize), HiQPercentCount[j][i], 100.0*tmpf, 100.0*tmpf2);

		if (i % 10 == 9)
			fprintf(fphtm, "\n");
	}
	
	fprintf(fphtm, "  [0");
	for (int k=0; k<HiQParamSize; k++)
		fprintf(fphtm, ",1");
	fprintf(fphtm, "] ] );\n");
	fprintf(fphtm, "  var chartH = new google.visualization.LineChart(document.getElementById('hq'));\n");
	fprintf(fphtm, "  chartH.draw(dH, optH);\n");
	fprintf(fphtm, "\n");
//	fprintf(fpout, "\n");


	// end of htm
	fprintf(fphtm, "}\n");
	fprintf(fphtm, "google.charts.load('current', {'packages': ['corechart', 'line']});\n");
	fprintf(fphtm, "google.charts.setOnLoadCallback(DrawDist);"); 
	fprintf(fphtm, "</script>\n");
	fprintf(fphtm, "</body>\n</html>\n");

//	fclose(fpout);
	fclose(fphtm);

	return true;
}



//=============================================================================
int main(int argc, char **argv)
{
//=============================================================================
	if(argc != 3)
	{
		printf("=== readQdist: Read a read FASTQ file and generate quality distribution and GC%% of the FASTQ file ===\n\n");
		printf("Usage: readQdist in.fq outPrjName\n");
//		printf("Output: outPrjName.htm, outPrjName.csv\n");
		printf("Output: outPrjName.htm\n");
		printf("Verson: 1.0 (2018/04) \n");
		printf("Author: Yu-Jung Chang\n\n");

		return 1;
	}

//=============================================================================
	ProbeFASTQ(argv[1], argv[2]);

//=============================================================================
	return 0;
}

