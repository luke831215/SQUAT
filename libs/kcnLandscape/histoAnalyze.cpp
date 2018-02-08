/*
== Author: Yu-Jung Chang 2016/03~
Read *.pfqc or *.pfqdc file and output records in text
*/
//=============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include <errno.h>

using namespace std;

//=============================================================================
#define LINE_BUF_SIZE 10000

class TangentIntervalRecord
{
public:
	char label; // -1, 0, +1
	size_t start, end;
	long inflect;
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

char GetTanLabel(long tangent)
{
	char curLabel;
	if (tangent > 0)
		curLabel = '+';
	else // <=0
		curLabel = '-';

	return curLabel;
}

//=============================================================================
bool histoAnalyze(char *input, char *outprefix)
{
	FILE *fpi, *fpo;
	char tmps[512], buf[LINE_BUF_SIZE];

	// Open files
	OpenFile(input, "rt", fpi);
	sprintf(tmps, "%s.html", outprefix);
	OpenFile(tmps, "wt", fpo);

/*	FILE *fdebug;
	sprintf(tmps, "%s-debug.txt", outprefix);
	OpenFile(tmps, "wt", fdebug);
*/
	// Load histo
	vector<size_t> histoAry;
	size_t NumOfAllKmer=0;
	while(fgets(buf, LINE_BUF_SIZE, fpi) != NULL)
	{
		strtok(buf, "\t\n"); // skip #copy 
		size_t num = (size_t)atol(strtok(NULL, "\t\n"));
		histoAry.push_back(num);
		NumOfAllKmer += num;
	}
	fclose(fpi);
	printf("In %s, #copy: from 1 to %lu\nTotal #DistinctKmer = %lu\n", input, histoAry.size(), NumOfAllKmer);
	
	// Output html
	fprintf(fpo, "<html>\n<head>\n<script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script></head>\n");

	fprintf(fpo, "<body>\n");
	fprintf(fpo, "<table border=1 width=100%%>\n");
	fprintf(fpo, "  <tr><td id=d1 style=\"width: 100%%; height: 600px\"></td></tr>\n");

	// script
	fprintf(fpo, "<script type=\"text/javascript\">\n");
	fprintf(fpo, "function DrawKmerHistogram() {\n");
	// option
	fprintf(fpo, "  var opt1 = {\n");
	fprintf(fpo, "    title: 'Kmer Histogram', hAxis: { title: 'log #Copy', scaleType: 'log' }, vAxis: { title: 'log #DistinctKmer', scaleType: 'log'}, colors: ['#a52714', '#097138']\n");
	fprintf(fpo, "  };\n");

	fprintf(fpo, "  var d1 = new google.visualization.DataTable();\n");
	fprintf(fpo, "  d1.addColumn('number', '#Copy');\n");
	fprintf(fpo, "  d1.addColumn('number', '#DistinctKmer');\n");
	fprintf(fpo, "  d1.addRows( [ "); 
	
	for (size_t i=0; i<histoAry.size(); i++)
	{
		fprintf(fpo, "[%lu,%lu], ", i+1, histoAry[i]); // i+1 to correct #copy as 1-based
		if (i % 10 == 9)
			fprintf(fpo, "\n    ");
	}
	
	fprintf(fpo, "  ] );\n");
	fprintf(fpo, "  var chart1 = new google.visualization.LineChart(document.getElementById('d1'));\n");
	fprintf(fpo, "  chart1.draw(d1, opt1);\n");
	fprintf(fpo, "}\n");
	
    fprintf(fpo, "google.charts.load('current', {packages: ['corechart', 'line']});\n");
    fprintf(fpo, "google.charts.setOnLoadCallback(DrawKmerHistogram);"); 
	fprintf(fpo, "</script>\n");

	// Detect the 1st foot, mountains & peaks, KCov 
	vector<long> tangentAry; // +/- 3; tangent = deltaY/deltaX = (y-y0)/(x-x0) = a[i+3] - a[i-3] / (i+3-(i-3)) = a[i+3] - a[i-3] / 6
	for (size_t i=0; i<3; i++) // just for sync of i without affecting the trend
	{
		tangentAry.push_back(((long)histoAry[6]-(long)histoAry[0])/6);
//		fprintf(fdebug, "%lu,%lu,%ld\n", i+1, histoAry[i], tangentAry[i]);
	}
	for (size_t i=3; i<histoAry.size()-4; i++) // -4: -3 and skip 10000
	{
		tangentAry.push_back( ((long)histoAry[i+3]-(long)histoAry[i-3])/6 );
//		fprintf(fdebug, "%lu,%lu,%ld\n", i+1, histoAry[i], tangentAry[i]);
	}

	vector<TangentIntervalRecord> tanIntervals;
	TangentIntervalRecord tiRec = {0};
	tiRec.label = GetTanLabel(tangentAry[3]);
	tiRec.start = 0;
	tiRec.end = 3;
	long idxMin=-1, idxMax=-1;
	long imin=0x7FFFFFFF, imax=0-imin;
	for (size_t i=0; i<tangentAry.size()-4; i++)
	{
		char lab = GetTanLabel(tangentAry[i]);
		if (lab == tiRec.label)
		{
			tiRec.end = i;
			if (tangentAry[i] < imin)
			{
				imin = tangentAry[i];
				idxMin = (long)i;
			}
			if (tangentAry[i] > imax)
			{
				imax = tangentAry[i];
				idxMax = (long)i;
			}
		}
		else
		{
			// add the old record
			if (tiRec.end > tiRec.start)
			{
				tiRec.inflect = -1;
				if (tiRec.label == '+' && idxMax != tiRec.end)
				{
					tiRec.inflect = idxMax;
				}
				if (tiRec.label == '-' && idxMin != tiRec.start)
				{
					tiRec.inflect = idxMin;
				}
				tanIntervals.push_back(tiRec);

				// new record
				tiRec.label = lab;
				tiRec.start = i;
				tiRec.end = i;
				imin=0x7FFFFFFF;
				imax=0-imin;
			}
			else // tiRec.end == tiRec.start
			{
				tiRec.end = i;
				if (tangentAry[i] < imin)
				{
					imin = tangentAry[i];
					idxMin = (long)i;
				}
				if (tangentAry[i] > imax)
				{
					imax = tangentAry[i];
					idxMax = (long)i;
				}
			}
		}
	}
/*
	fprintf(fdebug, "IntervalNo,Type,Start,End,Size,InflectionPoint\n");
	for (size_t i=0; i<tanIntervals.size(); i++)
	{
		tiRec = tanIntervals[i];
		size_t isize = tiRec.end - tiRec.start + 1;
		fprintf(fdebug, "%lu,%c,%lu,%lu,%lu,%ld\n", i+1, tiRec.label, tiRec.start+1, tiRec.end+1, isize, tiRec.inflect+1); // +1 to 1-based
	}
*/	
	fprintf(fpo, "  <tr><td id=d2>\n");
	if (tanIntervals.size() >= 2 && tanIntervals[0].label == '-' && tanIntervals[1].label == '+')
	{
		if (tanIntervals[0].inflect >= 0)
			fprintf(fpo, "foot point: %lu (slope decent from %lu; inflection point: %ld)<br>\n", tanIntervals[0].end+1, tanIntervals[0].start+1, tanIntervals[0].inflect+1);
		else
			fprintf(fpo, "foot point: %lu (slope decent from %lu; inflection point: none)<br>\n", tanIntervals[0].end+1, tanIntervals[0].start+1);
	}
	int cntM = 0;
	for (size_t i=1; i<tanIntervals.size()-1; i++)
	{
		tiRec = tanIntervals[i];
		size_t isize = tiRec.end - tiRec.start + 1;
		if (tiRec.label == '+' && isize >= 10)
		{
			size_t peakPos;
			if (histoAry[tiRec.end] > histoAry[tiRec.end+1])
				peakPos = tiRec.end;
			else
				peakPos = tiRec.end+1;

			fprintf(fpo, "mountain %d: up(%lu-%lu,%lu), peak %lu, down(%lu-%lu,%lu)<br>\n", ++cntM, tiRec.start+1, tiRec.end+1, tiRec.inflect+1,
				peakPos+1, tanIntervals[i+1].start+1, tanIntervals[i+1].end+1, tanIntervals[i+1].inflect+1);
		}
	}
	

	fprintf(fpo, "  </td></tr>\n");
	fprintf(fpo, "</table>\n");
	fprintf(fpo, "</body>\n</html>\n");

	fclose(fpo);
	
	return true;
}



//=============================================================================
int main(int argc, char **argv)
{
//=============================================================================
	if(argc != 3)
	{
		printf("=== histoAnalyze: Read *_histo.txt and ouput html(Google chart) ... ===\n\n");
		printf("Usage1: histoAnalyze Input_histo.txt OutputPrefix\n\n");
		printf("=== By Yu-Jung Chang, Last modified: 2017.04 ===\n\n");

		return 1;
	}

//=============================================================================
	histoAnalyze(argv[1], argv[2]);

//=============================================================================
	return 0;
}

