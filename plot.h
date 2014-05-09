// plot.h 

#ifndef PLOT_H_
#define PLOT_H_

#include <string>
#include <iostream>
#include <vector>
#include "gnuplot.h"
#include "output.h"

using namespace std;

class Plot : Gnuplot {
	public:
		Plot();
		~Plot();
		void Divide(int, int);
		void Cd(int);
		void Title(const string &);
		void Data(spectra*);
		void Draw();
	protected:
		int dims[2];
		int size;
		int index;
		vector<spectra*> data;
		vector<string> titles;
};

Plot::Plot(){
	index = 0;
	dims[0] = 1;
	dims[1] = 1;
	size = 1;
	data.resize(1);
	titles.resize(1, "\"Title\"");
};

Plot::~Plot(){
};

void Plot::Divide(int columns, int rows) {
	dims[0] = columns;
	dims[1] = rows;
	size = columns * rows;
	data.resize(size);
	titles.resize(size, "");
	for(int i=0; i<size; i++) {
		char buff[200];
		sprintf(buff, "Title_%d", i);
		titles[i] = buff;
	}
};

void Plot::Cd(int cell) {
	if (cell > (size) ) {
		cerr << "Cd called for cell too large: " << cell << " dimensions: " << dims[0] << ", " << dims[1] << endl;
	}
	else {
		index = cell;
	}
};

void Plot::Title(const string &title) {
	titles[index] = title;
};

void Plot::Data(spectra *S) {
	data[index] = S;
};

void Plot::Draw() {
	if (size > 1){
		fprintf(gnuplotpipe, "set multiplot layout %d, %d\n", dims[0], dims[1]);
	}
	for(int i=0; i<size; i++) {
		fprintf(gnuplotpipe, "plot \"-\" title \"%s\" with lines\n", titles[i].c_str());
		double Xwidth = (data[i]->getUB1()-data[i]->getLB1())/data[i]->getBins1();
		double X = data[i]->getLB1();
		int bins = data[i]->getBins1();
		for(int j=0; j<bins; j++) {
			fprintf(gnuplotpipe, "%f %f\n", X, data[i]->getVal(j));
			X += Xwidth;
		}
		fprintf(gnuplotpipe, "e\n#\n");
	}
	fprintf(gnuplotpipe, "unset multiplot\n");
	fflush(gnuplotpipe);
};

#endif