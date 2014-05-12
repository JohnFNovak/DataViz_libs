#ifndef GNUPLOT_H_
#define GNUPLOT_H_

#include <string>
#include <iostream>

using namespace std;

class Gnuplot {
	public:
		Gnuplot();
		~Gnuplot();
		void operator ()(const string & command);
		// send any command to gnuplot
	protected:
		FILE *gnuplotpipe;
};

Gnuplot::Gnuplot() {
	// with -persist option you will see the windows as your program ends
	gnuplotpipe=popen("gnuplot -persist","w");
	if (!gnuplotpipe) { cerr<< ("Gnuplot not found !"); }
	fprintf(gnuplotpipe,"set term x11 noraise\n");
}

Gnuplot::~Gnuplot() {
	fprintf(gnuplotpipe,"exit\n");
	pclose(gnuplotpipe);
}

void Gnuplot::operator()(const string & command) {
	fprintf(gnuplotpipe,"%s\n",command.c_str());
	fflush(gnuplotpipe);
	// flush is necessary, nothing gets plotted else
};

#endif
