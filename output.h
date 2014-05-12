// output.h

// This is a collection of classes designed to make creating output easier

// John F. Novak
// June 26, 2012

#ifndef _OUTPUT_
#define _OUTPUT_

#include <math.h>
#include <iostream>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
// #include <data.h>
// #include <cuts.h>
#include <string>
#include <vector>

using namespace std;

class spectra {
    public:
        spectra(double, double, int);  // takes: lower, upper, number of bins
        spectra(double, double, int, double, double, int);
        template<class TYPE>
        void addEntry(TYPE);  // takes: entry value

        template<class TYPE_A, class TYPE_B>
        void addEntry(TYPE_A, TYPE_B);

        void print(string);  // takes: output file name
        double getMean(void);
        double getXMean(int);
        double getYMean(int);
        double getStandardDeviation(void);
        double getXStandardDeviation(int);
        double getYStandardDeviation(int);
        void zero(void);
        double getVal(int);
        double getVal(int, int);
        int length(void);
        int height(void);
        double getLB(void);
        double getUB(void);
        int getBins(void);
        double getLB1(void);
        double getUB1(void);
        int getBins1(void);
        double getLB2(void);
        double getUB2(void);
        int getBins2(void);
        int D(void);
        void set(int, double);
        void set(int, int, double);
        void setCount(int);
        int Count(void);
        void sum(spectra, spectra);
        int count;

    private:
        bool OneD;
        double Xlb, Ylb, Xub, Yub;
        int Xbins, Ybins;
        vector<double> spectra1;
        vector< vector<double> > spectra2;
};

spectra::spectra(double l_bound, double u_bound, int bins) {
    Xlb = l_bound;
    Xub = u_bound;
    Xbins = bins;
    for (int i = 0; i < Xbins; i++) {
        spectra1.push_back(0);
    }
    count = 0;
    OneD = true;
}

spectra::spectra(double l_bound1, double u_bound1, int bins1,
                 double l_bound2, double u_bound2, int bins2) {
    Xlb = l_bound1;
    Xub = u_bound1;
    Xbins = bins1;
    Ylb = l_bound2;
    Yub = u_bound2;
    Ybins = bins2;
    OneD = false;
    vector<double> temp;
    for (int i = 0; i < Xbins; i++) {
        temp.push_back(0);
    }
    for (int i = 0; i < Ybins; i++) {
        spectra2.push_back(temp);
    }
    count = 0;
}

template<class TYPE>
void spectra::addEntry(TYPE val) {
    if ((val < Xub) && (val > Xlb)) {
        spectra1[int(((val-Xlb)/(Xub-Xlb))*Xbins)]++;
        count++;
    }
}

template<class TYPE_A, class TYPE_B>
void spectra::addEntry(TYPE_A Xval, TYPE_B Yval) {
    if ((Xval < Xub) && (Xval > Xlb) && (Yval < Yub) && (Yval > Ylb)) {
        spectra2[int(((Yval-Ylb)/(Yub-Ylb))*Ybins)][int(((Xval-Xlb)/(Xub-Xlb))*Xbins)]++;
        count++;
    }
}

void spectra::print(string outputfile) {
    fstream output;
            float Xwidth, Ywidth;
    output.open(outputfile.c_str(), fstream::out);
    if (OneD) {
                    Xwidth = (Xub-Xlb)/Xbins;
        for (int i = 0; i < Xbins; i++) {
            output << ((i*Xwidth)+Xlb+(Xwidth/2)) << " ";
        }
        output << endl;
        for (int i = 0; i < Xbins; i++) {
            output << spectra1[i] << " ";
        }
    }
    if (!OneD) {
        output << 0 << " ";
                    Xwidth = (Xub-Xlb)/Xbins;
                    Ywidth = (Yub-Ylb)/Ybins;
        for (int i = 0; i < Xbins; i++) {
            output << ((i*Xwidth)+Xlb+(Xwidth/2)) << " ";
        }
        output << endl;
        for (int j = 0; j < Ybins; j++) {
            output << ((j*Ywidth)+Ylb+(Ywidth/2)) << " ";
            for (int i = 0; i < Xbins; i++) {
                output << spectra2[j][i] << " ";
            }
            output << endl;
        }
    }
    output.close();
}

double spectra::getMean(void) {
    if (!OneD) {
        cout << "getMean() is only defined for 1-D spectra." << endl;
        return 0;
    }
    double sum = 0;
        double tcount = 0;
    for (int i = 0; i < Xbins; i++) {
        sum += i*spectra1[i];
                tcount += i;
    }

    return (sum/tcount)*(Xub-Xlb)+Xlb;
}

double spectra::getXMean(int bin) {
    if (OneD) {
        cout << "getXMean( int ) is only defined for 2-D spectra." << endl;
        return 0;
    }
    if ( bin > Xbins ) {
        cout << "bin number given is greater than the number of bins." << endl;
        return 0;
    }
    double sum = 0;
        double tcount = 0;
    for (int i = 0; i < Ybins; i++) {
        sum += i*spectra2[i][bin];
        tcount += i;
    }
    return (sum/tcount)*(Xub-Xlb)+Xlb;
}

double spectra::getYMean(int bin) {
    if (OneD) {
        cout << "getYMean( int ) is only defined for 2-D spectra." << endl;
        return 0;
    }
    if ( bin > Ybins ) {
        cout << "bin number given is greater than the number of bins." << endl;
        return 0;
    }
    double sum = 0;
    double tcount = 0;
    for (int i = 0; i < Xbins; i++) {
        sum += i*spectra2[bin][i];
        tcount += i;
    }

    return (sum/tcount)*(Yub-Ylb)+Ylb;
}

double spectra::getStandardDeviation(void) {
    if (!OneD) {
        cout << "getStandardDeviation() is only defined for 1-D spectra." << endl;
        return 0;
    }
    double mean = (getMean()-Xlb)/(Xub-Xlb);
    double sum = 0;
    double tcount = 0;
    for (int i = 0; i < Xbins; i++) {
        sum += (i*(spectra1[i] - mean)*(spectra1[i] - mean));
        tcount += i;
    }

    return sqrt(sum/tcount)*(Xub-Xlb);
}

double spectra::getXStandardDeviation(int bin) {
    if (OneD) {
        cout << "getXStandardDeviation( int ) is only defined for 2-D spectra." << endl;
        return 0;
    }
    double mean = (getXMean(bin) - Xlb) / (Xub - Xlb);
    double sum = 0;
        double tcount = 0;
    for (int i = 0; i < Ybins; i++) {
        sum += (i*(spectra2[i][bin] - mean)*(spectra2[i][bin] - mean));
        tcount += i;
    }

    return sqrt(sum/tcount)*(Xub-Xlb);
}

double spectra::getYStandardDeviation(int bin) {
    if (OneD) {
        cout << "getYStandardDeviation( int ) is only defined for 2-D spectra." << endl;
        return 0;
    }
    double mean = (getYMean(bin) - Ylb) / (Yub - Ylb);
    double sum = 0;
        double tcount = 0;
    for (int i = 0; i < Xbins; i++) {
        sum += (i*(spectra2[bin][i] - mean)*(spectra2[bin][i] - mean));
        tcount += i;
    }

    return sqrt(sum/tcount)*(Yub - Ylb);
}

void spectra::zero(void) {
    if (OneD) {
        for (int i = 0; i < Xbins; i++) {
            spectra1[i] = 0;
        }
    }
    if (!OneD) {
        for (int j = 0; j < Ybins; j++) {
            for (int i = 0; i < Xbins; i++) {
                spectra2[j][i] = 0;
            }
        }
    }
}

double spectra::getVal(int i) {
    return spectra1[i];
}

double spectra::getVal(int i, int j) {
    return spectra2[i][j];
}

double spectra::getLB(void) {
    return Xlb;
}

double spectra::getUB(void) {
    return Xub;
}

int spectra::getBins(void) {
    return Xbins;
}

double spectra::getLB1(void) {
    return Xlb;
}

double spectra::getUB1(void) {
    return Xub;
}

int spectra::getBins1(void) {
    return Xbins;
}

double spectra::getLB2(void) {
    if (OneD) {
        cout << "getLB2 is only defined for 2-D spectra" << endl;
        return 0;
    }
    return Ylb;
}

double spectra::getUB2(void) {
    if (OneD) {
        cout << "getUB2 is only defined for 2-D spectra" << endl;
        return 0;
    }
    return Yub;
}

int spectra::getBins2(void) {
    if (OneD) {
        cout << "getBins2 is only defined for 2-D spectra" << endl;
        return 0;
    }
    return Ybins;
}

int spectra::length(void) {
    if (!OneD) {
        return spectra2.size();
    }
    if (OneD) {
        return spectra1.size();
    }
        return 0;
}

int spectra::height(void) {
    if (OneD) {
        cout << "height is only defined for 2-D spectra" << endl;
        return 0;
    }
    return spectra2[0].size();
}

int spectra::D(void) {
    if (OneD) {
        return 1;
    }
    if (!OneD) {
        return 2;
    }
        return 0;
}

void spectra::set(int i, double val) {
    spectra1[i] = val;
}

void spectra::set(int i, int j, double val) {
    spectra2[i][j] = val;
}

void spectra::setCount(int Count) {
    count = Count;
}

int spectra::Count(void) {
    return count;
}

void spectra::sum(spectra A, spectra B) {
    if ((A.D() == B.D()) && (A.length() == B.length())) {
        if (A.D() == 1) {
            for (int i = 0; i < A.length(); i++) {
                spectra1[i]=(A.getVal(i)+B.getVal(i));
            }
            count=(A.Count()+B.Count());
        }
        if ((A.D() == 2) && (A.height() == B.height())) {
            for (int i = 0; i < A.length(); i++) {
                for (int j = 0; j < A.height(); j++) {
                    spectra2[i][j]=(A.getVal(i, j)+B.getVal(i, j));
                }
            }
            count=(A.Count()+B.Count());
        }
    }
}
#endif
