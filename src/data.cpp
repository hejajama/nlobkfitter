
/*
 * Loads HERA DIS data on reduced cross section
 * Keeps also track whether the data is for reduced cross section
 * or for only charm
 *
 * Datafile format is
 * Q^2 x y sigmared error
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 */
#include "data.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

double DataStrToReal(std::string str)
{
    std::stringstream buff(str);
    double tmp;
    buff >> tmp;
    return tmp;
}

// Load data from given file, return 0 if no errors
int Data::LoadData(string filename, DataType type)
{
    
    ifstream file(filename.c_str());
    if (!file.is_open())
    {
        cerr << "ERROR! Coudn't read file " << filename << endl;
        exit(-1);
    }
    
    int points=0;
    
    bool onlycharm, onlyFL = false;
    if (type == CHARM )
        onlycharm = true;
    else if (type == FL )
        onlyFL = true;

    while(!file.eof() )
    {
        string line;
        getline(file, line);
        if (line.substr(0, 1)=="#" or line.length() < 5)
        continue;
        string x,qsqr,y,sigmar,err;
        stringstream l(line);
        if (onlyFL == true){
            l >> qsqr; l>>x; l>>sigmar; l>>err;

            if (DataStrToReal(x)>maxx or DataStrToReal(x)<minx or DataStrToReal(qsqr)<minQ2 or DataStrToReal(qsqr)>maxQ2) continue;
            points++;
            Qsqrvals.push_back(DataStrToReal(qsqr)); xbjvals.push_back(DataStrToReal(x));
            sigmarvals.push_back(DataStrToReal(sigmar)); errors.push_back(DataStrToReal(err));
        }
        else{
            l >> qsqr; l>>x; l>>y; l>>sigmar; l>>err;

            if (DataStrToReal(x)>maxx or DataStrToReal(x)<minx or DataStrToReal(qsqr)<minQ2 or DataStrToReal(qsqr)>maxQ2) continue;
            points++;
            Qsqrvals.push_back(DataStrToReal(qsqr)); xbjvals.push_back(DataStrToReal(x));
            yvals.push_back(DataStrToReal(y));
            sigmarvals.push_back(DataStrToReal(sigmar)); errors.push_back(DataStrToReal(err));
            only_charm.push_back(onlycharm);
        }
	//cout << "Accept line " << line << endl;
    }
    
    cout << "# Loaded " << points << " datapoints from " << filename << " in Q2 range " << minQ2 << " - " << maxQ2 << " GeV^2, now we have in total " << sigmarvals.size() << " points " << endl;
    
    return 0;
}
int Data::NumOfPoints() const
{
    return sigmarvals.size();
}
                   
double Data::Qsqr(unsigned int n) const
{
    if (n >= NumOfPoints())
    {
        cerr << "Index " << n << " out of range!" << endl;
        exit(1);
    }
    return Qsqrvals[n];
}
 double Data::xbj(unsigned int n) const
{
    return xbjvals[n];
}
 double Data::y(unsigned int n) const
{
    return yvals[n];
}
 double Data::ReducedCrossSection(unsigned int n) const
{
    return sigmarvals[n];
}
double Data::ReducedCrossSectionError(unsigned int n) const
{
    return errors[n];
}

bool Data::OnlyCharm(unsigned int n) const
{
		    return only_charm[n];
}

Data::Data()
{
    minx=1e-99;
    maxx=0.01;
    minQ2=0.75;
    maxQ2=650;
    weight = 1.0;
}
