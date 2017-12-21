
/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */

#include <tools/tools.hpp>
#include "ic_datafile.hpp"
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

double IC_datafile::DipoleAmplitude(double r, double b)
{
	if (b>1e-10)
		cerr << "Impact parameter is not supported! " << LINEINFO << endl;
		
	return interpolator->Evaluate(r);
}

/*
 * Load amplitude data from a given file
 * Syntax: r amplitude 
 * Lines starting with "#" are comments
 * 
 * If an error occurs, returns -1, otherwise returns 0
 */
int IC_datafile::LoadFile(std::string file)
{
    fname=file;
    std::ifstream f(file.c_str());
    if (!f.is_open())
    {
        cerr << "ERROR! Coudn't read file " << file << " " << LINEINFO << endl;
        return -1;
    }
    std::vector<double> rvals;
    std::vector<double> nvals;	// amplitude
    while(!f.eof())
    {
        string line;
        getline(f, line);
        if (line.substr(0, 1)=="#" or line.length() < 3)
			continue;
		std::stringstream ss(line);
		std::string r,n;
		ss >> r;
		ss >> n;
		rvals.push_back(StrToReal(r));
		nvals.push_back(StrToReal(n));
    }
    f.close();

    // If there wasn't enough data points
    if (nvals.size() < 10)
    {
        cerr << "There was only " << nvals.size() << " datapoints in file " << file <<", exiting..." << endl;
        exit(1);
    }
    
    // Interpolator
    interpolator = new Interpolator(rvals, nvals);
    interpolator->Initialize();
    interpolator->SetFreeze(true);
    interpolator->SetUnderflow(0.0);
    interpolator->SetOverflow(1.0);
    
    return 0;
}

std::string IC_datafile::GetString()
{
    std::stringstream ss;
    ss << "datafile " << fname << ", number of points " << interpolator->GetNumOfPoints() << " minr " << interpolator->MinX() << " maxr " << interpolator->MaxX() << endl;
    return ss.str();
}

IC_datafile::IC_datafile()
{
	interpolator=NULL;
}

IC_datafile::~IC_datafile()
{
	if (interpolator!=NULL)
		delete interpolator;
	
}

double IC_datafile::MinR()
{
	return interpolator->MinX();
}

double IC_datafile::MaxR()
{
	return interpolator->MaxX();
}
