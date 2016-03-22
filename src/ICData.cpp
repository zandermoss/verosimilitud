#include "ICData.h"
#include <iostream>
#include <string>


void ICData::OpenCSV(std::string filename)
{
	datafile.open(filename.c_str()); 
}

void ICData::ReadCSV(std::vector<double>* cosz, std::vector<double>* eprox)
{
	std::string value;
	while ( datafile.good() )
	{
	     std::getline( datafile, value, '\n' );
	
		std::vector<std::string> tokens;
		std::string token;
		while(token != value){
		  token = value.substr(0,value.find_first_of(" "));
		  value = value.substr(value.find_first_of(" ") + 1);
		  tokens.push_back(token);
		}
	
		cosz->push_back(atof(tokens[tokens.size()-1].c_str()));
		eprox->push_back(atof(tokens[tokens.size()].c_str()));
	}
}


ICData::~ICData()
{	
	datafile.close();
}

