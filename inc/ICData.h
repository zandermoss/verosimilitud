#ifndef __ICDATA_H_INCLUDED__
#define __ICDATA_H_INCLUDED__



#include <string>
#include <vector>
#include <iostream>
#include <fstream>

class ICData
{
	public:
		ICData();
		~ICData();

		void ReadCSV(std::vector<double>* cosz, std::vector<double>* eprox);
		void OpenCSV(std::string filename);

	private:
		std::ifstream datafile;
};


#endif // __ICDATA_H_INCLUDED__

