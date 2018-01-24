#include <iostream>
#include <fstream>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
void trim(String inputFile, String outputFile, int start, int end)
{
    string inFile = inputFile;
    string outFile = outputFile;
    // check input file format
    if(inFile.substr(inFile.size()-5,5) != "fastq") 
    {
       Rcout << "Error! Input file should be .fastq format." << endl;
       return;
    }

    string line;
    ifstream myfile (inFile);
    ofstream outfile (outFile);
    
    int length = end - start + 1;
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
        outfile << line <<endl;
        getline (myfile,line);
        outfile << line.substr(start-1, length) <<endl;
        getline (myfile,line);
        outfile << line <<endl;
        getline (myfile,line);
        outfile << line.substr(start-1, length) <<endl;
      }
      myfile.close();
    }

    return;
}


