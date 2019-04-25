#include <iostream>
#include <fstream>
#include <unordered_set>
#include <Rcpp.h>


using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
void uniqueBar(String inputFile, String outputFile)
{
    string inFile = inputFile;
    string outFile = outputFile;
    // check input file format
    if(inFile.substr(inFile.size()-5,5) != "fasta" && inFile.substr(inFile.size()-5,5) != "fastq") 
    {
       Rcout << "Error! Input file should be .fasta or .fastq format." << endl;
       return;
    }

    int length = 0;
    std::unordered_set<std::string> myset;
    string line;
    ifstream myfile (inFile);
    ofstream outfile (outFile);

    if (myfile.is_open())
    {
      if(inFile.substr(inFile.size()-5,5) == "fastq")
      {
        while ( getline (myfile,line) )
        {
          getline (myfile,line);
          if(myset.find(line) == myset.end())
            myset.insert(line);
          getline (myfile,line);
          getline (myfile,line);
        }
        myfile.close(); 
      }
      if(inFile.substr(inFile.size()-5,5) == "fasta")
      {
        while ( getline (myfile,line) )
        {
          getline (myfile,line);
          if(myset.find(line) == myset.end())
            myset.insert(line);
        }
        myfile.close(); 
      }
    }

    int i = 0;
    for(const std::string& x: myset)
    {
      if(inFile.substr(inFile.size()-5,5) == "fastq")
      {
        outfile << "@fake ID" << x <<endl;
        outfile << x <<endl;
        outfile <<"+"<<endl;
        outfile << "fake phred score" <<endl;
      }
      if(inFile.substr(inFile.size()-5,5) == "fasta")
      {
        outfile << ">" << i <<endl;
        outfile << x <<endl;
        i ++;
      }
    }

    return;
}


