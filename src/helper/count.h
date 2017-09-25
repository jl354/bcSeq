//This is function to update count table for the reads
#ifndef __count_h__
#define __count_h__

#include "../structures/Trie.h"
#include "fastqCheck.h"
#include <fstream>
#include <iostream>
#include <map>
#include <string>

using namespace std;

static bool readFastq(istream& in,
    vector<string>& library,
    vector<string>& library_ids)
{
    string id;
    string line;

    while (getline(in, id) && getline(in, line)) {
        auto t = find(library.begin(),
            library.end(),
            line);

        if (t != library.end()) {
            Rcpp::Rcerr << "Library Error: Repeat entry in library file\n"
                 << line << endl;
            return false;
        } else {

            if (!checkRead(line)) {
                Rcpp::Rcout << "Library Error: invalid characters in read\n"
                     << "ID: " << id << "\n"
                     << "read: " << line << "\n"
                     << endl;
                return false;
            }
            id = id.substr(0, id.find_first_of(" "));

            library.push_back(move(line));
            library_ids.push_back(move(id));
        }

        if (!(getline(in, line) && getline(in, line))) {
            Rcpp::Rcout << "Error encountered reading library file" << endl;
            return false;
        }
    }
    return in.eof();
}

static bool readFasta(istream& in,
    vector<string>& library,
    vector<string>& library_ids)
{
    string id;
    string line;

    while (getline(in, id) && getline(in, line)) {
        auto t = find(library.begin(),
            library.end(),
            line);

        if (t != library.end()) {
            Rcpp::Rcerr << "Library Error: Repeat entry in library file\n"
                 << line << endl;
            return false;
        } else {

            if (!checkRead(line)) {
                Rcpp::Rcout << "Library Error: invalid characters in read\n"
                     << "ID: " << id << "\n"
                     << "read: " << line << "\n"
                     << endl;
                return false;
            }
            id = id.substr(1, id.length());

            library.push_back(move(line));
            library_ids.push_back(move(id));
        }
    }
    return in.eof();
}

static bool readLibrary(vector<string>& library,
    vector<string>& library_ids,
    const string& libFile)
{
    try {
        ifstream in(libFile);

        if (in.peek() == '@') {
            return readFastq(in, library, library_ids);
        } else if (in.peek() == '>') {
            return readFasta(in, library, library_ids);
        } else {
            return false;
        }

    } catch (Rcpp::exception& e) {
        Rcpp::Rcout << e.what() << endl;
        return false;
    }
}

static void count2CSV(const vector<double>& countTable,
    const vector<string>& library,
    const string& csvFileName,
    const ios::openmode& mode = ios::out | ios::app)
{
    ofstream myfile(csvFileName, mode);

    for (auto x = countTable.begin(); x < countTable.end(); x++) {
        myfile << library[x - countTable.begin()] << "," << *x << "\n";
    }

    if (!myfile) {
        Rcpp::Rcout << "Error writing results to " << csvFileName << endl;
    }
}

#endif
