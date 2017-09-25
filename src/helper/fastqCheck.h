#ifndef __fastqCheck_h__
#define __fastqCheck_h__

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <istream>
#include <memory>
#include <string>

using namespace std;

// true if no error
static auto noError(const char& x) -> bool
{
    return x == 'A' || x == 'C' || x == 'G' || x == 'T' || x == 'N';
}

// true if there are no errors found
static auto checkRead(const std::string& read) -> bool
{
    return all_of(read.begin(), read.end(), noError);
}

// true if file is read successfully
static auto readSamples(const string& sampleFile,
    vector<string>& sequence,
    vector<string>& sequence_ids,
    vector<string>& phredScore)
    -> bool
{
    try {
        ifstream in(sampleFile);
        string id;
        string line;
        string phred;

        while (getline(in, id) && getline(in, line) && getline(in, phred) && getline(in, phred)) {
            if (!(line.length() == phred.length())) {
                Rcpp::Rcout << "Sample Error: sequence length differs from Phred string length\n"
                     << "ID: " << id << "\n"
                     << "read: " << line << "\n"
                     << "Phred: " << phred
                     << endl;
                return false;
            }

            if (!checkRead(line)) {
                Rcpp::Rcout << "Sample Error: invalid characters in read\n"
                     << "id: " << id << "\n"
                     << "read: " << line << "\n"
                     << endl;
                return false;
            }
            id = id.substr(0, id.find_first_of(" "));

            sequence.push_back(move(line));
            sequence_ids.push_back(move(id));
            phredScore.push_back(move(phred));
        }
        return in.eof();

    } catch (Rcpp::exception& e) {
        Rcpp::Rcout << e.what() << endl;
        return false;
    }
}

static auto phred2err(vector<double>& cErr, string phred) -> void
{
    cErr.clear();
    transform(phred.begin(),
        phred.end(),
        back_inserter(cErr),
        [](const char x) {
            return pow(10.0, (33 - x) / 10.0);
        });
}

#endif
