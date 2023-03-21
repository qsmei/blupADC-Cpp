#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <future>
#include <cstring>
#include <sys/resource.h>
#include <omp.h>
#include "ProgressBar.hpp"

using namespace std;

vector<vector<string>> hapmap_to_ped(const vector<vector<string>>& data_hmp) {
    int n_ind = data_hmp[0].size() - 11, n_snp = data_hmp.size()-1, i, j;
    vector<string> temp_vec1, IND_name=data_hmp[0];
	cout<<"n_SNP="<<n_snp<<endl;
    vector<vector<string>> data_ped(n_ind, vector<string>(n_snp * 2 + 6));
    cout << "Start Hapmap to Plink data format conversion......" << endl;
	
	// Initialize the progress bar
    progresscpp::ProgressBar progressBar(n_ind,50);
    #pragma omp parallel for shared(data_ped, data_hmp, n_ind, n_snp, progressBar) private(i, j)
    for (i = 0; i < n_ind; i++) {
        for (j = 0; j < n_snp; j++) {
            data_ped[i][6 + 2 * j] = data_hmp[j+1][i + 11][0];
            data_ped[i][6 + 2 * j + 1] = data_hmp[j+1][i + 11][1];
        }
        data_ped[i][0] = IND_name[i+11];
        data_ped[i][1] = IND_name[i+11];
        data_ped[i][2] = "0";
        data_ped[i][3] = "0";
        data_ped[i][4] = "0";
        data_ped[i][5] = "0";

         // record the tick
        ++progressBar;

        // display the bar
        progressBar.display();
    }

    // Finish the progress bar
     progressBar.done();

    cout << "Complete Hapmap to Plink data format conversion!" << endl;
    return data_ped;
}



// Function to read Hapmap format data from file
vector<vector<string>> read_hapmap_file(string filename) {
    vector<vector<string>> data;
    ifstream file(filename);
    if (file.is_open()) {
		cout << "Start read hapmap data......" << endl;	
		auto start = chrono::high_resolution_clock::now(); // measure the start time	
        string line;
        while (getline(file, line)) {
            vector<string> row;
            stringstream ss(line);
            string cell;
            while (getline(ss, cell, '\t')) {
                row.push_back(cell);
            }
            data.push_back(row);
        }
        file.close();
		auto end = chrono::high_resolution_clock::now(); // measure the end time
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);		
		cout << "Complete read hapmap data!" << endl;	
		cout << "Time taken to read file: " << duration.count() << "ms" << endl; // print the time taken to read the file			
    }
    else {
        cerr << "Unable to open file " << filename << endl;
    }
    return data;
}

// Function to write Plink format data to file
void write_ped_file(const vector<vector<string>>& data_ped, string filename) {
    ofstream file(filename);
    if (file.is_open()) {
        for (auto row : data_ped) {
            for (auto cell : row) {
                file << cell << "\t";
            }
            file << endl;
        }
        file.close();
    }
    else {
        cerr << "Unable to open file " << filename << endl;
    }
}


int main(int argc, char* argv[]) {
  int num_threads = 1;
  string hapmap_filename;
  string ped_filename;

  // Parse command line arguments
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-threads") {
      if (i + 1 < argc) {
        num_threads = atoi(argv[i+1]);
        omp_set_num_threads(num_threads);
      } else {
        cerr << "Error: -threads option requires a number argument" << endl;
        return 1;
      }
    } else if (string(argv[i]) == "-hmp") {
      if (i + 1 < argc) {
        hapmap_filename = string(argv[i+1]);
      } else {
        cerr << "Error: -hmp option requires a filename argument" << endl;
        return 1;
      }
    } else if (string(argv[i]) == "-ped") {
      if (i + 1 < argc) {
        ped_filename = string(argv[i+1]);
      } else {
        cerr << "Error: -ped option requires a filename argument" << endl;
        return 1;
      }
    }
  }

  if (hapmap_filename.empty() || ped_filename.empty()) {
    cerr << "Usage: " << argv[0] << " -hmp hapmap_file -ped ped_file [-threads num_threads]" << endl;
    return 1;
  }

  rusage usage_start, usage_end;
  getrusage(RUSAGE_SELF, &usage_start);

  vector<vector<string>> data_hmp = read_hapmap_file(hapmap_filename);
  vector<vector<string>> data_ped = hapmap_to_ped(data_hmp);
  write_ped_file(data_ped, ped_filename);

  getrusage(RUSAGE_SELF, &usage_end);
  long memory_usage = usage_end.ru_maxrss - usage_start.ru_maxrss;
  cout << "Memory usage: " << memory_usage << " kilobytes" << endl;

  return 0;
}
	