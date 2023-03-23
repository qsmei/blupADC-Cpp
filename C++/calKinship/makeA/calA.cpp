#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

arma::Mat<double> makeA_cpp(arma::Mat<int>& Pedigree){  

 cout<<"Start constructing pedigree additive relationship matrix......"<<endl;    
 arma::Mat<double> A_matrix(Pedigree.n_rows,Pedigree.n_rows,fill::eye); 
 arma::Col<int> Sire=Pedigree.col(1);
 arma::Col<int> Dam=Pedigree.col(2);
 int IND=Pedigree.n_rows;
 bool status1,status2;
 
 for(int i=0; i < IND; i++){
 status1=(Sire[i]==0);
 status2=(Dam[i]==0);

 if(!status1&&!status2){  // both parents known
   A_matrix(i,i)=1+0.5*(A_matrix(Sire[i]-1,Dam[i]-1)); 
 for(int j=0;j<i;j++){
 if(i>j) A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Sire[i]-1)+A_matrix(j,Dam[i]-1));
   }
  }else if(!status1&&status2){  //known sire，unknown dam 

 for(int j=0;j<i;j++){
 if(i>j) A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Sire[i]-1));
   }
  }else if(status1&&!status2){  //known dam，unknown sire 

 for(int j=0;j<i;j++){
 if(i>j) A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Dam[i]-1));
   }
  }
 }
 cout<<"Complete constructing pedigree additive relationship matrix!"<<endl;   
 return A_matrix;
 }
	
int main(int argc, char* argv[]) {
    if (argc != 5 || string(argv[1]) != "--ped" || string(argv[3]) != "--out") {
        cerr << "Usage: " << argv[0] << " --ped <input_ped_file> --out <output_Amatrix_file>" << endl;
        return 1;
    }

    // Read pedigree file
    string pedigree_file = argv[2];
    ifstream in_pedigree(pedigree_file);

    vector<vector<int>> pedigree_data;
    int id, sire, dam;
    while (in_pedigree >> id >> sire >> dam) {
        pedigree_data.push_back({id, sire, dam});
    }

    //Mat<int> Pedigree = conv_to<Mat<int>>::from(pedigree_data);

    // Create a Mat with the same dimensions as the vector
    int rows = pedigree_data.size();
    int cols = pedigree_data[0].size();
    Mat<int> Pedigree(rows, cols);

    // Copy the data from the vector to the matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            Pedigree(i, j) = pedigree_data[i][j];
        }
    }

    // Compute A matrix
    mat A_matrix = makeA_cpp(Pedigree);
    // Write A matrix to file
    string A_matrix_file = argv[4];
    A_matrix.save(A_matrix_file, raw_ascii);

    return 0;
}
