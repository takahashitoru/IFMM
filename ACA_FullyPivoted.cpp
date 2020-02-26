#include "ACA_FullyPivoted.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

void ACA_FullyPivoted(MatrixXd &A,MatrixXd &U, MatrixXd &V,double &epsilon,int &rank,int minRank){

  int maxRank = min(A.rows(),A.cols());
  // cout << "A.rows()" << A.rows() << endl;
  // cout << "A.cols()" << A.cols() << endl;

  // bool minmax=false;
  // if (minRank>maxRank){
  // 	cout << "minmax" << endl;
  // 	printf("rank, maxRank, minRank: %d %d %d \n",rank, maxRank, minRank);
  // 	minmax=true;
  // 	// minRank=maxRank;
  // }
  // bool testnan=false;

  MatrixXd R = A;
  U = MatrixXd(A.rows(),0);
  V = MatrixXd(A.cols(),0);

  VectorXi Max_abs_col_index(A.cols());
  VectorXd Max_abs_col(A.cols());

  double Max;
  double gamma;

  double Anorm=A.norm();
  double Rnorm=R.norm();

  // double MinPivot = 1e-17;
  const double machine_eps = std::numeric_limits<double>::epsilon();
  double MinPivot = machine_eps;

   // if (isnan(Anorm)){
	  // cout << "A" << A << endl;
	// }

  // cout << "test" << endl;

  if (rank==0){

  	 // printf("rank, maxRank, Rnorm,Anorm,minRank: %d %d %f %f %d \n",rank, maxRank, Rnorm,Anorm,minRank);
  	 
  	 // printf("isnan(Anorm):%d\n",isnan(Anorm));
  	 // if (isnan(Anorm)){
  	 	// cout << "A" << A << endl;
  	 // }


  	 // printf("Rnorm,Anorm: %20.5e %20.5e \n",Rnorm,Anorm);

	  // while ((rank < maxRank && Rnorm>epsilon*Anorm && Rnorm > machine_eps) || rank<minRank){
	  while ((rank < maxRank && Rnorm>epsilon && Rnorm > machine_eps) || rank<minRank){


		// printf("Rank: %d\n",rank);


	  	U.conservativeResize(Eigen::NoChange,rank+1);
	  	V.conservativeResize(Eigen::NoChange,rank+1);



		for (int iCol = 0; iCol < A.cols(); iCol++){
		  Max_abs_col(iCol) = R.col(iCol).cwiseAbs().maxCoeff(&Max_abs_col_index(iCol));
		}
		
	    int sel_row, sel_col;
		double Max_abs = Max_abs_col.maxCoeff(&sel_col);
		sel_row = Max_abs_col_index(sel_col);
		Max = R(sel_row,sel_col);
		
		// if (abs(Max_abs) <){
		  // break;
		// }

		gamma = 1/Max;
		// printf("gamma: %26.20e \n",gamma);
		// printf("Max: %26.20e \n",Max);

		// printf("abs(Max) : %26.20e %26.20e \n",abs(Max),MinPivot);
		// printf("abs(Max) > MinPivot : %d \n ", abs(Max) > MinPivot);


	    // UPDATE U AND V
	    if (abs(Max) > MinPivot && !isnan(Max)){
	    	// cout << "OK" << endl;
			U.col(rank) = gamma*R.col(sel_col);
		    V.col(rank) = R.row(sel_row);
		}
		else{
	    	// cout << "problem" << endl;

	    	// testnan=true;

			// cout << "gamma*R.col(sel_col)" << gamma*R.col(sel_col) << endl;
			// cout << "R.row(sel_row)" << R.row(sel_row) << endl;

			U.col(rank) = VectorXd::Zero(A.rows());
			U(rank,rank) = 1;
		    V.col(rank) = VectorXd::Zero(A.cols());
		    // U(rank,rank) = MinPivot/10;
		    // V.col(rank) = VectorXd::Ones(A.cols())*MinPivot;

		}

	    // UPDATE R
	    R -= U.col(rank)*V.col(rank).transpose();
	  	rank++;

	  	
	  	Rnorm=R.norm();
	  	// printf("Rnorm,Anorm: %20.5e %20.5e \n",Rnorm,Anorm);

	  	if (isnan(Rnorm)){
	  		cout << "NaN" << endl;

			printf("abs(Max) : %26.20e %26.20e \n",abs(Max),MinPivot);

			
			cout << "A" << endl << A << endl;
	  		cout << "R" << endl << R << endl;

		  	cout << "U" << U << endl;
		  	cout << "V" << V << endl;

	  		if(A.rows() < A.cols()){
				// FAT
	      		U = MatrixXd::Identity(A.rows(),A.rows());
			      V = A.transpose();
			      rank = A.rows();
			 } 
		  	else {
		  		// THIN
			      U = A;
			      V = MatrixXd::Identity(A.cols(),A.cols());
			      rank = A.cols();
		    }
		    break;
	  	}



	  }

	  // if (testnan){
	  // 	cout << "U" << U << endl;
	  // 	cout << "V" << V << endl;
	  // }

  }
  else{

	U.conservativeResize(Eigen::NoChange,rank);
  	V.conservativeResize(Eigen::NoChange,rank);


	for(int i=0;i<rank;i++){
		for (int iCol = 0; iCol < A.cols(); iCol++){
		  Max_abs_col(iCol) = R.col(iCol).cwiseAbs().maxCoeff(&Max_abs_col_index(iCol));
		}

		int sel_row, sel_col;
		double Max_abs = Max_abs_col.maxCoeff(&sel_col);
		sel_row = Max_abs_col_index(sel_col);
		Max = R(sel_row,sel_col);
		// if (abs(Max_abs) <){
		  // break;
		// }

		gamma = 1/Max;

		// UPDATE U AND V
		U.col(i) = gamma*R.col(sel_col);
		V.col(i) = R.row(sel_row);

		// UPDATE R
		R -= U.col(i)*V.col(i).transpose();



	}
}

// printf("Rank: %d\n",rank);


  // cout << "A.rows()" << A.rows() << endl;
  // cout << "A.cols()" << A.cols() << endl;
  // printf("Rank: %d\n",rank);

  // NOT A LOW RANK MATRIX
  // if (rank*(A.rows()+A.cols()) > A.rows()*A.cols()){
 //  if(rank==maxRank){
 //  	cout << "nlr" << endl;
	// if(A.rows() < A.cols()){
	// 	// FAT
 //      U = MatrixXd::Identity(A.rows(),A.rows());
 //      V = A.transpose();
 //      rank = A.rows();
 //    } 
 //  	else {
 //  		// THIN
 //      U = A;
 //      V = MatrixXd::Identity(A.cols(),A.cols());
 //      rank = A.cols();
 //    }
 //  }
}