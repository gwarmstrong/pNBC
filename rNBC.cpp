/* 
 * File:   rNBC.cpp
 * Author: George Armstrong
 *
 * Created on August 14, 2017
 */


#include <Rcpp.h>
#include "proposed_classification.h"
using namespace Rcpp;
using namespace std;



void quickSort(int* arr, int left, int right){
	int i = left, j = right;
	int tmp;
	int pivot = arr[(left+right)/2];

	// partition
	while(i <= j){
		while(arr[i] < pivot){
			i++;
		}
		while(arr[j] > pivot){
			j--;
		}
		if(i <= j){
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
			i++;
			j--;
		}
	}
	// recursion
	if (left < j){
		quickSort(arr, left, j);
	}
	if (i < right){
		quickSort(arr, i, right);
	}
}

int index_for_gene_id(int gene_id, int* gene_ids, int T){
	//Takes a gene id, a list of gene_ids, and the number of genes as parameters and returns corresponding index
	for (int i = 0; i < T; i++){
		if(gene_id == gene_ids[i]){
			return i;
		}
	}
	return -1;
}


int removeDuplicates(int* arr, int n){
	if(n==0 || n==1){
		return n;
	}
	int j = 0;
	for (int i = 0; i < n-1; i++){
		if(arr[i] != arr[i+1]){
			arr[j++] = arr[i];
		}
	}
	arr[j++] = arr[n-1];

	return j;
}

int cholesky(double *A, int n) { //cholesky decomposition for symmetric PD matrix A.
    int i, j, k, in;
    for (i = 0; i < n; i++) {
        in = i*n;
        for (k = 0; k < i * n; k += n)
            A[in + i] -= A[k + i] * A[k + i];
        if (A[in + i] <= 0)
            return 0; //ERROR: non-positive definite matrix!
        A[in + i] = sqrt(A[in + i]);
        for (j = i + 1; j < n; j++) {
            A[in + j] = A[in + j];
            for (k = 0; k < i * n; k += n)
                A[in + j] -= A[k + i] * A[k + j];
            A[in + j] /= A[in + i];
        }
    }
    return 1;
}

void luEvaluate(double* U, double* b, double* x, int n) {
    // Ax = b -> LUx = b. Then y is defined to be Ux  [L = U']
    int i = 0;
    int j = 0;
    // Forward solve Lx = b
    for (i = 0; i < n; i++) {
        x[i] = b[i];
        for (j = 0; j < i; j++) {
            x[i] -= U[j * n + i] * x[j];
        }
        x[i] /= U[i * n + i];
    }
    // Backward solve Ux = x
    for (i = n - 1; i >= 0; i--) {
        for (j = i + 1; j < n; j++) {
            x[i] -= U[i * n + j] * x[j];
        }
        x[i] /= U[i * n + i];
    }
}

void printMatrix(double* A, int rows, int columns){
    int i = 0;
    int j = 0;
    printf("\n\n");
    string curLine = "";
    while(i < rows){
        while(j < columns){
            if(j==columns-1){
                curLine = curLine + to_string(A[i*columns + j]);
            }
            else{
                curLine = curLine + to_string(A[i*columns + j]) + ",";
            }
            j++;
        }
        printf("%s\n\n",curLine.c_str());
        curLine = "";
        j = 0;
        i++;
    }
    printf("\n\n");
}

bool cholesky_least_squares(double*AA, double*x, double*Ab, int n) {
    //solve Ax = b problem, return x. n is the size of b.
    if (!cholesky(AA, n)) return false; //cholesky decomposition.
    luEvaluate(AA, Ab, x, n); //forward-backward substitution.
    return true;
}



void fprintMatrix(FILE* fp, int* A, int rows, int columns){
    int i = 0;
    int j = 0;
    string curLine = "";
    while(i < rows){
        while(j < columns){
            if(j==columns-1){
                curLine = curLine + to_string(A[i*columns + j]);
            }
            else{
                curLine = curLine + to_string(A[i*columns + j]) + ",";
            }
            j++;
        }
        fprintf(fp,"%s\n",curLine.c_str());
        curLine = "";
        j = 0;
        i++;
    }
    fprintf(fp,"\n");
}

double* regression(double*A, double*b, int ns, int nc){
	//compute A'A
	double* AA = new double [nc*nc];
	double* Ab = new double [nc]; 
	memset(AA,0,sizeof (double) * nc * nc);
	double* ptr_AA = 0;
	for (int d1 = 0; d1 < nc; d1++) {
    	for (int d2 = d1; d2 < nc; d2++) {
            ptr_AA = AA + d1 * nc + d2;
            for (int n = 0; n < ns * nc; n += nc){
                *ptr_AA += A[n + d1] * A[n + d2];
            }
            AA[d2 * nc + d1] = *ptr_AA; //by symmetric.
        }
    }

    for (int d = 0; d < nc; d++){
        AA[d * nc + d] += 1e-3; //regularization.
    }

	//compute Ab
	memset(Ab ,0,sizeof(double)*nc);
    for(int d1 = 0; d1<nc;d1++){
    	for(int n = 0;n < ns; n++){
            Ab[d1] += A[n*nc+d1]*b[n];
        }
    }

	//solve least quares
    double* x = new double [nc];
	cholesky_least_squares(AA, x, Ab, nc);
	delete[] AA;
	delete[] Ab;
	return x;
}

double* Multiply(double*A, double*x, int n, int m, int c){
    // A is n x m
    // x is m x c
    // then Ax is n x c
    double* Ax = new double [n*c];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < c; j++){
            Ax[c * i + j] = 0;
            for (int k = 0; k < m; k++){
                Ax[c * i + j] += A[m * i + k] * x[c * k + j];
            }
        }
    }
    return Ax;
}

double* quadratic_nosq_regression(double*A, double*b, int ns, int nc){
	//compute A with extra coefficient terms
	int rowLen = nc + ((nc - 1)*(nc)/2) + 1;//nc + ((nc + 1)*(nc)/2);
	int i;
	int j;
	double* Asq = new double[rowLen*ns];
	memset(Asq,0,sizeof(double)*(rowLen) * ns);
	for(int r = 0; r < ns; r++){
		for(int c1 = 0; c1 < nc; c1++){
			Asq[r*rowLen + c1] = A[r*nc + c1];
		}
		i = 0;
		j = 1;
		for(int c2 = nc; c2 < rowLen - 1; c2++){
			Asq[r*rowLen + c2] = A[r*nc + i]*A[r*nc + j];
			j++;
			if(j >= nc){
				i++;
				j = i + 1;
			}
		}

        Asq[r*rowLen + rowLen - 1] = 1;

	}

	double* x = regression(Asq, b, ns, rowLen);

	delete[] Asq;

	return x;
}

double* quadratic_regression(double*A, double*b, int ns, int nc){
    //compute A with extra coefficient terms
    int rowLen = nc + ((nc + 1)*(nc)/2) + 1;//nc + ((nc + 1)*(nc)/2);
    int i;
    int j;

    double* Asq = new double[rowLen*ns];
    memset(Asq,0,sizeof(double)*(rowLen) * ns);
    for(int r = 0; r < ns; r++){
        for(int c1 = 0; c1 < nc; c1++){
            Asq[r*rowLen + c1] = A[r*nc + c1];
        }
        i = 0;
        j = 0;
        for(int c2 = nc; c2 < rowLen - 1; c2++){
            Asq[r*rowLen + c2] = A[r*nc + i]*A[r*nc + j];

            j++;
            if(j >= nc){
                i++;
                j = i;//+ 1; //this line
            }
        }
        Asq[r*rowLen + rowLen - 1] = 1;
    }

    double* x = regression(Asq, b, ns, rowLen);

    delete[] Asq;

    return x;
}

double* exponential_regression(double* A, double* b, int ns, int nc){
	double* logb = new double[ns];
	for(int r = 0; r < ns; r++){
		logb[r] = log(b[r]);
	}
	double* x = regression(A, logb, ns, nc);

	delete[] logb;

	return x;
}




int gl_union(int* gene_list, int* l1, int* l2, int m){
	int i = 0, j = 0, k = 0; 
	while(i < m && j < m){
		if(l1[i] < l2[j]){
			gene_list[k++] = l1[i++];
		}
		else if(l2[j] < l1[i]){
			gene_list[k++] = l2[j++];
		}
		else{
			gene_list[k++] = l2[j++];
			i++;
		}
	}
	while(i < m){
		gene_list[k++] = l1[i++];
	}
	while(j < m){
		gene_list[k++] = l2[j++];
	}
	return removeDuplicates(gene_list,k);
}

GRAPH* graph_construction(int num_edges, int* edge1, int* edge2, int* gene_list, int T){
	// T is the number of features
	GRAPH* G = new GRAPH;
	G->T = T;
	G->vertices = new VERTEX[T];
	for (int i = 0; i <  T; i++){
		G->vertices[i].nn = 0;
		G->vertices[i].neighbors = new int [T-1]; // allocate space so that it could be connected to all nodes except itself
	}

	for (int i = 0; i < num_edges; i++){
		int gene_a = edge1[i];
		int gene_b = edge2[i];
		gene_a = index_for_gene_id(gene_a, gene_list, T);

		gene_b = index_for_gene_id(gene_b, gene_list, T);

		G->vertices[gene_a].neighbors[G->vertices[gene_a].nn] = gene_b; // add gene b to neighbors of gene a
		G->vertices[gene_a].nn++;
	 	G->vertices[gene_b].neighbors[G->vertices[gene_b].nn] = gene_a; // add gene a to neighbors of gene b
	 	G->vertices[gene_b].nn++;
	}

	return G;
}

Rcpp::NumericMatrix learn_coefficients(GRAPH* G, Rcpp::NumericMatrix B, Reg reg_type){
	int T = G->T;
	int ns = B.nrow();
	Rcpp::NumericMatrix coeffs(T + 1, T); // need an extra row for the coefficient of the constant
	double* ptr_smpl = 0; // a free point for more efficient memory access
	double*b = new double [ns];
    double* A = new double [ns * T]; //temporary memory for holding data matrix. //w = (A'A)^-1*A'*b
    double* A_strip = new double [ns * T];

    for (int j = 0; j < T; j++) { //for each vertex
        int nn = G->vertices[j].nn; //#vertices connecting to it.
        int nc = nn + 1; //number of coefficients.
        if (nc > T) {
            puts("#Coefficients cannot exceed #Attributes.");
            exit(-1);
        }
        //compute A
        for (int p = 0; p < nc; p++) { //for each vertex connecting to it. [+1 because we need extra constant input 1]
            if (p == nn) { //the 1.0 pseudo input.
                for (int n = 0; n < ns * nc; n += nc)
                    A[n + p] = 1; //pseudo input.
            } else { //take samples from neighbors.
                int id = G->vertices[j].neighbors[p]; //the id of the neighbor node.
                for (int n = 0; n < ns; n++)
                    A[n * nc + p] = B(n, id);
            }
        }
        for (int n = 0; n < ns; n++){
        	b[n] = B(n, j);
        }
        G->vertices[j].coef = regression(A, b, ns, nc);

            
 //            if (reg_type == Reg::exponential){
 //            	Gs[i]->vertices[j].coef = exponential_regression(A, b, Gs[i]->ns, nc); // changed nn to Gs[i]->ns to change from nc-1 to number of samples
 //            }
 //            else if (reg_type == Reg::quadratic_nosq){
 //                //strip 1's column from A since it messes up 
 //                for (int p = 0; p < nn; p++){
 //                    for (int n = 0; n < Gs[i]->ns; n++){
 //                        A_strip[n * nn + p] = A[n * nc + p];
 //                    }
 //                }
 //            	Gs[i]->vertices[j].coef = quadratic_nosq_regression(A_strip, b, Gs[i]->ns, nn);
 //            }
 //            else if (reg_type == Reg::quadratic){
 //                //strip 1's column from A since it messes up 
 //                for (int p = 0; p < nn; p++){
 //                    for (int n = 0; n < Gs[i]->ns; n++){
 //                        A_strip[n * nn + p] = A[n * nc + p];
 //                    }
 //                }
 //                Gs[i]->vertices[j].coef = quadratic_regression(A_strip, b, Gs[i]->ns, nn);
 //            }
 //            else {
 //            	Gs[i]->vertices[j].coef = regression(A, b, Gs[i]->ns, nc);
 //            }

        //iterate through each coefficient
        for (int p = 0; p < nc; p++){
        	if(p==nn){
        		coeffs(T,j) = G->vertices[j].coef[nn]; 
        	}
        	else{
        		int id = G->vertices[j].neighbors[p];
        		coeffs(id,j) = G->vertices[j].coef[p]; 
        		//cout << "vertex: " << j << ", neighbor: " << G->vertices[j].neighbors[p] << ", coefficient: " << G->vertices[j].coef[p] << "\n";
        	}
        }       
        
    }
    delete[] A;
    delete[] b;
    delete[] A_strip;


	return coeffs;
} 

int dropEdges(int* g1, int* g2, NumericMatrix A, int m){
	double r;
	int k = 0;
	int g1c[m];
	int g2c[m];

	// make decisions about edge inclusion
	bool* decisions = new bool[m];
	for(int i = 0; i < m; i++){
		r = ((double) rand() / (RAND_MAX));
		if(r<A(i,2)){
			decisions[i] = false;
		}
		else{
			decisions[i] = true;
		}
	}

	// make copies of the original pointers
	for (int i = 0; i<m; i++){
		g1c[i] = g1[i];
		g2c[i] = g2[i];
	}


	// rewrite original pointers to have only the included edges
	for(int i = 0; i < m; i++){
		if(decisions[i]){
			g1[k] = g1c[i];
			g2[k] = g2c[i];
			k++;
		}
	}
	delete[] decisions;
	return k;
}

//[[Rcpp::export]]
Rcpp::NumericMatrix nbcTrain(NumericMatrix A, NumericMatrix B){
	int m = A.nrow(), n = A.ncol();
  	double total = 0;

  	// find the set of genes
  	int* gl1 = new int[m];
  	int* gl2 = new int[m];
  	int* g1 = new int[m];
  	int* g2 = new int[m];
  	for (int k = 0; k < m; k++){
  		gl1[k] = (int)A(k,0);
  		gl2[k] = (int)A(k,1);
  		g1[k] = gl1[k];
  		g2[k] = gl2[k];
  	}
  	quickSort(gl1, 0, m-1);
  	quickSort(gl2, 0, m-1);
  	int* new_gene = new int[2*m];
  	int unique_genes = gl_union(new_gene, gl1, gl2, m);
  	delete[] gl1;
  	delete[] gl2;


  	int* gene_list = new int[unique_genes];
  	
  	int numEdges = dropEdges(g1, g2, A, m);

  	for (int i = 0; i < unique_genes; i++){
  		gene_list[i] = new_gene[i];
  	}
  	delete[] new_gene;

  	
  	// construct the deterministc network based on the pearson coefficients probabilities
  	GRAPH* G = graph_construction(numEdges, g1, g2, gene_list, B.ncol());

  	//learn the coefficients of the network
  	Reg reg_type = Reg::linear;
  	NumericMatrix reg_coeffs = learn_coefficients(G, B, reg_type);

  	for (int i = 0; i < m; i++){
  		for (int j = 0; j < n; j++){
  			total += A(i,j);
  		}
  	}

  	delete[] gene_list;
  
    return reg_coeffs;
  	//return List::create(Named("coefficients") =  reg_coeffs);
	//return List::create(Named("UG") = unique_genes, Named("g1") = out, Named("g2") = out2, Named("coefficients") = reg_coeffs);
}

//[[Rcpp::export]]
Rcpp::List nbcTrainN(NumericMatrix A, NumericMatrix B, int N){

    std::list<Rcpp::NumericMatrix> coeffs_mat = {};
    for(int i = 0; i < N; i++){
        coeffs_mat.push_back(nbcTrain(A,B));
    }
    Rcpp::List ans = wrap(coeffs_mat);
    return ans;
    // Rcpp::List ans(N);
    // for (int i = 0; i < N; i++){
    //     ans[i].push_back(nbcTrain(A,B));
    // }
    // return ans;
    //NumericMatrix* coeff_mats = new NumericMatrix[N];
}








