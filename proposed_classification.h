/* 
 * File:   common.h
 * Author: dihong
 *
 * Created on February 15, 2014, 9:36 PM
 */

#ifndef PROPOSED_CLASSIFICATION_H
#define	PROPOSED_CLASSIFICATION_H

#include <cstdlib>
//#include "svm.h"
#include <vector>
#include <cstdio>
#include<string>
#include<memory.h>
#include<cmath>
#include<algorithm>
#include<time.h>
#include<iostream>
#include <sys/time.h>
#include <unistd.h>
#include<set>
//#include"parallel.h"
#include <fstream> // SAM added this
#include <sstream> // SAM added this

enum class Reg {linear, quadratic, exponential};

typedef struct MODEL_PARAM { //This is used to describe the training/testing dataset, where dimension has been reduced.
    float* feat1; //N1xT matrix (containing reduced-dimension training features) where N1 is the number of training samples in total, T is the number of selected attributes.
    float* feat2; //N2xT matrix where N2 is the number of testing samples, containing testing samples.
    int N1, N2; //#training and #testing, respectively.
    int T; //#selected attributes.
    int* L1; //1xN1: the class label for the training samples. the label must be 1,2,...,C.
    int* L2; //1xN2: similarly to L1, for testing samples.
    float eps; //correlation threshold.
    int* conmat; //the confusion matrix for testing samples.
    int num_class; //number of classes.
    int* gene_ids; // stores the gene ids
    bool valid; //true if the GRAPH is determined [coefficients are well determined], false otherwise.
    int set_num; // training/test set number
    int fs; // feature selection algorithm
	float dpi; // dpi value (used in Aracne to remove indirect associations)
    char* suffix;
    //Reg reg_type; // type of regression to be used for classification

    MODEL_PARAM() {
        valid = true;
        feat1 = 0;
        feat2 = 0;
        L1 = 0;
        L2 = 0;
        //reg_type = Reg::quadratic; 
    }

    ~MODEL_PARAM() {
        if (feat1) delete feat1;
        if (feat2) delete feat2;
        if (L1) delete L1;
        if (L2) delete L2;
        if (conmat) delete conmat;
    }
} MODEL_PARAM;

typedef struct VERTEX { //the vertex of the graph.
    int nn; //#neighbors, i.e. vertices connecting to this vertex. 0,1,2,...
    double* coef; //the coefficients of the vertex.
    int * neighbors; //index for neighbors of the vertex.

    VERTEX() {
        coef = 0;
        neighbors = 0;
        nn = 0;
    }

    ~VERTEX() {
        if (coef) delete coef;
        if (neighbors) delete neighbors;
    }
} VERTEX;

typedef struct GRAPH { //the graph, one graph for each class.

    GRAPH() {
        vertices = 0;
        samples = 0;
        min_deg = -1;
        max_deg = -1;
        num_edge = -1;
        ns = -1;
        T = -1;
    }
    VERTEX* vertices; //1xT marix containing the the vertices of the graph, where D is the #selected attributes (equal to #vertices).
    int min_deg; //minimum degree of the vertex.
    int max_deg;
    int num_edge; //number of edges in the graph.
    int ns; //number of training samples associated with the graph
    int T; //the dimension of the features.
    int num_isolated; //#vertices isolated.
    float* samples; //ns x T matrix for selected training samples for the graph [Note that it contains training samples only for one class.]

    ~GRAPH() {
        delete samples, vertices;
    }
} GRAPH;

std::vector<GRAPH*> graph_construction(MODEL_PARAM* param);

std::vector<GRAPH*> aracne_graph_construction(MODEL_PARAM* param);

std::vector<GRAPH*> load_aracne_graphs(int train_num, int num_features, int num_classes); //Sam added this

void load_filtered_data(const std::string path, int set_num, float*& feat, int*& classes, int num_samples, int num_class, int*& gene_ids, bool train, int sel, int fs, char* suffix_c);

void load_samples_per_part(const std::string path, int*& samples_per_part, int K, char* suffix_c);

void learn_coefficients(std::vector<GRAPH*>& Gs);

int* predict(const std::vector<GRAPH*>& Gs, float* test, int nt);

bool solve_least_squares(double*A, float*x, float*b, int n);

int index_for_gene_id(int gene_id, int gene_ids[], int T);

extern struct timeval begin_ts, end_ts;  //for time elapsed measurement.

extern void tic();

extern void toc();
#endif	/* COMMON_H */
