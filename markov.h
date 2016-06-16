#ifndef MARKOV_H
#define MARKOV_H
#include <sequence.h>
#include <math.h>
#include <iostream>
#include <set>
#include <iterator>
class Result
{
public:
    Result(vector <int>* original_Intron_ids, vector <int>* result_Intron_ids);
    vector <int>* original_Intron_ids;
    vector <int>* result_Intron_ids;


};

class Markov
{
public:
    Markov();
    Markov(vector<sequence *> *data);
    vector<sequence *> *data;
    void compute_matrixes(int window);
    double matrix_plus [4][4];
    double matrix_minus [4][4];
    void add_to_matrix(double matrix[][4],char nuc1, char nuc2);
    void finalize_matrix(double matrix[][4]);
    void set_data(vector<sequence *> *data);
    void search_for_cut_placement(int window, int window2);
    bool check_if_belong_to_intron(vector <char>* seq_nuc,int start, int end);
    bool check_if_belong_to_intron2(vector <char>* seq_nuc,int start, int end);
    double prob_plus;
    double prob_minus;
    void add_to_prob(double &prob, char nuc1, char nuc2, double matrix[][4]);
    void reset_probs();
    vector <Result*> *result;
    set <vector<char>> taught_cut_beginnings;
    set <vector<char>> taught_cut_endings;


};

#endif // MARKOV_H
