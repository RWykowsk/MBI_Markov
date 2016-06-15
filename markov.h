#ifndef MARKOV_H
#define MARKOV_H
#include <sequence.h>
#include <math.h>
class Markov
{
public:
    Markov();
    Markov(vector<sequence *> *data);
    vector<sequence *> *data;
    void compute_matrixes();
    double matrix_plus [4][4];
    double matrix_minus [4][4];
    void add_to_matrix(double matrix[][4],char nuc1, char nuc2);
    void finalize_matrix(double matrix[][4]);
    void set_data(vector<sequence *> *data);
    double prob_plus;
    double prob_minus;


};

#endif // MARKOV_H
