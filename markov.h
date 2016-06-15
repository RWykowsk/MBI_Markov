#ifndef MARKOV_H
#define MARKOV_H
#include <sequence.h>

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


};

#endif // MARKOV_H
