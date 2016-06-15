#include "markov.h"



Markov::Markov()
{

}

Markov::Markov(vector<sequence *> *data)
{
    this->data=data;
}

void Markov::compute_matrixes()
{
    for(int i=0;i<data->size();i++){
        vector <char>* seq_nuc=data->at(i)->seq_nuc;
        vector <pair<int,int>>* Intron_ids=data->at(i)->Intron_ids;
        int intron_counter=0;
        int case_number=1;
        for(int j=0;j<seq_nuc->size()-1;j++){
            switch( case_number )
            {
            case 1:
                add_to_matrix(matrix_minus,seq_nuc->at(j),seq_nuc->at(j+1));
                if(Intron_ids->at(intron_counter).first==j+1)
                    case_number=2;
                break;

            case 2:
                add_to_matrix(matrix_plus,seq_nuc->at(j),seq_nuc->at(j+1));
                if(Intron_ids->at(intron_counter).second==j+1){
                    case_number=1;
                    intron_counter++;
                }
            }

        }
    }
}

void Markov::add_to_matrix(double matrix[][4], char nuc1, char nuc2)
{
    if(nuc1=='A'&&nuc2=='A')
        matrix[0][0]++;
    else if(nuc1=='A'&&nuc2=='C')
        matrix[0][1]++;
    else if(nuc1=='A'&&nuc2=='G')
        matrix[0][2]++;
    else if(nuc1=='A'&&nuc2=='T')
        matrix[0][3]++;
    else if(nuc1=='C'&&nuc2=='A')
        matrix[1][0]++;
    else if(nuc1=='C'&&nuc2=='C')
        matrix[1][1]++;
    else if(nuc1=='C'&&nuc2=='G')
        matrix[1][2]++;
    else if(nuc1=='C'&&nuc2=='T')
        matrix[1][3]++;
    else if(nuc1=='G'&&nuc2=='A')
        matrix[2][0]++;
    else if(nuc1=='G'&&nuc2=='C')
        matrix[2][1]++;
    else if(nuc1=='G'&&nuc2=='G')
        matrix[2][2]++;
    else if(nuc1=='G'&&nuc2=='T')
        matrix[2][3]++;
    else if(nuc1=='T'&&nuc2=='A')
        matrix[3][0]++;
    else if(nuc1=='T'&&nuc2=='C')
        matrix[3][1]++;
    else if(nuc1=='T'&&nuc2=='G')
        matrix[3][2]++;
    else if(nuc1=='T'&&nuc2=='T')
        matrix[3][3]++;

}


