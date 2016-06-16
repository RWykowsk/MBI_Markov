#include "markov.h"

using namespace std;



Markov::Markov()
{

}

Markov::Markov(vector<sequence *> *data)
{
    this->data=data;
    for(int i=0;i<4;i++)
        for(int j=0;j<4;j++){
            matrix_minus[i][j]=0;
            matrix_plus[i][j]=0;
        }
    prob_plus=0;
    prob_minus=0;;
    result=new vector <Result*>();

}

void Markov::compute_matrixes(int window)
{
    for(int i=0;i<data->size();i++){
        vector <char>* seq_nuc=data->at(i)->seq_nuc;
        vector <int>* Intron_ids=data->at(i)->Intron_ids;
        int intron_counter=0;
        int case_number=1;
        for(int j=0;j<seq_nuc->size()-1;j++){
            switch( case_number )
            {
            case 1:
                add_to_matrix(matrix_minus,seq_nuc->at(j),seq_nuc->at(j+1));
                if(intron_counter<Intron_ids->size())
                    if(Intron_ids->at(intron_counter)==j+1){
                        case_number=2;
                        intron_counter++;
                        vector<char> beginnings;
                        for(int k=j+1;k<j+window+1;k++)
                            beginnings.push_back(seq_nuc->at(k));
                        taught_cut_beginnings.insert(beginnings);
                    }
                break;

            case 2:
                add_to_matrix(matrix_plus,seq_nuc->at(j),seq_nuc->at(j+1));
                if(Intron_ids->at(intron_counter)==j+window-1){
                    case_number=1;
                    intron_counter++;
                    vector<char> endings;
                    for(int k=j+1;k<j+window+1;k++){
                        endings.push_back(seq_nuc->at(k));
                        add_to_matrix(matrix_plus,seq_nuc->at(k-1),seq_nuc->at(k));
                    }
                    taught_cut_endings.insert(endings);
                    j=j+window;
                }
                break;
            }

        }
    }
    finalize_matrix(matrix_minus);
    finalize_matrix(matrix_plus);

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

void Markov::finalize_matrix(double matrix[][4])
{
    for(int i=0;i<4;i++){
        int sum_row=matrix[i][0]+matrix[i][1]+matrix[i][2]+matrix[i][3];
        for(int j=0;j<4;j++){
            matrix[i][j]=matrix[i][j]/sum_row;
        }
    }
}

void Markov::set_data(vector<sequence *> *data)
{
    this->data=data;
}

void Markov::search_for_cut_placement(int window, int window2)
{
    //this->window=window;
    this->result->clear();
    for(int i=0;i<data->size()-1;i++){
        vector <char>* seq_nuc=data->at(i)->seq_nuc;
        vector <int>* original_Intron_ids=data->at(i)->Intron_ids;
        vector <int>* result_Intron_ids=new vector <int>();
        vector <int> Intron_beginnings;
        vector <int> Intron_endings;
        for(int j=0;j<seq_nuc->size()-1;j++){
            //if(seq_nuc->at(j)=='G',seq_nuc->at(j+1)=='T'){
            vector <char> seq;
            if(j+window2<seq_nuc->size()){
                for(int n=j;n<j+window2;n++)
                    seq.push_back(seq_nuc->at(n));
                set<vector <char>>::iterator it;
                for (it = taught_cut_beginnings.begin(); it != taught_cut_beginnings.end(); ++it) {
                    if(seq==*it){
                        if(check_if_belong_to_intron(seq_nuc,j,j+window)){
                            Intron_beginnings.push_back(j);
                        }
                        reset_probs();
                    }
                }
                for (it = taught_cut_beginnings.begin(); it != taught_cut_beginnings.end(); ++it) {
                    if(seq==*it){
                        if(check_if_belong_to_intron2(seq_nuc,j+window2-window,j+window2)){
                            Intron_endings.push_back(j);
                        }
                        reset_probs();
                    }
                }

            }
            }



        for(int k=0;k<Intron_beginnings.size();k++){
            for(int l=0;l<Intron_endings.size();l++){
                if(Intron_endings[l]-Intron_beginnings[k]>0){
                    if(Intron_endings[l]-Intron_beginnings[k]<=2*window){
                        result_Intron_ids->push_back(Intron_beginnings[k]);
                        result_Intron_ids->push_back(Intron_endings[l]);
                    }
                else if(Intron_endings[l]-Intron_beginnings[k]>2*window){
                        break;
                    }
                }
            }
        }
        result->push_back(new Result(original_Intron_ids, result_Intron_ids));


    }
}

bool Markov::check_if_belong_to_intron(vector<char> *seq_nuc, int start, int end)
{
    if(end>seq_nuc->size())
        return false;
        for(int j=start;j<end-1;j++){
            add_to_prob(prob_plus, seq_nuc->at(j),seq_nuc->at(j+1),matrix_plus);
            add_to_prob(prob_minus, seq_nuc->at(j),seq_nuc->at(j+1),matrix_minus);
        }
    if(prob_plus<prob_minus)
        return false;
    else
        return true;
}

bool Markov::check_if_belong_to_intron2(vector<char> *seq_nuc, int start, int end)
{
    if(start<0)
        return false;

        for(int j=start;j<end-1;j++){
            add_to_prob(prob_plus, seq_nuc->at(j),seq_nuc->at(j+1),matrix_plus);
            add_to_prob(prob_minus, seq_nuc->at(j),seq_nuc->at(j+1),matrix_minus);
        }
    if(prob_plus<prob_minus)
        return false;
    else
        return true;
}

void Markov::add_to_prob(double &prob, char nuc1, char nuc2, double matrix[][4])
{
//    for(int i=0;i<4;i++){
//        for(int j=0;j<4;j++){
//            cout<<log(matrix_minus[i][j])<<'\t';
//        }
//        cout<<endl;
//    }
    if(nuc1=='A'&&nuc2=='A')
        prob=prob+log(matrix[0][0]);
    else if(nuc1=='A'&&nuc2=='C')
        prob=prob+log(matrix[0][1]);
    else if(nuc1=='A'&&nuc2=='G')
        prob=prob+log(matrix[0][2]);
    else if(nuc1=='A'&&nuc2=='T')
        prob=prob+log(matrix[0][3]);
    else if(nuc1=='C'&&nuc2=='A')
        prob=prob+log(matrix[1][0]);
    else if(nuc1=='C'&&nuc2=='C')
        prob=prob+log(matrix[1][1]);
    else if(nuc1=='C'&&nuc2=='G')
        prob=prob+log(matrix[1][2]);
    else if(nuc1=='C'&&nuc2=='T')
        prob=prob+log(matrix[1][3]);
    else if(nuc1=='G'&&nuc2=='A')
        prob=prob+log(matrix[2][0]);
    else if(nuc1=='G'&&nuc2=='C')
        prob=prob+log(matrix[2][1]);
    else if(nuc1=='G'&&nuc2=='G')
        prob=prob+log(matrix[2][2]);
    else if(nuc1=='G'&&nuc2=='T')
        prob=prob+log(matrix[2][3]);
    else if(nuc1=='T'&&nuc2=='A')
        prob=prob+log(matrix[3][0]);
    else if(nuc1=='T'&&nuc2=='C')
        prob=prob+log(matrix[3][1]);
    else if(nuc1=='T'&&nuc2=='G')
        prob=prob+log(matrix[3][2]);
    else if(nuc1=='T'&&nuc2=='T')
        prob=prob+log(matrix[3][3]);
}

void Markov::reset_probs()
{
    prob_minus=0;
    prob_plus=0;
}




Result::Result(vector<int> *original_Intron_ids, vector<int> *result_Intron_ids)
{
    this->original_Intron_ids=original_Intron_ids;
    this->result_Intron_ids=result_Intron_ids;
}
