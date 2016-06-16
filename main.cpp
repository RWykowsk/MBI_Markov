#include <QCoreApplication>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <sequence.h>
#include <markov.h>
using namespace std;
template <typename T>
ostream& operator<<(ostream& o, const vector<T>& v) {
    copy(v.begin(), v.end(), std::ostream_iterator<T>(o, " "));

    return o;
}
//---------------------------------------------------------------------------------
class Seq {
public:
    vector<int> introns;
    vector<int> exons;
    vector<char> data;

};
istream &operator>> (istream &in, Seq &s);
//---------------------------------------------------------------------------------
vector<Seq> seq;
int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);


//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------

    int window2=8;
    int window=200;
    if (argc < 2) {
                // Tell the user how to run the program
                cerr << "Usage: " << argv[0] << " error" << endl;

                return 1;
        }
    window2=atoi(argv[1]);
    window=atoi(argv[2]);
    ifstream fdat("araclean.dat");

    if (!fdat.is_open()) {
        cerr << "Can't open file!" << endl;
        return 0;
    }
    vector <sequence*>* data_in=new vector<sequence*>();
    vector <sequence*>* data_test=new vector<sequence*>();
    int N=142;
    seq.resize(N);
    for (int i = 0; i <= N*3/4; ++i) {
        cout << "SEQ:  " << i << endl;
        fdat >> seq[i];
        vector <int>* introns_ids=new vector<int>(seq[i].introns);
        vector <char>* sec_nec=new vector <char>(seq[i].data) ;
        data_in->push_back(new sequence(sec_nec,introns_ids));
        cout << "INT:  " << seq[i].introns << endl;
//        cout << "EXO:  " << seq[i].exons << endl;
        cout << "DLEN: " << seq[i].data.size() << endl;
        cout << endl;
    }

    for (int i = N*3/4; i <= N-1; ++i) {
        cout << "SEQ:  " << i << endl;
        fdat >> seq[i];
        vector <int>* introns_ids=new vector<int>(seq[i].introns);
        vector <char>* sec_nec=new vector <char>(seq[i].data) ;
        data_test->push_back(new sequence(sec_nec,introns_ids));
        cout << "INT:  " << seq[i].introns << endl;
//        cout << "EXO:  " << seq[i].exons << endl;
        cout << "DLEN: " << seq[i].data.size() << endl;
        cout << endl;
    }

        fdat.close();
    Markov mark(data_in);
    mark.compute_matrixes(window2);
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            cout<<mark.matrix_plus[i][j]<<'\t';
        }
        cout<<endl;
    }

    cout<<endl;

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            cout<<log(mark.matrix_minus[i][j])<<'\t';
        }
        cout<<endl;
    }

    mark.set_data(data_test);
    mark.search_for_cut_placement(window,window2);
    int n=mark.result->size();
    for(int j=0;j<n;j++){
        cout<<"Wynik "<<j<<endl;
        cout<<*(mark.result->at(j)->original_Intron_ids)<<endl;
        cout<<"Rezultat"<<endl;
        cout<<*(mark.result->at(j)->result_Intron_ids)<<endl;
    }
    mark.set_data(data_in);
    mark.search_for_cut_placement(window,window2);
    n=mark.result->size();
    for(int j=0;j<n;j++){
        cout<<"Wynik "<<j<<endl;
        cout<<*(mark.result->at(j)->original_Intron_ids)<<endl;
        cout<<"Rezultat "<<j<<endl;
        cout<<*(mark.result->at(j)->result_Intron_ids)<<endl;
    }




    return 0;
}
//---------------------------------------------------------------------------------
istream & operator >> (istream & in, Seq & s)
{
    int tmp = 0;
    char sign;
    //Ignore lines of sequence header...
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //	>Seq 0 Len:
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //	2701
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //	5UTR
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //	0 884
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //	Intergenic
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //	empty-line
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //	Introns

    //read introns positions to vector
    s.introns.clear();
    while (in.peek() == ' ') {
        in >> tmp;
        s.introns.push_back(tmp);
    }
    in.clear();
    in.ignore();

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //	Exons

    //read exons
    s.exons.clear();
    while (in.peek() == ' ') {
        in >> tmp;
        s.exons.push_back(tmp);
    }
    //fix ... delete last '3' value form vector
    //s.exons.pop_back();
    in.clear();
    in.ignore();

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //	3UTR
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //	2481 2700
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //	Data


    //read data
    s.data.clear();
    while (true) {
        in >> sign;
        if (in.eof())
            break;
        if (sign == 'A' || sign == 'G' || sign == 'T' || sign == 'C') {
            s.data.push_back(sign);
            continue;
        }
        if (sign == '>') {
            break;
        }
        s.data.push_back('A'); //to fix for example 'N' sign in seq 27 at 1460 pos in data
        //cerr << ">>>>> incorrect letter fix <<<<<" << endl;
    }

    return in;
}

