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
    QCoreApplication a(argc, argv);


//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------


    ifstream fdat("araclean.dat");

    if (!fdat.is_open()) {
        cerr << "Can't open file!" << endl;
        return 0;
    }
    vector <sequence*>* data=new vector<sequence*>();
    int N=142;
    seq.resize(N);
    for (int i = 0; i <= N-1; ++i) {
        cout << "SEQ:  " << i << endl;
        fdat >> seq[i];
        vector <int>* introns_ids=new vector<int>(seq[i].introns);
        vector <char>* sec_nec=new vector <char>(seq[i].data) ;
        data->push_back(new sequence(sec_nec,introns_ids));
        cout << "INT:  " << seq[i].introns << endl;
//        cout << "EXO:  " << seq[i].exons << endl;
        cout << "DLEN: " << seq[i].data.size() << endl;
        cout << endl;
    }
        fdat.close();
    Markov mark(data);
    mark.compute_matrixes();
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            cout<<mark.matrix_plus[i][j]<<'\t';
        }
        cout<<endl;
    }

    cout<<endl;

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            cout<<mark.matrix_minus[i][j]<<'\t';
        }
        cout<<endl;
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

