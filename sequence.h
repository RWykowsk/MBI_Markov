#ifndef SEQUENCE_H
#define SEQUENCE_H
#include <vector>
#include <map>
using namespace std;
class sequence
{
public:
    sequence(vector <char>* seq_nuc,vector <pair<int,int>>* Intron_ids);
    vector <char>* seq_nuc;//sekwencja nukleotydow
    vector <pair<int,int>>* Intron_ids;//pary poczatku i konca intronu

};

#endif // SEQUENCE_H
