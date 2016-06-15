#include "file.h"


//TODO - wiekszosc jest, dodac glownie obsluge pliku
vector<sequence*> *read_from_file()
{
    vector <sequence*>* data=new vector<sequence*>();

    vector <int>* intron_row=new vector<int>();
    vector<pair<int, int> >* Intron_ids= new vector<pair<int, int> >();
    vector <char>* sec_nec_row=new vector <char>() ;
    int file_row_counter=0;
    while(1){//while (!end_of_file)

        file_row_counter=file_row_counter+8;
        intron_row=get_introns_row_from_file();
        for(int i=0;i<intron_row->size();i=i+2){
            Intron_ids->push_back(make_pair(intron_row->at(i), intron_row->at(i+1)));

        }
        file_row_counter=file_row_counter+6;
        sec_nec_row=get_sec_nec_row_from_file();
        data->push_back(new sequence(sec_nec_row,Intron_ids));
    }
    return data;



}

//TODO
vector<int> *get_introns_row_from_file()
{

}


vector<char> *get_sec_nec_row_from_file()
{

}
