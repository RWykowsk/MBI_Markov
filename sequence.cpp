#include "sequence.h"



sequence::sequence(vector<char> *seq_nuc, vector<pair<int, int> > *Intron_ids)
{
    this->seq_nuc=seq_nuc;
    this->Intron_ids=Intron_ids;
}
