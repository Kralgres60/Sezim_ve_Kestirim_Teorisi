
#include "markov.h"



#define MC_ITERATION            10000



MARKOV::MARKOV() : 
m_NumberOfNA(4)
{
    m_NucleicAcids[0] = 'A';
    m_NucleicAcids[1] = 'G';
    m_NucleicAcids[2] = 'T';
    m_NucleicAcids[3] = 'C';

    memset(m_TrMatrix,0x00,sizeof(m_TrMatrix));
    memset(m_NAHist,0x00,sizeof(m_NAHist));
    memset(m_StateHist,0x00,sizeof(m_StateHist));
    memset(StPrDist,0x00,sizeof(StPrDist));
}

void 	MARKOV::reinitializeParams()
{

}

void MARKOV::estimateNucleotidsInverseTransform()
{

}


double 	MARKOV::get_random(void)
{
	double  r_num = (double)rand()/RAND_MAX;
	return 	r_num;
}















