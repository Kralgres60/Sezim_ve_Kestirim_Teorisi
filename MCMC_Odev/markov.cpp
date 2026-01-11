
#include "markov.h"







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


















