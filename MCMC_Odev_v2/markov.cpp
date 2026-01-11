
#include "markov.h"
#include "util.h"






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
    memset(m_NAHist,0x00,sizeof(m_NAHist));
    for (int i = 0; i < m_NumberOfNA; ++i)
    {
       for (int j = 0; j < m_NumberOfNA; ++j)
       {
           m_TrMatrix [i][j] = 0.0f;
           m_StateHist[i][j] = 0.0f;
       }
    }
}

void    MARKOV::calculate_transition_matrix(const char* array,uint32_t length)
{
    reinitializeParams();

    char     CurrChar   = 0;
    char     prevChar   = 0;
    uint32_t tmpInputCounter = 0;
    uint32_t prevNucleoidIndex = 0;


    for (int i = 0; i < length; ++i)
    {
        CurrChar = array[i];

        /*We could not do a transition for first character*/
        if (0 == tmpInputCounter)
        {
            for (int j = 0; j < m_NumberOfNA; ++j)
            {
                if (CurrChar == m_NucleicAcids[j])
                {
                    m_NAHist[j] += 1.0f;
                    break;
                }
            }
            prevChar = CurrChar;
            ++tmpInputCounter;
        }

        else
        {
            for (int j = 0; j < m_NumberOfNA; ++j)
            {
                if (CurrChar == m_NucleicAcids[j])
                {
                    m_NAHist[j] += 1.0f;
                    break;
                }
            }

            //prevChar --> CurrChar 
            for (int j = 0; j < m_NumberOfNA; ++j)
            {
                if (prevChar == m_NucleicAcids[j])
                    prevNucleoidIndex = j;
            }

            for (int j = 0; j < m_NumberOfNA; ++j)
            {
                if (CurrChar == m_NucleicAcids[j])
                {
                    m_StateHist[prevNucleoidIndex][j] += 1.0f;
                }
            }

            prevChar = CurrChar;
        }
    }
    bool     tmpNoTransition = false;
    uint32_t tmpIndex = 0;
    double   tmpResidue = 0.0f;
    double   transtionSum = 0.0f;

    // normalization
    for (int i = 0; i < m_NumberOfNA; ++i)
    {
        for (int j = 0; j < m_NumberOfNA; ++j)
        {
            if (0 == m_StateHist[i][j])
            {
                tmpNoTransition = true;
                tmpIndex = j;
            }

            transtionSum += m_StateHist[i][j];
        }

        if (transtionSum != m_NAHist[i])
        {
            tmpResidue = (double)(m_NAHist[i] - transtionSum);

            if (tmpNoTransition)
                tmpResidue /= 3.0f;
            else
                tmpResidue /= 4.0f;

            for (int k = 0; k < m_NumberOfNA; ++k)
            {
                if((k == tmpIndex) && (tmpNoTransition == true))
                    m_StateHist[i][k] = 0.0f;     
                else
                    m_StateHist[i][k] = m_StateHist[i][k] + tmpResidue; 
            }
        }

        tmpNoTransition = false;
        transtionSum = 0.0f;
    }     


    /*Transition Matrix Calculation*/
    double t_cnt = 0.0f, s_cnt = 0.0f;

    for (int i = 0; i < m_NumberOfNA; ++i)
    {
        for (int j = 0; j < m_NumberOfNA; ++j)
        {
            t_cnt           = m_StateHist[i][j];
            s_cnt           = m_NAHist[i];

            if (s_cnt == 0)
                m_TrMatrix[i][j]  = 0.0f;
            else
                m_TrMatrix[i][j]  = static_cast<double>(t_cnt / s_cnt);
        }
    }


#ifdef PRINT_TRANSITION_MATRIX
    printf("\n\r");
    printf("Transition Matrix : \n\r");

    for (int i = 0; i < m_NumberOfNA; ++i)
    {
        for (int j = 0; j < count; ++j)
        {
            printf("%.3f ", m_TrMatrix[i][j]);
        }
        printf("\n\r");
    }
#endif

}
int     MARKOV::it_sampler(int NACode,double TrMatrix[4][4])
{
    int m, tr_value = 0;
    double u = 0.0f;
    double lc,hc; 

    u = get_random(); //0-1 arasi rastgele sayi
    
    if(u <= TrMatrix[NACode][0]){ //A
        tr_value = 0;   
    }
    
    lc = TrMatrix[NACode][0];
    hc = TrMatrix[NACode][0] + TrMatrix[NACode][1];
    if(lc<=u && u<hc){ //G
        tr_value = 1;   
    }
    
    lc = TrMatrix[NACode][0] + TrMatrix[NACode][1];
    hc =TrMatrix[NACode][0] + TrMatrix[NACode][1] + TrMatrix[NACode][2];
    if(lc<=u && u<hc){ // T
        tr_value = 2;   
    }
    
    lc = TrMatrix[NACode][0] + TrMatrix[NACode][1] + TrMatrix[NACode][2];
    hc = TrMatrix[NACode][0] + TrMatrix[NACode][1] + TrMatrix[NACode][2] + TrMatrix[NACode][3];
    if(lc<=u && u<hc){ // C
        tr_value = 3;   
    }
    
    return tr_value;
}
void MARKOV::estimateNucleotidsInverseTransform(const char* fileData,
                                                const uint32_t fileDataSize,
                                                const char* silece,
                                                const uint32_t sileceNumber,
                                                const uint32_t mcIteration,
                                                const uint32_t unknownWordNumbers)
{
    int         row    = 0;
    int         column = 0;
    char        prev = 0;
    
    char        predictedWords[unknownWordNumbers][sileceNumber];
    memset(predictedWords,0x00,sizeof(predictedWords));

    for (int i = 0; i < fileDataSize; ++i)
    {
        char c = fileData[i];

        if (c == 'X')
        {
            column = 0;

            for (int s = 0; s < sileceNumber; ++s)
            {
                uint32_t silence  = silece[s];
                if (i >= silence)
                {
                    uint32_t index = i - silence;
                    // Transition matrix hesapla
                    calculate_transition_matrix(&fileData[index], silence);


                    /*Find the last letter*/

                    for (int j = 0; j < m_NumberOfNA; ++j)
                    {
                        if (prev == m_NucleicAcids[j])
                        {
                            int mc_count[4] = {0};

                            for (int iter  = 0; iter  < mcIteration; ++iter )
                            {
                                int m = it_sampler(j,m_TrMatrix);
                                mc_count[m]++;
                            }

                            int best = 0;

                            for (int n = 0; n < 4; ++n)
                            {
                                if (mc_count[n] > mc_count[best])
                                    best = n;
                            }

                            char predicted = m_NucleicAcids[best];

                            predictedWords[row][column++] = predicted;
                        }
                    }
                }
            }
            row++;
        }

        prev = c;
    }

}


















