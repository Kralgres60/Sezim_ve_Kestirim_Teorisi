#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

/*print out the character from input file*/
//#define PRINT_CHAR_FROM_FILE

#define FILE_NAME "gene.txt"



//Nucleic Acids
const char NucleicAcids[4] = {'A','G','T','C'};
const int NumberOfNA = 4;

const char *NADef[]={
    "Ade",
    "Gua",
    "Thy",
    "Cyt",
};


double 	TrMatrix[4][4] 		= {0};
double 	NAHist[4] 	  		= {0};
double 	StateHist[4][4] 	= {0};
// stationary probability distribution after monte carlo
int 	StPrDist[4] 		= {0}; 


/* 
Nucleic Acid Codes in FASTA File
A   Ade Adenine
G   Gua Guanine
T   Thy Thymine
C   Cyt Cytosine
*/

float get_random(void)
{
    float  r_num = (float)rand()/RAND_MAX;
    return r_num;
}







void calculate_transition_matrix(const char* array,uint32_t length)
{

    /*reset Elemenets*/
    memset(NAHist,0x00,sizeof(NAHist));
    for (int i = 0; i < NumberOfNA; ++i)
    {
       for (int j = 0; j < NumberOfNA; ++j)
       {
           TrMatrix [i][j] = 0.0f;
           StateHist[i][j] = 0.0f;
       }
    }
    /*reset Elemenets*/

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
    		for (int j = 0; j < NumberOfNA; ++j)
    		{
    			if (CurrChar == NucleicAcids[j])
    			{
    				NAHist[j] += 1.0f;
    				break;
    			}
    		}
    		prevChar = CurrChar;
    		++tmpInputCounter;
    	}

    	else
    	{
    		for (int j = 0; j < NumberOfNA; ++j)
    		{
    			if (CurrChar == NucleicAcids[j])
    			{
    				NAHist[j] += 1.0f;
    				break;
    			}
    		}

    		//prevChar --> CurrChar 
    		for (int j = 0; j < NumberOfNA; ++j)
    		{
    			if (prevChar == NucleicAcids[j])
    				prevNucleoidIndex = j;
    		}

    		for (int j = 0; j < NumberOfNA; ++j)
    		{
    			if (CurrChar == NucleicAcids[j])
    			{
    				StateHist[prevNucleoidIndex][j] += 1.0f;
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
    for (int i = 0; i < NumberOfNA; ++i)
    {
        for (int j = 0; j < NumberOfNA; ++j)
        {
            if (0 == StateHist[i][j])
            {
                tmpNoTransition = true;
                tmpIndex = j;
            }

            transtionSum += StateHist[i][j];
        }

        if (transtionSum != NAHist[i])
        {
            tmpResidue = (double)(NAHist[i] - transtionSum);

            if (tmpNoTransition)
                tmpResidue /= 3.0f;
            else
                tmpResidue /= 4.0f;

            for (int k = 0; k < NumberOfNA; ++k)
            {
                if((k == tmpIndex) && (tmpNoTransition == true))
                    StateHist[i][k] = 0.0f;     
                else
                    StateHist[i][k] = StateHist[i][k] + tmpResidue; 
            }
        }

        tmpNoTransition = false;
        transtionSum = 0.0f;
    }     

    /*Transition Matrix Calculation*/
    printf("\n\r");
    printf("Transition Matrix : \n\r");

    double t_cnt = 0.0f, s_cnt = 0.0f;

    for (int i = 0; i < NumberOfNA; ++i)
    {
        for (int j = 0; j < NumberOfNA; ++j)
        {
            t_cnt           = StateHist[i][j];
            s_cnt           = NAHist[i];
            TrMatrix[i][j]  = t_cnt / s_cnt;
            printf("%.3f ", TrMatrix[i][j]);
        }
        printf("\n\r");
    }

}




int main(int argc, char** argv) 
{
	
    // Gene Squences FASTA Code from NCBI Database
    FILE *fptr = fopen(FILE_NAME, "r");

    if (fptr == NULL)
    {
    	printf("%s\n","File is not found");
    	return 1;
    }

    char 	 CurrChar 	= 0;
    char 	 prevChar 	= 0;
    uint32_t tmpInputCounter = 0;
    uint32_t prevNucleoidIndex = 0;

    while((CurrChar = fgetc(fptr)) != EOF)
    {
    	#ifdef PRINT_CHAR_FROM_FILE
    		printf("%c", CurrChar);
    	#endif

    	/*We could not do a transition for first character*/
    	if (0 == tmpInputCounter)
    	{
    		for (int i = 0; i < NumberOfNA; ++i)
    		{
    			if (CurrChar == NucleicAcids[i])
    			{
    				NAHist[i] += 1.0f;
    				break;
    			}
    		}
    		prevChar = CurrChar;
    		++tmpInputCounter;
    	}
    	else
    	{
    		for (int i = 0; i < NumberOfNA; ++i)
    		{
    			if (CurrChar == NucleicAcids[i])
    			{
    				NAHist[i] += 1.0f;
    				break;
    			}
    		}

    		//prevChar --> CurrChar 
    		for (int i = 0; i < NumberOfNA; ++i)
    		{
    			if (prevChar == NucleicAcids[i])
    				prevNucleoidIndex = i;
    		}

    		for (int i = 0; i < NumberOfNA; ++i)
    		{
    			if (CurrChar == NucleicAcids[i])
    			{
    				StateHist[prevNucleoidIndex][i] += 1.0f;
    			}
    		}

    		prevChar = CurrChar;
    	}
    }

    // Close File
    fclose(fptr);

    bool 	 tmpNoTransition = false;
    uint32_t tmpIndex = 0;
    double 	 tmpResidue = 0.0f;
    double   transtionSum = 0.0f;
    // normalization
    for (int i = 0; i < NumberOfNA; ++i)
    {
    	for (int j = 0; j < NumberOfNA; ++j)
    	{
    		if (0 == StateHist[i][j])
    		{
    			tmpNoTransition = true;
				tmpIndex = j;
    		}

    		transtionSum += StateHist[i][j];
    	}

    	if (transtionSum != NAHist[i])
    	{
    		tmpResidue = (double)(NAHist[i] - transtionSum);

    		if (tmpNoTransition)
    			tmpResidue /= 3.0f;
    		else
    			tmpResidue /= 4.0f;

    		for (int k = 0; k < NumberOfNA; ++k)
    		{
                if((k == tmpIndex) && (tmpNoTransition == true))
                    StateHist[i][k] = 0.0f;     
                else
                    StateHist[i][k] = StateHist[i][k] + tmpResidue; 
    		}
    	}

    	tmpNoTransition = false;
    	transtionSum = 0.0f;
    }

    /*Transition Matrix Calculation*/
    printf("\n\r");
    printf("Transition Matrix : \n\r");

	double t_cnt = 0.0f, s_cnt = 0.0f;

    for (int i = 0; i < NumberOfNA; ++i)
    {
    	for (int j = 0; j < NumberOfNA; ++j)
    	{
    		t_cnt 			= StateHist[i][j];
    		s_cnt 			= NAHist[i];
			TrMatrix[i][j] 	= t_cnt / s_cnt;
			printf("%.3f ", TrMatrix[i][j]);
    	}
    	printf("\n\r");
    }


   	printf("\n\r");
    printf("*****************************Sekans Sonu*******************************\n\r");


	return 0;
}
