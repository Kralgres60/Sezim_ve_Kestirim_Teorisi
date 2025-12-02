#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

/*print out the character from input file*/
//#define PRINT_TRANSITION_MATRIX

#define FILE_NAME "gene_rev.txt"



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


long int findSize(const char file_name[])
{
    // opening the file in read mode
    FILE* fp = fopen(file_name, "r");

    // checking if the file exist or not
    if (fp == NULL) {
        printf("File Not Found!\n");
        return -1;
    }

    fseek(fp, 0L, SEEK_END);

    // calculating the size of the file
    long int res = ftell(fp);

    // closing the file
    fclose(fp);

    return res;
}

/*Random Number Generator*/
float get_random(void)
{
    float  r_num = (float)rand()/RAND_MAX;
    return r_num;
}

/** 
* Inverse Transform
*/
int it_sampler(int NACode){
    
    int m, tr_value = 0;
    float u = 0.0f;
    float lc,hc;

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
    double t_cnt = 0.0f, s_cnt = 0.0f;

    for (int i = 0; i < NumberOfNA; ++i)
    {
        for (int j = 0; j < NumberOfNA; ++j)
        {
            t_cnt           = StateHist[i][j];
            s_cnt           = NAHist[i];
            TrMatrix[i][j]  = t_cnt / s_cnt;
        }
    }


#ifdef PRINT_TRANSITION_MATRIX
    printf("\n\r");
    printf("Transition Matrix : \n\r");

    for (int i = 0; i < NumberOfNA; ++i)
    {
        for (int j = 0; j < count; ++j)
        {
            printf("%.3f ", TrMatrix[i][j]);
        }
        printf("\n\r");
    }
#endif

}




int main(int argc, char** argv) 
{
	
    uint32_t fileSize = findSize(FILE_NAME); 

    if (fileSize == -1)
    {
        printf("%s\n","File is not found");
        return 1;
    }
    else
    {
        printf("File Size = %i \r\n",fileSize);
    }
    // Gene Squences FASTA Code from NCBI Database
    FILE *fptr = fopen(FILE_NAME, "r");

    if (fptr == NULL)
    {
    	printf("%s\n","File could not open");
    	return 1;
    }



    char* fileData = new char[fileSize];

    if (!fileData)
    {
        printf("insufficient Memory \r\n");
        return 1;
    }

    /*Copy Process*/
    for (int i = 0; i < fileSize; ++i)
    {
       char ch = fgetc(fptr);

       if (ch == EOF)
         break;
       

       fileData[i] = ch;
    }
    
    // Close File
    fclose(fptr);


    /*Find the X and Calculate Transition Matrix*/

    const uint32_t SliceLength = 20; 

    char prev = 0;
   
    for (int i = 0; i < fileSize; ++i)
    {
        char c = fileData[i];

        if (c == 'X')
        {
            uint16_t index = static_cast<uint16_t>(i - SliceLength);

            if (index != 0)
            {

               calculate_transition_matrix(reinterpret_cast<const char*>(&fileData[index]),SliceLength - 1);

             
               /*Find the last letter*/
               for (int j = 0; j < NumberOfNA; ++j)
               {
                   if (prev == NucleicAcids[j])
                   {
                       int m = it_sampler(j);

                       if ((0 <= m) && (m < (NumberOfNA)))
                       {
                            fileData[i] = NucleicAcids[m];
                       }   
                   }
               }

               printf("Predicted Word = %c \r\n",fileData[i]);
            }
        }

        prev = c;
    }

    /*Create a file*/
    fptr = fopen("gene_fixed.txt", "w+"); 

    /*write fixed data*/
    if (fptr != NULL)
    {
        for (int i = 0; i < fileSize; ++i)
        {
           fputc(fileData[i],fptr);
        }
    }
    else
    {
        printf("!!!! Fixed File couldn't create !!!! \r\n");
    }

    /*deallocate memory*/
    delete[] fileData;
	return 0;
}
