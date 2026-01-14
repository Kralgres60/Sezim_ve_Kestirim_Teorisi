#include <iostream>
#include "util.h"
#include "markov.h"


#define ORIGINAL_FILE_NAME               "gene.txt"
#define MODIFIED_FILE_NAME               "gene_rev.txt"
#define MC_ITERATION                      10000
#define SLICE_NUMBER                      4

char*       g_fileData          = NULL;
char*       g_unknownNums       = NULL;
uint32_t    g_fileSize          = 0;
uint32_t    g_unknownNumSize    = 0;

const uint32_t sliceNums[SLICE_NUMBER]     = {10,20,30,40};




int main() 
{
	srand(0);


    g_fileData = parseGeneFile(MODIFIED_FILE_NAME,&g_fileSize);

    if (g_fileData == NULL)
    {
       printf("file is NULL\r\n");
       return 1;
    }

    g_unknownNumSize = calculateUnknownNumbers(MODIFIED_FILE_NAME);

    printf("File Size = %i Unknown Nums = %i \r\n",g_fileSize,g_unknownNumSize);

    if (!g_unknownNumSize)
    {
        return 1;
    }

    g_unknownNums = getUnknownNums(ORIGINAL_FILE_NAME,MODIFIED_FILE_NAME);


    MARKOV markov;

    char predicted_words[g_unknownNumSize][SLICE_NUMBER];
    markov.estimateNucleotidsInverseTransform((const char*)g_fileData,g_fileSize,(const char*)sliceNums,SLICE_NUMBER,MC_ITERATION,g_unknownNumSize,(char*)g_unknownNums);



    delete[] g_unknownNums;
    delete[] g_fileData;
	return 0;
}
