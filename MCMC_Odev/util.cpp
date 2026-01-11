#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <time.h>



uint32_t findSize(const char file_name[])
{
    // opening the file in read mode
    FILE* fp = fopen(file_name, "r");

    // checking if the file exist or not
    if (fp == NULL) {
        printf("File Not Found!\n");
        exit(1);
    }

    fseek(fp, 0L, SEEK_END);

    // calculating the size of the file
    uint32_t res = ftell(fp);

    // closing the file
    fclose(fp);

    return res;
}



/*Random Number Generator*/
double get_random(void)
{
    double  r_num = (double)rand()/RAND_MAX;
    return r_num;
}

/*returns real size*/
char* parseGeneFile(const char originalFile_name[],uint32_t* size)
{
    if (buf == NULL)
    {
        printf("\r\n ERROR :: buffer is NULL !!!\r\n");
        exit(1);
    }
    else
    {
    	memset(buf,0x00,size);
    }

	FILE *fptr = fopen(originalFile_name, "r");

    if (fptr == NULL)
    {
    	printf("%s\n","ERROR ::File could not open !!!");
    	exit(1);
    }

    char* fileData = new char[size];

    if (!fileData)
    {
        printf("ERROR :: insufficient Memory \r\n");
        exit(1);
    }


    /*Copy Process*/
    for (int i = 0; i < size; ++i)
    {
       char ch = fgetc(fptr);

       if (ch == EOF)
         break;
       

       fileData[i] = ch;
    }

    // Close File
    fclose(fptr);



    char *tmpBuf = new char[size];
    memset(tmpBuf,0x00,size);


    uint32_t tmpRealSize = 0;

    for (int i = 0; i < size; ++i)
    {
        const char c = fileData[i];

        if (c == 'X' || c == 'A' || c == 'T' || c == 'G' || c == 'C')
            tmpBuf[tmpRealSize++] = c;
        
    }

    //findCorrectLetter(ORIGINAL_FILE,fileParsedData,realSize);
    printf("Real File Size = %i Unknown = %i \r\n",tmpRealSize,calculateUnknownNumbers(tmpBuf,tmpRealSize));

    memcpy(buf,tmpBuf,tmpRealSize);
    delete[] tmpBuf;

    return tmpRealSize;
}

/*this function calculates the number of unknown letter numbers*/
uint32_t calculateUnknownNumbers(const char* buf,const uint32_t fileSize)
{
    if (buf == NULL)
    {
        printf("!!! buffer is NULL !!! ");
        return 0;
    }

    uint32_t tmpCnt = 0;

    for (int i = 0; i < fileSize; ++i)
    {
        if (buf[i] == 'X')
            ++tmpCnt;
    }

    return tmpCnt;
}