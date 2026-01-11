#include "util.h"


/*Random Number Generator*/
double get_random(void)
{
    double  r_num = (double)rand()/RAND_MAX;
    return r_num;
}

uint32_t findSize(const char file_name[])
{
    // opening the file in read mode
    FILE* fp = fopen(file_name, "r");

    // checking if the file exist or not
    if (fp == NULL) {
        printf("File Not Found!\n");
        return 0;
    }

    fseek(fp, 0L, SEEK_END);

    // calculating the size of the file
    uint32_t res = ftell(fp);

    // closing the file
    fclose(fp);

    return res;
}


uint32_t    calculateUnknownNumbers(const char file_name[])
{
    uint32_t originalFileSize = findSize(file_name);

    if (originalFileSize == 0)
        return 0;

    char* fileData = new char[originalFileSize];

    if (fileData == NULL)
        return 0;

    if (!getFileData(file_name,originalFileSize,fileData))
    {
        delete[] fileData;
        return 0;
    }

    uint32_t unknownNum = 0;
    for (int i = 0; i < originalFileSize; ++i)
    {
        if (fileData[i] == 'X')
            ++unknownNum;
    }

    delete[] fileData;
    return unknownNum;
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



bool getFileData(const char file_name[],uint32_t fileSize,char*const buf)
{
    if (buf == NULL)
        return false;

    // opening the file in read mode
    FILE* fp = fopen(file_name, "r");

    if (fp == NULL)
        return false;


    for (int i = 0; i < fileSize; ++i)
    {
       char ch = fgetc(fp);

       if (ch == EOF)
         break;

        buf[i] = ch;
    }

    fclose(fp);

    return true;
}

uint32_t getRealData(char* buf,uint32_t size)
{
    if (buf == NULL)
        return 0;
    
    char* temp = new char[size];
    if (temp == NULL)
        return 0;
    

    uint32_t realSize = 0;

    for (int i = 0; i < size; ++i)
    {
        const char c = buf[i];

        if (c == 'X' || 
            c == 'A' || 
            c == 'T' || 
            c == 'G' || 
            c == 'C')
        {
            temp[realSize++];
        }
    }

    
    memcpy(buf,temp,realSize);

    delete[] temp;

    return realSize;
}

char* parseGeneFile(const char file_name[],uint32_t* size)
{
    
    uint32_t originalFileSize = findSize(file_name);

    if (originalFileSize == 0)
        return NULL;

    char* fileData = new char[originalFileSize];


    if (!getFileData(file_name,originalFileSize,fileData))
       return NULL;
    

    *size = getRealData(fileData,originalFileSize);

    if (*size == 0)
    {
        return NULL;
    }

    return fileData;
}








char*  getUnknownNums(const char original_file[], const char modified_file[])
{

    uint32_t original_file_size = findSize(original_file);
    if (original_file_size == 0)
        return NULL;

    uint32_t modified_file_size = findSize(modified_file);
    if (modified_file_size == 0)
        return NULL;


    if (original_file_size != modified_file_size)
    {
       printf("Original file and Modified file size are not same \r\n");
       return NULL;
    }

    char* original_fileData = new char[original_file_size];
    char* modified_fileData = new char[modified_file_size];

    if (!getFileData(original_file,original_file_size,original_fileData) || !getFileData(modified_file,modified_file_size,modified_fileData))
    {
       if (original_fileData != NULL)
        delete[] original_fileData;
    
       if (modified_fileData != NULL)
        delete[] modified_fileData;
       return NULL;
    }


    const uint32_t unknownNumSize = calculateUnknownNumbers(modified_file);
    char* unknownNums = new char[unknownNumSize];
    uint32_t writeCntr = 0;

    for (int i = 0; i < modified_file_size; ++i)
    {
        if (modified_fileData[i] == 'X')
        {
            if (writeCntr < unknownNumSize)
                unknownNums[writeCntr++] = original_fileData[i];
        }
    }


    delete[] original_fileData;
    delete[] modified_fileData;


    return unknownNums;
}