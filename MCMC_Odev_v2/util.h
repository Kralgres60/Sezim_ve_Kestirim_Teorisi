#ifndef __UTIL_H__
#define __UTIL_H__


#include <stdint.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>


double 	 	get_random(void);
uint32_t 	findSize(const char file_name[]);
uint32_t 	getRealData(char* buf,uint32_t size);
bool 		getFileData(const char file_name[],uint32_t fileSize,char*const buf);
char* 		parseGeneFile(const char originalFile_name[],uint32_t* size);

uint32_t 	calculateUnknownNumbers(const char originalFile_name[]);
uint32_t 	calculateUnknownNumbers(const char* buf,const uint32_t fileSize);

char* 		getUnknownNums(const char original_file[], const char modified_file[]);




#endif




