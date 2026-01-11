#ifndef __UTIL_H__
#define __UTIL_H__


#include <stdint.h>


double 	 	get_random(void);
uint32_t 	findSize(const char file_name[]);

char* 		parseGeneFile(const char originalFile_name[],uint32_t* size);

uint32_t 	calculateUnknownNumbers(const char* buf,const uint32_t fileSize);




void 		findCorrectLetter(const char originalFile_name[],const char* buf,const uint32_t tmpfileSize)












#endif




