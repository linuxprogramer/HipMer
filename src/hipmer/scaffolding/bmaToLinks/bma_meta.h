#ifndef BMA_H
#define BMA_H

#include <assert.h>

/* Pair Categories */
#define REASONS 9
#define TRUNCA 0
#define SMALL 1
#define INORM 2
#define ISHIRT 3
#define INNIE 4
#define SELFL 5
#define EDIST 6
#define REDUN 7
#define ACCPT 8
#define MY_LINESIZE 1024             // linesize for bmaDataFile

typedef struct bma_info {
   char bmaFile[MAX_FILE_PATH];
   int lib_id;
   int insertSize;
   int stdDev;
   int readLength;
   int64_t libPairSummary[REASONS];
   int innieRemoval;
} bma_info;


/* Split a bma line */
int split_bmaLine(char *line, struct bma_info *bma, int lib_id){
   char *token;
   char *aux;
   int len;

   token = strtok_r(line, "\t", &aux);
   assert(token != NULL);
   assert(strlen(token) > 0);
   bma->insertSize = atoi(token);

   token = strtok_r(NULL, "\t", &aux);
   assert(token != NULL);
   assert(strlen(token) > 0);
   bma->stdDev = atoi(token);

   token = strtok_r(NULL, "\t", &aux);
   assert(token != NULL);
   assert(strlen(token) > 0);
   bma->innieRemoval = atoi(token);

   token = strtok_r(NULL, "\t", &aux);
   assert(token != NULL);
   assert(strlen(token) > 0);
   bma->readLength = atoi(token);

   token = strtok_r(NULL, "\t", &aux);
   assert(token != NULL);
   assert(strlen(token) > 0);
   len = strlen(token);
   strncpy(bma->bmaFile, token, len-1);
   bma->bmaFile[len-1] = '\0';
   bma->lib_id = lib_id;
   return 1;
}


#endif
