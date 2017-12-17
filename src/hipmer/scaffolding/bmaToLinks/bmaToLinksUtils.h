#ifndef __BMA_TO_LINKS_UTILS_H
#define __BMA_TO_LINKS_UTILS_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <upc.h>

#include "../../common/common.h"
#include "../../common/upc_compatibility.h"
#include "bma_meta.h"

#define MY_LINESIZE 1024            // linesize for bmaDataFile
#define MAX_END_SEPARATION 8000     // Maximum end separation
#define OFFSET 4000                 // default offset to accomodate negative separation values
#define MAX_LIB_NUM 30              // At most MAX_LIB_NUM libraries
#define SPLINT_UNDEFINED (-999999999)
#define SPLINT 0
#define PAIR 1

#define EXP_FACTOR 2

#define EPSILON 1e-15

#define UNDEFINED 0
#define CONVERGENT 1
#define DIVERGENT 2


int64_t __atoi64 (const char *nptr)
{
   int c;
   int64_t value;
   int sign;
   
   while (isspace((int)*nptr))
      ++nptr;
   
   c = (int)*nptr++;
   sign = c;
   if (c == '-' || c == '+')
      c = (int)*nptr++;
   
   value = 0;
   
   while (isdigit(c))
   {
      value = 10 * value + (c - '0');
      c = (int)*nptr++;
   }
   
   if (sign == '-')
      return -value;
   else
      return value;
}

/* Split a bma entry */
int split_bma_entry(char *line, char *col0, char *col1, char *col2, char *col3, char *col4, char *col5){
   char *token;
   char *aux;

   token = strtok_r(line, "\t", &aux);
   assert(token != NULL);
   strcpy(col0, token);
   token = strtok_r(NULL, "\t", &aux);
   assert(token != NULL);
   strcpy(col1, token);
   token = strtok_r(NULL, "\t", &aux);
   assert(token != NULL);
   strcpy(col2, token);
   token = strtok_r(NULL, "\t", &aux);
   assert(token != NULL);
   strcpy(col3, token);
   token = strtok_r(NULL, "\t", &aux);
   assert(token != NULL);
   strcpy(col4, token);
   token = strtok_r(NULL, "\t", &aux);
   assert(token != NULL);
   strcpy(col5, token);
   
   return 1;
}

/* Split blastmapped info existent in bma files */
int split_bmap_info(char *info, char *readStatus, char *readName, int *readStart, int *readEnd, int *readLength, char *contig_name, int *contigStart, int *contigEnd, int *contigLength, int *strand) {
   char *token;
   char *aux;
   
   token = strtok_r(info, " ", &aux);
   assert(token != NULL);
   strcpy(readStatus, token+1);
   token = strtok_r(NULL, " ", &aux);
   assert(token != NULL);
   strcpy(readName, token);
   
#ifdef ILLUMINA2
   token = strtok_r(NULL, " ", &aux);
   assert(token != NULL);
   strcat(readName, " ");
   strcat(readName, token);
#endif
   
   token = strtok_r(NULL, " ", &aux);
   assert(token != NULL);
   (*readStart) = atoi(token);
   token = strtok_r(NULL, " ", &aux);
   assert(token != NULL);
   (*readEnd) = atoi(token);
   token = strtok_r(NULL, " ", &aux);
   assert(token != NULL);
   (*readLength) = atoi(token);
   token = strtok_r(NULL, " ", &aux);
   assert(token != NULL);
   strcpy(contig_name, token);
   token = strtok_r(NULL, " ", &aux);
   assert(token != NULL);
   (*contigStart) = atoi(token);
   token = strtok_r(NULL, " ", &aux);
   assert(token != NULL);
   (*contigEnd) = atoi(token);
   token = strtok_r(NULL, " ", &aux);
   assert(token != NULL);
   (*contigLength) = atoi(token);
   token = strtok_r(NULL, " ", &aux);
   assert(token != NULL);
   (*strand) = ( memcmp(token, "Plus", 4 * sizeof(char) ) == 0  ) ? PLUS : MINUS;

   return 0;
}

/* Splits a read status */
int split_read_status(char *readIn, char *leftStat, char *rightStat, char *locStat){
   char *token;
   char *aux;
   
   token = strtok_r(readIn, ".", &aux);
   assert(token != NULL);
   strcpy(leftStat, token);
   token = strtok_r(NULL, ".", &aux);
   assert(token != NULL);
   strcpy(rightStat, token);
   token = strtok_r(NULL, ".", &aux);
   assert(token != NULL);
   strcpy(locStat, token);

   return 0;
}

/* Splits a splint */
int processSplint(char *end1, char *align1, char *end2, char *align2, bma_info *bma, int64_t *local_index, link_t *local_buffs, link_heap_t link_heap) {

   char contig1[40], contig2[40], contig_extracted[40], link[100];
   char rStat[30], readName[100];
   char leftStat[10], rightStat[10], locStat[10];
   int len1, len2, r0, r1, rL, c0, c1, cL, strand, rightCoord, leftCoord, endSeparation;
   link_t new_entry;
   int64_t hashval;
   char tail1;
   char tail2;
   char compressedTails;
   int64_t *libPairSummary = bma->libPairSummary;
   
   len1 = strlen(end1);
   memcpy(contig1, end1, (len1 - 2) * sizeof(char));
   contig1[len1 - 2] = '\0';
   len2 = strlen(end2);
   memcpy(contig2, end2, (len2 - 2) * sizeof(char));
   contig2[len2 - 2] = '\0';
   
   tail1 = *(end1+len1-1);
   tail2 = *(end2+len2-1);
   
   split_bmap_info(align1, rStat, readName, &r0, &r1, &rL, contig_extracted, &c0, &c1, &cL, &strand);
   
   /* TODO: Make sure we have set the read length for the library we are working on */
   
   if (strstr(rStat, "TRUNC") != NULL) {
      return 1;
   }
   
   split_read_status(rStat, leftStat, rightStat, locStat);
   
   if ( strcmp(rightStat, "GAP") != 0 ) {
      return 1;
   }
   
   leftCoord = ( strand == PLUS ) ? r1+cL-c1 : r1+c0-1 ;
   
   split_bmap_info(align2, rStat, readName, &r0, &r1, &rL, contig_extracted, &c0, &c1, &cL, &strand);

   if (strstr(rStat, "TRUNC") != NULL) {
      return 1;
   }
   
   split_read_status(rStat, leftStat, rightStat, locStat);
   
   if ( strcmp(leftStat, "GAP") != 0 ) {
      return 1;
   }

   rightCoord = ( strand == PLUS ) ? r0-(c0-1) : r0-(cL-c1) ;
   
   endSeparation = rightCoord - leftCoord - 1;
   
   /* In splinter we operate only on contigs */
   if (strcmp(end1, end2) < 0 ) {
      //strcpy(link, contig1);
      //strcat(link, "<=>");
      //strcat(link, contig2);
      sprintf(link, "%s<=>%s", end1, end2);
      new_entry.end1_id = atoi(contig1+6);
      new_entry.end2_id = atoi(contig2+6);
      compressedTails = compressTails(tail1, tail2);

   } else {
      //strcpy(link, contig2);
      //strcat(link, "<=>");
      //strcat(link, contig1);
      sprintf(link, "%s<=>%s", end2, end1);
      new_entry.end2_id = atoi(contig1+6);
      new_entry.end1_id = atoi(contig2+6);
      compressedTails = compressTails(tail2, tail1);

   }
   
   new_entry.link_type = SPLINT;
   new_entry.compressedTails = compressedTails;
   /* FXIME: Name conventions for libraries */
   new_entry.lib_id = bma->lib_id;
   new_entry.endSeparation = endSeparation;
   hashval = hashstr(THREADS, link);
   add_link_to_shared_heaps(&new_entry, hashval, local_index, local_buffs, link_heap);
   
   return 0;
}

/* Process a pair */
int processPair(char *end1, char *align1, char *end2, char *align2, bma_info *bma, int64_t *local_index, link_t *local_buffs, link_heap_t link_heap, int *previousMaxInsertSize, int minNetEndDistance) {
   int insertSize = bma->insertSize;
   int insertStdDev = bma->stdDev;
   int endDistance = insertSize + 3 * insertStdDev;
   char contig1[40], contig2[40], contig_extracted[40], link[100];
   char rStat[30], readName[100];
   char leftStat[10], rightStat[10], locStat[10];
   int len1, len2, r0, r1, rL, c0, c1, cL, strand, rightCoord, leftCoord, strand1, sStart1, sStop1, obj1_len, obj2_len, d1, d2, leftend1, leftend2, rightend1, rightend2, strand2, sStart2, sStop2, orientation, separation, maxSepZ, endSeparation;
   double sepZ;
   link_t new_entry;
   int64_t hashval;
   int objectNameOffset;
   char tail1;
   char tail2;
   char compressedTails, rCompressedTails;
   int64_t *libPairSummary = bma->libPairSummary;
   
   int innieMaxSep = 1000;
   int innieLib = bma->innieRemoval;
   
   len1 = strlen(end1);
   memcpy(contig1, end1, (len1 - 2) * sizeof(char));
   contig1[len1 - 2] = '\0';
   len2 = strlen(end2);
   memcpy(contig2, end2, (len2 - 2) * sizeof(char));
   contig2[len2 - 2] = '\0';
   
   tail1 = *(end1+len1-1);
   tail2 = *(end2+len2-1);
   compressedTails = compressTails(tail1, tail2);
   rCompressedTails = compressTails(tail2, tail1);
   
   split_bmap_info(align1, rStat, readName, &r0, &r1, &rL, contig_extracted, &c0, &c1, &cL, &strand);
   /* readLengths can vary, always reset it if it is smaller or much larger */
   assert(rL > 0);
   bma->readLength = rL;
   if (bma->readLength < rL || bma->readLength > rL + 50)
      bma->readLength = rL;
   
   /* TODO: Add objectLengths relative instructions */
   strand1 = strand;
   sStart1 = c0;
   sStop1 = c1;
   obj1_len = cL;
   
   if (strstr(rStat, "TRUNC") != NULL) {
      libPairSummary[TRUNCA] += 1;
      return 1;
   }
   
   if (cL < (*previousMaxInsertSize)/2) {
      libPairSummary[SMALL] += 1;
      return 1;
   }
   
   d1 = -1;
   leftend1 = c0;
   rightend1 = c1;
   if (strand == PLUS) {
      d1 = cL-c0+1+r0-1;
      leftend1 -= r0-1;
      rightend1 += rL - r1;
   } else {
      d1 = c1+r0-1;
      leftend1 -= rL-r1;
      rightend1 += r0 -1;
   }
   
   split_bmap_info(align2, rStat, readName, &r0, &r1, &rL, contig_extracted, &c0, &c1, &cL, &strand);
   // read lenghts can vary, choose the larger of the two paired reads
   if (bma->readLength < rL)
      bma->readLength = rL;
   strand2 = strand;
   sStart2 = c0;
   sStop2 = c1;
   obj2_len = cL;

   if (strstr(rStat, "TRUNC") != NULL) {
      libPairSummary[TRUNCA] += 1;
      return 1;
   }
   
   if (cL < (*previousMaxInsertSize)/2) {
      libPairSummary[SMALL] += 1;
      return 1;
   }
   
   d2 = -1;
   leftend2 = c0;
   rightend2 = c1;
   if (strand == PLUS) {
      d2 = cL-c0+1+r0-1;
      leftend2 -= r0-1;
      rightend2 += rL - r1;
   } else {
      d2 = c1+r0-1;
      leftend2 -= rL-r1;
      rightend2 += r0 -1;
   }
   
   /* Check for inappropriate self-linkage */
   if (strcmp(end1, end2) == 0) {
      orientation = UNDEFINED;
      separation = 0;
      sepZ = 0.0;
      maxSepZ = 5;
      
      if (strand1 != strand2) {
         if ( ((strand1 == PLUS) && (sStart1 <= sStart2) && (sStop1 <= sStop2)) || ((strand2 == PLUS && (sStart2 <= sStart1) && (sStop2 <= sStop1)) )) {
            orientation = CONVERGENT;
            separation = (strand1 == PLUS) ? (rightend2 - leftend1 + 1) : (rightend1 - leftend2 + 1);
            sepZ = fabs((separation-insertSize)/(1.0*insertStdDev));
         } else if ( ((strand1 == PLUS) && (sStart1 > sStart2)) || ((strand2 == PLUS) && (sStart2 > sStart1)) ) {
            orientation = DIVERGENT;
            separation = (strand1 == PLUS) ? (rightend1 - leftend2 + 1) : (rightend2 - leftend1 + 1);
         }
      }
      
      if ( (orientation == CONVERGENT) && (sepZ < maxSepZ) ) {
         libPairSummary[INORM] += 1;
         return 1;
      } else if ( (orientation == CONVERGENT) && (separation < minNetEndDistance) ) {
         libPairSummary[ISHIRT] += 1;
         return 1;
      } else {
         if (innieLib && (orientation == DIVERGENT) && (separation < innieMaxSep) ){
            libPairSummary[INNIE] += 1;
            return 1;
         }
         /* Report self-link */
         libPairSummary[SELFL] += 1;
         //printf("SELF-LINK\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n", end1, strand1, sStart1, sStop1, end2, strand2, sStart2, sStop2, orientation, sepZ);
      }
   }
   
   if ( d1 == -1 || d2 == -1 ) {
      return 1;
   }
   
   if ( (d1 >= endDistance) || (d2 >= endDistance) || ((d1+d2) <= minNetEndDistance)) {
      libPairSummary[EDIST] += 1;
      return 1;
   }
   
   endSeparation = insertSize - (d1+d2);
   objectNameOffset = ((*contig1) == 'C') ? 6 : 8 ;
   scaffoldSpans = ((*contig1) == 'C') ? 0 : 1;
   
   if (strcmp(end1, end2) < 0 ) {
      //strcpy(link, contig1);
      //strcat(link, "<=>");
      //strcat(link, contig2);
      sprintf(link, "%s<=>%s", end1, end2);
      new_entry.end1_id = atoi(contig1+objectNameOffset);
      new_entry.end2_id = atoi(contig2+objectNameOffset);
      //sprintf(link, "%d.%c<=>%d.%c", new_entry.end1_id, tail1 , new_entry.end2_id, tail2);
      //printf("Example link %s\n", link);
      new_entry.d1 = d1;
      new_entry.d2 = d2;
      new_entry.o1_length = obj1_len;
      new_entry.o2_length = obj2_len;
      new_entry.nature = (*contig1);
      new_entry.compressedTails = compressedTails;
   } else {
      //strcpy(link, contig2);
      //strcat(link, "<=>");
      //strcat(link, contig1);
      sprintf(link, "%s<=>%s", end2, end1);
      new_entry.end2_id = atoi(contig1+objectNameOffset);
      new_entry.end1_id = atoi(contig2+objectNameOffset);
      //sprintf(link, "%d.%c<=>%d.%c", new_entry.end1_id, tail2 , new_entry.end2_id, tail1);
      new_entry.d1 = d2;
      new_entry.d2 = d1;
      new_entry.o1_length = obj2_len;
      new_entry.o2_length = obj1_len;
      new_entry.nature = (*contig1);
      new_entry.compressedTails = rCompressedTails;
   }
   
   new_entry.link_type = PAIR;
   new_entry.lib_id = bma->lib_id;
   new_entry.endSeparation = endSeparation;
   hashval = hashstr(THREADS, link);
   add_link_to_shared_heaps(&new_entry, hashval, local_index, local_buffs, link_heap);
   
   return 0;
}

#endif
