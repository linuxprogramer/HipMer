#ifndef __DHTABLE_H
#define __DHTABLE_H

typedef struct {
    char key[MAX_READ_NAME_LEN + 1];
    int value;
} entry_t;

void dhtable_init(long tot_reads, double load_factor);
void dhtable_free(void);
int dhtable_put(char *key, int value);
int dhtable_get(char *key, int *values, int max_values);
void dhtable_check(void);

#endif
