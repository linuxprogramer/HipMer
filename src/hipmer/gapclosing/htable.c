#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <inttypes.h>

#include "utils.h"
#include "htable.h"

#if 0
#define DEBUG_MSG(fmt, ...)                                          \
    do {                                                         \
        printf("[%s:%d] " fmt, __FILE__, __LINE__, __VA_ARGS__); \
    } while (0)
#else
#define DEBUG_MSG(...)
#endif


static const unsigned int primes[] = {
    53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
    196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
    50331653, 100663319, 201326611, 402653189, 805306457, 1610612741
};

static const unsigned int num_primes = 26;// = sizeof(primes) / sizeof(primes[0]);

struct hentry {
    char *key;
    char *val;
    size_t hash;
    struct hentry *next;
};

struct htable {
    size_t capacity;
    size_t num_entries;
    size_t max_entries;
    size_t load_threshold;
    size_t collisions;
    size_t threshold_crossings;
    struct hentry **entries;
    char descriptor[255];
};

struct htable_iter {
    size_t index;
    struct hentry *next;
    int end;
};


/*
static unsigned int fnv1a_hash(const char *key, size_t len) {
    unsigned int hash = 2166136261;
    for (unsigned int i = 0; i < len; ++i)
        hash = 16777619 * (hash ^ key[i]);
    return hash ^ (hash >> 16);
}

static unsigned int larson_hash(const char *key, size_t len) {
    unsigned int hash = 0;
    for (unsigned int i = 0; i < len; ++i)
        hash = 101 * hash + key[i];
    return hash ^ (hash >> 16);
}

// http://murmurhash.googlepages.com/MurmurHash2.cpp
static unsigned int murmur2_hash(const char * key, size_t len)
{
	// 'm' and 'r' are mixing constants generated offline.
	// They're not really 'magic', they just happen to work well.
	const unsigned int m = 0x5bd1e995;
	const int r = 24;
	// Initialize the hash to a 'random' value
    unsigned int seed = 0x3FB0BB5F;
	unsigned int h = seed ^ len;
	// Mix 4 bytes at a time into the hash
	const unsigned char * data = (const unsigned char *)key;
	while(len >= 4)	{
		unsigned int k = *(unsigned int *)data;
		k *= m; 
		k ^= k >> r; 
		k *= m; 
		
		h *= m; 
		h ^= k;

		data += 4;
		len -= 4;
	}
	// Handle the last few bytes of the input array
	switch(len)	{
	case 3: h ^= data[2] << 16;
	case 2: h ^= data[1] << 8;
	case 1: h ^= data[0];
	        h *= m;
	};
	// Do a few final mixes of the hash to ensure the last few
	// bytes are well-incorporated.
	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;
	return h;
} 

// CRC-32
#define CRCPOLY 0xEDB88320
#define CRCINIT 0xFFFFFFFF

static unsigned int g_crc_precalc[4][256];

static void crc32_init(void) {
    for(unsigned int i = 0; i <= 0xFF; i++) {
        unsigned int x = i;
        for(unsigned int j = 0; j < 8; j++)
            x = (x>>1) ^ (CRCPOLY & (-(int)(x & 1)));
        g_crc_precalc[0][i] = x;
    }

    for(unsigned int i = 0; i <= 0xFF; i++) {
        unsigned int c = g_crc_precalc[0][i];
        for(unsigned int j = 1; j < 4; j++) {
            c = g_crc_precalc[0][c & 0xFF] ^ (c >> 8);
            g_crc_precalc[j][i] = c;
        }
    }
}

static unsigned int crc32_hash(const char* key, size_t len) {
    unsigned int crc = CRCINIT;
    size_t ndwords = len / sizeof(uint32_t);
    for(; ndwords; ndwords--) {
        crc ^= *(uint32_t*)key;
        crc =
            g_crc_precalc[3][(crc      ) & 0xFF] ^
            g_crc_precalc[2][(crc >>  8) & 0xFF] ^
            g_crc_precalc[1][(crc >> 16) & 0xFF] ^
            g_crc_precalc[0][(crc >> 24) & 0xFF];
        key += sizeof(uint32_t);
    }
    if (len & sizeof(uint16_t)) {
        unsigned int c = crc ^ *(uint16_t*)key;
        crc = g_crc_precalc[1][(c     ) & 0xFF] ^
            g_crc_precalc[0][(c >> 8) & 0xFF] ^ (crc >> 16);
        key += sizeof(uint16_t);
    }
    if (len & sizeof(uint8_t))
        crc = g_crc_precalc[0][(crc ^ *key) & 0xFF] ^ (crc >> 8);
    return ~crc;
}

// simple DJB2 hash
static size_t djb2_hash(char *key)
{
    size_t h = 5381;
    int c;
    while ((c = *key++))
        h = ((h<< 5) + h) + c; // h * 33 + c 
    return h;
}
*/

size_t htable_hash(char *key)
{
    //return djb2_hash(key);
    //return murmur2_hash(key, strlen(key));
    return murmur_hash2_64(key, strlen(key));
	//return MurmurHash3_x64_64(key, strlen(key));
}

#define GET_INDEX(hash) ((hash) % (h->capacity))

struct htable *create_htable(size_t min_capacity, const char *descriptor)
{
    //crc32_init();
    size_t capacity;
    int prime_found = 0;
    // get the prime just above the min capacity
    for (int i = 0; i < num_primes; i++) {
        if (primes[i] > min_capacity) { 
            capacity = primes[i]; 
            prime_found = 1;
            break; 
        }
    }
    if (!prime_found) {
        WARN("prime not found > %lu\n", min_capacity);
        errno = ENOMEM;
        return NULL;
    }
    DEBUG_MSG("Chose prime %lu for htable capacity\n", capacity);
    struct htable *h = (struct htable *)malloc(sizeof(struct htable));
    if (!h) {
        WARN("No memory for malloc of htable %s\n", descriptor);
        errno = ENOMEM;
        return NULL;
    }
    h->entries = (struct hentry **)malloc(sizeof(struct hentry *) * capacity);
    if (!h->entries) {
        WARN("No memory for malloc of htable %s entries\n", descriptor);
        free(h); 
        errno = ENOMEM;
        return NULL; 
    }
    memset(h->entries, 0, sizeof(struct hentry *) * capacity);
    h->capacity  = capacity;
    h->num_entries = 0;
    // for warnings about heavily loaded table
    h->load_threshold = capacity * 0.65;
    h->collisions = 0;
    h->threshold_crossings = 0;
    h->max_entries = 0;
    strcpy(h->descriptor, descriptor);
    return h;
}

void destroy_htable(struct htable *h, int to_free)
{
    if (h->num_entries >= h->load_threshold) {
        serial_printf("\n");
		WARN("Hash table %s load threshold %ld exceeded: entries %ld, capacity %ld\n"
             "[Results will be correct, but efficiency will be reduced because of higher collision rates]\n", 
             h->descriptor, h->load_threshold, h->num_entries, h->capacity);
    }
	for (int i = 0; i < h->capacity; i++) {
        struct hentry *entry = h->entries[i];
        while (entry) {
			struct hentry *prev_entry = entry;
			entry = entry->next;
            if (to_free & FREE_KEYS)
                free(prev_entry->key);
            if (to_free & FREE_VALS)
                free(prev_entry->val);
			free(prev_entry);
		}
	}
    free(h->entries);
    free(h);
}

int htable_num_entries(struct htable *h) 
{
    return h->num_entries;
}

int htable_capacity(htable_t h)
{
	return h->capacity;
}

int htable_put(struct htable *h, char *key, void *val, int check_dups)
{
    // make it unique by first checking for a duplicate
    size_t hash = htable_hash(key);
    size_t index = GET_INDEX(hash);
    if (check_dups == CHECK_DUPS) {
        for (struct hentry *entry = h->entries[index]; entry; entry = entry->next) {
            if (hash == entry->hash && strcmp(key, entry->key) == 0)
                return EEXIST;
        }
    }
    struct hentry *entry = (struct hentry *)malloc(sizeof(struct hentry));
    if (!entry) 
        return ENOMEM; 

    entry->key = key;
    entry->val = val;
    entry->hash = hash;
    entry->next = h->entries[index];
    if (entry->next) {
        h->collisions++;
        DEBUG_MSG("collision at %d, tot %lu\n", index, h->collisions);
    }
    h->entries[index] = entry;
    if (h->num_entries++ == h->load_threshold) {
        //WARN("Hit load threshold %lu on hash table %s\n", h->load_threshold, 
        //     h->descriptor);
        h->threshold_crossings++;
    }
    if (h->num_entries > h->max_entries)
        h->max_entries = h->num_entries;
    return 0;
}

void *htable_get(struct htable *h,  char *key)
{
    size_t hash = htable_hash(key);
    for (struct hentry *entry = h->entries[GET_INDEX(hash)]; entry; entry = entry->next) {
        if (hash == entry->hash && strcmp(key, entry->key) == 0)
            return entry->val;
    }
    return NULL;
}

void *htable_del(struct htable *h,  char *key)
{
    size_t hash = htable_hash(key);
    struct hentry **prev_link = &(h->entries[GET_INDEX(hash)]);
    for (struct hentry *entry = h->entries[GET_INDEX(hash)]; entry; entry = entry->next) {
        if (hash == entry->hash && strcmp(key, entry->key) == 0) {
            *prev_link = entry->next;
            h->num_entries--;
            void *val = entry->val;
            free(entry);
            return val;
        }
        prev_link = &(entry->next);
    }
    return NULL;
}

struct htable_iter *htable_get_iter(struct htable *h)
{
    struct htable_iter *iter = malloc(sizeof(struct htable_iter));
    if (!iter)
        return NULL;
    iter->index = 0;
    iter->next = NULL;
    return iter;
}

void *htable_get_next(struct htable *h, struct htable_iter *iter)
{
    while (iter->index < h->capacity) {
        struct hentry *entry = (iter->next ? iter->next : h->entries[iter->index]);
        if (entry) {
            void *val = entry->val;
            iter->next = entry->next;
            if (!entry->next)
                iter->index++;
            return val;
        } else {
            iter->index++;
        }
    }
    return NULL;
}

long mem_htable(struct htable *h)
{
    return sizeof(*h);
}

void test_htable(void)
{
    const int min_capacity = 500;
    const int num_elems = 1000;
    const int num_rounds = 10000;
    unsigned int seed = 31;
    
    struct htable *h = create_htable(min_capacity, "test");
    if (!h) {
        WARN("cannot create htable\n");
        return;
    }
    
    struct elem {
        char name[8];
        char data[8];
        int in_ht;
    };
    
    struct elem *elems = malloc(num_elems * sizeof(struct elem));

    for (int i = 0; i < num_elems; i++) {
        sprintf(elems[i].name, "node%d", i);
        sprintf(elems[i].data, "data%d", i);
        elems[i].in_ht = 0; 
    }
    struct elem *e;
    int puts = 0, gets = 0, dels = 0;
    for (int i = 0; i < num_rounds; i++) {
        int k = rand_r(&seed) % num_elems;
        int op = rand_r(&seed) % 3;
        if (op == 0) {
            if (!elems[k].in_ht) {
                int err = htable_put(h, elems[k].name, &elems[k], CHECK_DUPS);
                if (err)
                    WARN("Could not put %s into table\n", elems[k].name);
                elems[k].in_ht = 1;
                puts++;
            }
        } else {
            if (op == 1) {
                e = htable_get(h, elems[k].name);
                gets++;
            } else {
                e = htable_del(h, elems[k].name);
                dels++;
            }
            if (elems[k].in_ht) {
                if (!e) {
                    WARN("Could not find %s in table\n", elems[k].name);
                } else {
                    if (strcmp(elems[k].name, e->name) != 0) {
                        WARN("Found wrong elem %s should be %s\n", 
                             e->name, elems[k].name);
                    }
                }
            } else if (e) {
                WARN("Found %s in table when it should not be there\n", e->name);
            }
            if (op == 2)
                elems[k].in_ht = 0;
        }
        if (!(i % (num_rounds / 10)))
            printf("htable test round %d\n", i);
    }
    printf("Htable test concluded, capacity %lu, max elems: %lu (%.3f)\n"
           "%d puts, %d gets, %d dels, %lu collisions (%.3f), %lu crossings\n",
           h->capacity, h->max_entries, (double)h->max_entries / h->capacity,
           puts, gets, dels, h->collisions, (double)h->collisions / puts, 
           h->threshold_crossings);

    for (int i = 0; i < num_elems; i++) {
        if (elems[i].in_ht)
            htable_del(h, elems[i].name);
    }

    for (int i = 19; i >= 0; i--) {
        htable_put(h, elems[i].name, &elems[i], CHECK_DUPS);
    }

    struct htable_iter *iter = htable_get_iter(h);
    printf("Listing table via iterator: %d elems\n", htable_num_entries(h));
    int i = 0;
    while ((e = htable_get_next(h, iter)) != NULL) {
        printf("%d: name \"%s\" data \"%s\"\n", i, e->name, e->data);
        i++;
    }

}

