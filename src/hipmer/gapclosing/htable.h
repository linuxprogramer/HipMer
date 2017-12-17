#ifndef _HTABLE_H
#define _HTABLE_H

#include <sys/types.h>

#define NO_FREE 0
#define FREE_KEYS 1
#define FREE_VALS 2

#define NO_CHECK_DUPS 0
#define CHECK_DUPS 1

typedef struct htable *htable_t;
typedef struct htable_iter *htable_iter_t;

// returns NULL on error and sets errno
htable_t create_htable(size_t min_capacity, const char *descriptor);
void destroy_htable(htable_t h, int to_free);

// returns 0 on success, errno otherwise
int htable_put(htable_t h, char *key, void *val, int check_dups);
// returns NULL if not found
void *htable_get(htable_t h, char *key);
void *htable_del(htable_t h, char *key);
htable_iter_t htable_get_iter(htable_t h);
void *htable_get_next(htable_t h, htable_iter_t iter);

// Type safety. Always define and use these when possible
#define DEFINE_HTABLE_TYPE(value_type)                                  \
    typedef struct {char *key; value_type##_t *val;} hentry_##value_type##_t; \
    int htable_put_##value_type(htable_t h, char *key, value_type##_t *val, int check_dups) { \
        return htable_put(h, key, val, check_dups);}                    \
    value_type##_t *htable_get_##value_type(htable_t h, char *key) {    \
        return htable_get(h, key);}                                     \
    value_type##_t *htable_del_##value_type(htable_t h, char *key) {    \
        return htable_del(h, key);}                                     \
    value_type##_t *htable_get_next_##value_type(htable_t h, htable_iter_t iter) { \
        return htable_get_next(h, iter);}                                     
    
int htable_num_entries(htable_t h);
int htable_capacity(htable_t h);
long mem_htable(htable_t h);
void test_htable(void);
size_t htable_hash(char *key);

#endif /* _HASHTABLE_H_ */

