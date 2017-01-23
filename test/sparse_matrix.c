#define _XOPEN_SOURCE 600
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <stdlib.h>

#include <pdht.h>

#define MATRIX_SIZE 8
#define BLOCK_SIZE 2
#define BLOCKS_PER_ROW (MATRIX_SIZE/BLOCK_SIZE)

extern pdht_context_t *c;
int eprintf(const char *format, ...);

int main(int argc, char **argv);

#define START_TIMER(TMR) TMR.last = pdht_get_wtime();
#define STOP_TIMER(TMR) TMR.total += pdht_get_wtime() - TMR.last;
#define READ_TIMER(TMR) TMR.total

#define KEY_A(x) (x + 1)
#define KEY_B(x) ((x + 1) << 8)
#define KEY_OUT(x) ((x + 1) << 16)


typedef struct _Block {
    float elements[BLOCK_SIZE][BLOCK_SIZE];
} Block_t;

void localhash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  (*rank).rank = 0;
  *mbits = *(unsigned long *)key;
  //*ptindex = *(unsigned long *)key % dht->nptes;
  *ptindex = 1;
}

void remotehash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  (*rank).rank = 1;
  *mbits = *(unsigned long *)key;
  //*ptindex = *(unsigned long *)key % dht->nptes;
  *ptindex = 1;
}

double rand_double() {
    return rand() / RAND_MAX;
}

void printBlock(Block_t *b) {
    for (int x = 0; x < BLOCK_SIZE; x++) {
        for (int y = 0; y < BLOCK_SIZE; y++) {
            if (b != NULL) {
                printf("%5.5f ", b->elements[x][y]);
            } else {
                printf("%5.5f ", 0.0);
            }
        }
        printf("\n");
    }
}

Block_t* filledBlock(int seed) {
    Block_t *newBlock = malloc(sizeof(Block_t));
    for (int x = 0; x < BLOCK_SIZE; x++) {
        for (int y = 0; y < BLOCK_SIZE; y++) {
            newBlock->elements[x][y] = x + y * seed;
        }
    }
    return newBlock;
}

Block_t* multiplyBlocks(Block_t *a, Block_t *b) {
    Block_t *newBlock = malloc(sizeof(Block_t));
    float sum;
    for (int x = 0; x < BLOCK_SIZE; x++) {
        for (int y = 0; y < BLOCK_SIZE; y++) {
            for (int z = 0; z < BLOCK_SIZE; z++) {
                sum += a->elements[x][z] * b->elements[z][y];
            }
        }
        newBlock->elements[x][y] = sum;
    }
    return newBlock;
}

int main(int argc, char **argv) {
    srand(0);
    pdht_t *ht;
    pdht_status_t ret;
    size_t elemsize = sizeof(Block_t);
    unsigned long key = 10;
    void *val = NULL;
    pdht_timer_t total;
    unsigned long keyA, keyB, KeyOut;
    Block_t *outBlockA, *outBlockB, *resultBlock;

    val = malloc(elemsize);
    memset(val,0,elemsize);


    // create hash table
    ht = pdht_create(sizeof(unsigned long), elemsize, PdhtModeStrict);



    pdht_sethash(ht, localhash);

    // Store rows in hash table
    for (int i=0; i < pow(BLOCKS_PER_ROW, 2); i++) {
        if (c->rank == 0) {
            //if (rand_double() > 0.5) {
                keyA = KEY_A(i);
                keyB = KEY_B(i);
                pdht_put(ht, &keyA, filledBlock(i));
                pdht_put(ht, &keyB, filledBlock(-i));
            //}
        }
    }
    
    // Print matrices for debugging
    printf("Block A:\n");
    for (int x = 0; x < BLOCKS_PER_ROW; x++) {
        for (int y = 0; y < BLOCKS_PER_ROW; y++) {
            keyA = KEY_A(x * BLOCKS_PER_ROW + y);
            pdht_get(ht, &keyA, outBlockA);
            printBlock(outBlockA);
        }
    }
    
    printf("\n\nBlock B:\n");
    for (int x = 0; x < BLOCKS_PER_ROW; x++) {
        for (int y = 0; y < BLOCKS_PER_ROW; y++) {
            keyB = KEY_B(x * BLOCKS_PER_ROW + y);
            pdht_get(ht, &keyB, outBlockB);
            printBlock(outBlockB);
        }
    }

    START_TIMER(total);
    
    // Not parallel
    for (int i = 0; i < pow(BLOCKS_PER_ROW, 2); i++) {
        if (c->rank == 0) {
            keyA = KEY_A(i);
            keyB = KEY_B(i);
            KeyOut = KEY_OUT(i);
            pdht_get(ht, &keyA, outBlockA);
            pdht_get(ht, &keyB, outBlockB);
            if (outBlockA && outBlockB) {
                resultBlock = multiplyBlocks(outBlockA, outBlockB);
                pdht_put(ht, &KeyOut, resultBlock);
            }
        }
    }
    STOP_TIMER(total);

    pdht_barrier();

    pdht_print_stats(ht);
    eprintf("%12.7f ns\n", READ_TIMER(total) * 1000000000);

done:
  pdht_free(ht);
  free(val);
}






