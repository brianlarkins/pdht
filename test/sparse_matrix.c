#define _XOPEN_SOURCE 600
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <stdlib.h>

#include <pdht.h>

// Options
#define FILL_RATE 1.0  // Pecentage of blocks to be filled with data
#define MATRIX_SIZE 8 // total number of elements
#define BLOCK_SIZE 2 // Number of elements in each row/col in block elements

// Helpful macros
#define BLOCKS_PER_ROW (MATRIX_SIZE/BLOCK_SIZE)
#define TOTAL_BLOCKS (BLOCKS_PER_ROW * BLOCKS_PER_ROW)

extern pdht_context_t *c;
int eprintf(const char *format, ...);

int main(int argc, char **argv);

pdht_t *ht;

#define START_TIMER(TMR) TMR.last = pdht_get_wtime();
#define STOP_TIMER(TMR) TMR.total += pdht_get_wtime() - TMR.last;
#define READ_TIMER(TMR) TMR.total

// Debug log
#if 0
#define DEBUG_LOG(...) do{ printf(  __VA_ARGS__ ); } while( 0 )
#else
#define DEBUG_LOG(...)
#endif

typedef struct _Block {
    float elements[BLOCK_SIZE][BLOCK_SIZE];
} Block_t;


// No portals test
//#define test
#ifdef test

#define pdht_get test_get
#define pdht_put test_put

Block_t *block_map[TOTAL_BLOCKS * 3] = {0};

pdht_status_t test_get(pdht_t *dht, void *key, void *value) {
    uint64_t k = *(uint64_t *)key;
    
    Block_t *outBlock = block_map[k];
    if (outBlock == 0) {
        return PdhtStatusNotFound;
    }
    
    *(Block_t *)value = *outBlock;
    return PdhtStatusOK;
}

pdht_status_t test_put(pdht_t *dht, void *key, void *value) {
    uint64_t k = *(uint64_t *)key;
    
    block_map[k] = (Block_t *)value;
    return PdhtStatusOK;
}

int KEY_A(int x) {
    return x;
}

int KEY_B(int x) {
    return x + TOTAL_BLOCKS;
}

int KEY_OUT(int x) {
    return x + (TOTAL_BLOCKS * 2);;
}

int KEY2_A(int x, int y) {
    return KEY_A(y * BLOCKS_PER_ROW + x);
}

int KEY2_B(int x, int y) {
    return KEY_B(y * BLOCKS_PER_ROW + x);
}

int KEY2_OUT(int x, int y) {
    return KEY_OUT(y * BLOCKS_PER_ROW + x);
}

#else

int KEY_A(int x) {
    return (x + 1);
}

int KEY_B(int x) {
    return (x + 1) << 8;
}

int KEY_OUT(int x) {
    return (x + 1) << 16;
}

int KEY2_A(int x, int y) {
    return KEY_A(y * BLOCKS_PER_ROW + x);
}

int KEY2_B(int x, int y) {
    return KEY_B(y * BLOCKS_PER_ROW + x);
}

int KEY2_OUT(int x, int y) {
    return KEY_OUT(y * BLOCKS_PER_ROW + x);
}

#endif

void localhash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
    (*rank).rank = 0;
    *mbits = *(unsigned long *)key;
    //*ptindex = *(unsigned long *)key % dht->nptes;
    *ptindex = 0;
}

void remotehash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
    (*rank).rank = 1;
    *mbits = *(unsigned long *)key;
    //*ptindex = *(unsigned long *)key % dht->nptes;
    *ptindex = 1;
}

static double rand_double() {
    return rand() / (double)RAND_MAX;
}


// Helper functions

Block_t* filledBlock(int seed) {
    Block_t *newBlock = malloc(sizeof(Block_t));
    for (int x = 0; x < BLOCK_SIZE; x++) {
        for (int y = 0; y < BLOCK_SIZE; y++) {
            newBlock->elements[x][y] = x * 3 + y + 10 * seed + 1;
        }
    }
    return newBlock;
}

Block_t* multiplyBlocks(int row, int col) {
    DEBUG_LOG("\nCalculating block: (%d, %d)\n", row, col);
    unsigned long keyA, keyB;
    Block_t *newBlock = calloc(sizeof(Block_t), 1);
    Block_t blockA;
    Block_t blockB;
    // Loop over blocks to multiply
    for (int i = 0; i < BLOCKS_PER_ROW; i++) {
        
        
        keyA = KEY2_A(row, i);
        keyB = KEY2_B(i, col);

        if (pdht_get(ht, &keyA, &blockA) == PdhtStatusOK && pdht_get(ht, &keyB, &blockB) == PdhtStatusOK) {
            DEBUG_LOG("Multiplying blocks: (%d, %d) and (%d, %d)\n", i, row, col, i);
            // Loop over each element in the blocks
            for (int x = 0; x < BLOCK_SIZE; x++) {
                for (int y = 0; y < BLOCK_SIZE; y++) {
                    DEBUG_LOG("Element (%d, %d): %f\n", x, y, newBlock->elements[x][y]);
                    for (int z = 0; z < BLOCK_SIZE; z++) {
                        DEBUG_LOG("%f * %f\n", blockA.elements[z][y], blockB.elements[x][z]);
                        newBlock->elements[x][y] += blockA.elements[z][y] * blockB.elements[x][z];
                    }
                    DEBUG_LOG("\n");
                }
            }
        } else {
            DEBUG_LOG("Skipping blocks: (%d, %d) and (%d, %d)\n", i, row, col, i);
            DEBUG_LOG("Reason: %d %d\n", pdht_get(ht, &keyA, &blockA) == PdhtStatusOK, PdhtStatusOK && pdht_get(ht, &keyB, &blockB) == PdhtStatusOK);
        }
    }
    
    return newBlock;
}

// pass in the 2 value keying function of the matrix to print
void printMatrix(int (*f)(int, int)) {
    unsigned long key;
    Block_t outBlock;
    for (int row = 0; row < BLOCKS_PER_ROW; row++) {
        for (int y = 0; y < BLOCK_SIZE; y++) {
            for (int col = 0; col < BLOCKS_PER_ROW; col++) {
                key = f(row, col);
                pdht_status_t status = pdht_get(ht, &key, &outBlock);
                for (int x = 0; x < BLOCK_SIZE; x++) {
                    if (status == PdhtStatusOK) {
                        printf("%6.1f ", outBlock.elements[x][y]);
                    } else {
                        //printf("%6.1f ", 0.0);
                        printf("   x.x ");
                    }
                }
                printf(" |");
            }
            printf("\n");
        }
        printf("\n");
    }
}


int main(int argc, char **argv) {
    printf("Starting\n");
    srand(0);
    pdht_status_t ret;
    size_t elemsize = sizeof(Block_t);
    unsigned long key = 10;
    void *val = NULL;
    pdht_timer_t total;
    unsigned long keyA, keyB, KeyOut;
    Block_t outBlockA;
    Block_t outBlockB;
    Block_t *resultBlock;
    
    val = malloc(elemsize);
    memset(val,0,elemsize);
    
    
    //setenv("PTL_IGNORE_UMMUNOTIFY", "1", 1);
    setenv("PTL_PROGRESS_NOSLEEP", "1", 1);
    
    // setup experimental configuration
    pdht_config_t cfg;
    cfg.nptes        = 1;
    cfg.pendmode     = PdhtPendingPoll;
    //cfg.pendmode     = PdhtPendingTriggered;
    cfg.maxentries   = 250000;
    cfg.pendq_size   = 100000;
    cfg.ptalloc_opts = 0;
    //cfg.ptalloc_opts = PTL_PT_MATCH_UNORDERED;
    pdht_tune(PDHT_TUNE_ALL, &cfg);
    
    // create hash table
    ht = pdht_create(sizeof(Block_t), elemsize, PdhtModeStrict);
    
    
    pdht_sethash(ht, localhash);
    
    printf("Filling blocks\n");
    
    // Store rows in hash table
    for (int row = 0; row < BLOCKS_PER_ROW; row++) {
        for (int col = 0; col < BLOCKS_PER_ROW; col++) {
            if (c->rank == 0) {
                if (rand_double() < FILL_RATE) {
                    keyA = KEY2_A(row, col);
                    pdht_put(ht, &keyA, filledBlock(row * BLOCKS_PER_ROW + col));
                }
                if (rand_double() < FILL_RATE) {
                    keyB = KEY2_B(row, col);
                    pdht_put(ht, &keyB, filledBlock(-((row * BLOCKS_PER_ROW + col) + 1)));
                }
            }
        }
    }
    
    printf("Matrix A\n");
    printMatrix(&KEY2_A);
    
    printf("Matrix B\n");
    printMatrix(&KEY2_B);
    
    pdht_barrier();
    
    PDHT_START_ATIMER(total);
    
    // Not parallel
    for (int row = 0; row < BLOCKS_PER_ROW; row++) {
        for (int col = 0; col < BLOCKS_PER_ROW; col++) {
            KeyOut = KEY2_OUT(row, col);
            resultBlock = multiplyBlocks(row, col);
            pdht_put(ht, &KeyOut, resultBlock);
        }
    }
    
    pdht_barrier();
    
    printf("Result Matrix\n");
    printMatrix(&KEY2_OUT);
    
    PDHT_STOP_ATIMER(total);
    
    pdht_print_stats(ht);
    
    eprintf("total elapsed time: %12.7f ns\n", (double)PDHT_READ_ATIMER(total));
    
done:
    pdht_free(ht);
    free(val);
}






