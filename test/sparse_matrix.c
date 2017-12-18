#define _XOPEN_SOURCE 600
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <stdlib.h>

#include <pdht.h>
#include <city.h>

// Options
#define FILL_RATE .5  // Pecentage of blocks to be filled with data
#define MATRIX_SIZE 1000 // total number of elements in a row
#define BLOCK_SIZE 2 // Number of elements in each row/col in blocks

// Helpful Constants
#define BLOCKS_PER_ROW (MATRIX_SIZE/BLOCK_SIZE)
#define TOTAL_BLOCKS (BLOCKS_PER_ROW * BLOCKS_PER_ROW)

extern pdht_context_t *c;
int eprintf(const char *format, ...);

int main(int argc, char **argv);



// Debug log
#if 0
#define DEBUG_LOG(...) do{ printf(  __VA_ARGS__ ); } while( 0 )
#else
#define DEBUG_LOG(...)
#endif

//#define PRINT_MATRICES

typedef struct _Block {
    float elements[BLOCK_SIZE][BLOCK_SIZE];
} Block_t;


// No portals test
//#define NO_PORTALS
//#ifdef NO_PORTALS

//#define pdht_get test_get
//#define pdht_put test_put

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

//#endif

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


void localhash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
    (*rank).rank = 0;
    *mbits = CityHash64((char *)key,dht->keysize);
    //*mbits = *(unsigned long *)key;
    *ptindex = *mbits % dht->ptl.nptes;
    //*ptindex = 0;
}

static double rand_double() {
    return rand() / (double)RAND_MAX;
}


// Helper functions

Block_t* filledBlock(int seed,Block_t *newBlock) {
    //Block_t *newBlock = malloc(sizeof(Block_t));
    for (int x = 0; x < BLOCK_SIZE; x++) {
        for (int y = 0; y < BLOCK_SIZE; y++) {
            newBlock->elements[x][y] = x * 3 + y + 10 * seed + 1;
        }
    }
    return newBlock;
}

void multiplyBlocks(int row, int col,pdht_t *ht,Block_t *newBlock) {
    DEBUG_LOG("\nCalculating block: (%d, %d)\n", row, col);
    unsigned long keyA, keyB;
    //Block_t newBlock;
    Block_t blockA;
    Block_t blockB;
    memset(newBlock,0,sizeof(Block_t));
    //resultBlock = calloc(sizeof(Block_t),1);
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
    
    //return newBlock;
}

// pass in the 2 value keying function of the matrix to print
void printMatrix(int (*f)(int, int),pdht_t *ht) {
    unsigned long key;
    Block_t outBlock;
    pdht_status_t status;
    for (int row = 0; row < BLOCKS_PER_ROW; row++) {
        for (int y = 0; y < BLOCK_SIZE; y++) {
            for (int col = 0; col < BLOCKS_PER_ROW; col++) {
                key = f(row, col);
                status = pdht_get(ht, &key, &outBlock);
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
    
    pdht_t *ht;
	  //printf("number of blocks %d\n",TOTAL_BLOCKS);
    //printf("Starting\n");
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
    int mine;
    double *zero = calloc(sizeof(double), BLOCK_SIZE * BLOCK_SIZE);

    val = malloc(elemsize);
    memset(val,0,elemsize);
    
    
    //setenv("PTL_IGNORE_UMMUNOTIFY", "1", 1);
    setenv("PTL_PROGRESS_NOSLEEP", "1", 1);
    
    // setup experimental configuration
    pdht_config_t cfg;
    cfg.nptes        = 1;
    cfg.pendmode     = PdhtPendingTrig;
    //cfg.pendmode     = PdhtPendingTriggered;
    cfg.maxentries   = 25000;
    cfg.pendq_size   = 10000;
    cfg.ptalloc_opts = 0;
    cfg.ptalloc_opts = PTL_PT_MATCH_UNORDERED;
    cfg.local_gets = PdhtSearchLocal;
    pdht_tune(PDHT_TUNE_ALL, &cfg);
    
    // create hash table
    ht = pdht_create(sizeof(unsigned long), elemsize, PdhtModeStrict);

    //resultBlock = 
//    pdht_sethash(ht, localhash);
    pdht_fence(ht);
    Block_t *filledBlockHolder;
    filledBlockHolder = calloc(sizeof(Block_t),1);
    

    // Store rows in hash table
    if (c->rank == 0) {
        for (int row = 0; row < BLOCKS_PER_ROW; row++) {
            for (int col = 0; col < BLOCKS_PER_ROW; col++) {
                if (rand_double() < FILL_RATE) {
                    keyA = KEY2_A(row, col);
                    filledBlock(row * BLOCKS_PER_ROW + col,filledBlockHolder);
                    pdht_put(ht, &keyA, filledBlockHolder);
                }
                if (rand_double() < FILL_RATE) {
                    keyB = KEY2_B(row, col);
                    filledBlock(-((row * BLOCKS_PER_ROW + col) + 1),filledBlockHolder);
                    pdht_put(ht, &keyB, filledBlockHolder);
                }
            }
        }
    }
    free(filledBlockHolder);
    pdht_fence(ht);
    if (c->rank == 0) {
        printf("Matrix A\n");
#ifdef PRINT_MATRICES
        printMatrix(&KEY2_A,ht);
#endif
        printf("Matrix B\n");
#ifdef PRINT_MATRICES
        printMatrix(&KEY2_B,ht);
#endif
    }
    
    pdht_fence(ht);
    PDHT_START_ATIMER(total);
       
    resultBlock = calloc(sizeof(Block_t),1);
    for (int row = 0; row < BLOCKS_PER_ROW; row++) {
        for (int col = 0; col < BLOCKS_PER_ROW; col++) {
            mine = (row * BLOCKS_PER_ROW + col) % c->size;
#ifdef MPI
            if(mine % 2 == 1)
              mine--;
#endif
            if (mine == c->rank) {
                KeyOut = KEY2_OUT(row, col);

                multiplyBlocks(row, col,ht,resultBlock);
                if(memcmp(resultBlock->elements, zero, BLOCK_SIZE * BLOCK_SIZE) != 0){
                  pdht_put(ht, &KeyOut, resultBlock);
                }

                //printf("Processor: %d calculating block: (%d, %d)\n", c->rank, row, col);
            }
        }
    }
    
    free(resultBlock);
    PDHT_STOP_ATIMER(total);
    pdht_fence(ht);
    

    
    if (c->rank == 0) {
        printf("Result Matrix\n");
#ifdef PRINT_MATRICES        
        printMatrix(&KEY2_OUT,ht);
#endif    
    }
    
    
    //pdht_print_stats(ht);
    
    eprintf("total elapsed time: %12.7f ns\n", (double)PDHT_READ_ATIMER(total));
    
done:
    pdht_free(ht);
    free(val);
}
