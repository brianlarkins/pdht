/**
 * @file interface.h
 * @author Mark Hoemmen
 * @since 09 Jun 2006
 * @date 09 Jun 2006
 */


struct sparse_matrix_t;

struct sparse_matrix_t*
sp_load (const char* path, const char* fmt);

int
sp_save (struct sparse_matrix_t* A, const char* path, 
      const char* fmt);

void
sp_format (struct sparse_matrix_t* A);

int
sp_convert (struct sparse_matrix_t* A, const char* type);

struct sparse_matrix_t* 
sp_mult (struct sparse_matrix_t* B, struct sparse_matrix_t* A);

