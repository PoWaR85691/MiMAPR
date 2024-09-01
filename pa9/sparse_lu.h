#ifndef SPARSE_LU
#define SPARSE_LU

enum Sparse_type { CSR, CSC };

struct Sparse
{
    enum Sparse_type type;
    int size;
    int *line_end;
    int *index;
    double *value;
};

void matrix_to_sparse(double *matrix, int size, struct Sparse *sparse, enum Sparse_type type);
void sparse_to_matrix(struct Sparse *sparse, double *matrix);
void free_sparse(struct Sparse *sparse);
void change_sparse_type(struct Sparse *old_sparse, struct Sparse *new_sparse);

void sparse_lu(struct Sparse *a, struct Sparse *l, struct Sparse *u);
void solve_sparse_lu(struct Sparse *l, struct Sparse *u, double *b, double *x);

#endif