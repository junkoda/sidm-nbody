#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

void   force_costevaluate(void);
int    force_getcost_single(void);
int    force_getcost_quadru(void);
void   force_resetcost(void);
void   force_setupnonrecursive(int no);
void   force_treeallocate(int maxnodes, int maxpart);  
int    force_treebuild(void);
int    force_treebuild_single(int startnode, int *typelist, int *creatednodes);
void   force_treeevaluate(int target, double one_over_s_of_a);
void   force_treeevaluate_direct(int target, double one_over_s_of_a);
void   force_treeevaluate_single(int tree, int targetpart, double epsilon);
void   force_treeevaluate_single_BH(int tree, int targetpart, double epsilon);
void   force_treeevaluate_potential(int target);
void   force_treeevaluate_potential_single(int tree, int targetpart, double epsilon);
void   force_treeevaluate_potential_single_BH(int tree, int targetpart, double epsilon);
void   force_treefree(void);
void   force_update_node(int no, int flag);
void   force_update_node_recursive(int no);
void   force_update_size_of_parent_node(int no);


static float  INLINE_FUNC ngb_periodic(float x);
float  ngb_select_closest(int k, int n, float *arr, int *ind);
void   ngb_treeallocate(int npart);
void   ngb_treebuild(void);
float  ngb_treefind(float xyz[3], int desngb, float hguess, int parttype, int **ngblistback, float **r2listback);
int    ngb_treefind_pairs(float xyz[3], float hsml, int **ngblistback, float **r2listback);
int    ngb_treefind_variable(float xyz[3], float hguess, int parttype, int **ngblistback, float **r2listback);
void   ngb_treefree(void);
void   ngb_treesearch(int);
void   ngb_treesearch_pairs(int);
void   ngb_update_nodes(void);






