


double combinedWork(int type, int level, int axis, double xsplit, 
		      double *xmaxleftofguess, double *xminrightofguess,int *flagmemoryimbalance,
		      int    *pabove, int *pbelow, int *plimit);
void   decomposeType(int type);
void   DomainDecomposition(void); 
int    exchangeParticles_A(int type, int recvTask, int axis, double xsplit);
int    exchangeParticles_B(int type, int recvTask, int axis, double xsplit);
void   findExtent(int type, int level, int axis,double *xleft, double *xright);  
void   findFinalExtent(void);
double findSplitPoint(int type, int level, double xleft, double xright, int axis,int *imbalanceflags);
void   getWork(int type, int axis, double xsplit, double *workabove, double *workbelow, 
		 int *particlesabove,int *particlesbelow,double *xmaxleft, double *xminright);






