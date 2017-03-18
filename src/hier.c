/*
Minimax linkage agglomerative clustering.
*/

#include <R.h>
#include <stdlib.h>
#define BIGGEST 1e200
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define lt(i,j,n) n*j - j*(j+1)/2 + i - j - 1 // note: must have i > j

typedef struct node {
  int i;
  struct node *next;
} Cluster;

void printLT(double *d,int n,int* clustLabel);
void printCluster(Cluster * G);
void printMatrix(double *dmax,int n,int *clustLabel);
double maxlink(double * d, int n, Cluster * G, Cluster * H);
int findNN(double *dd, int n, int *clustLabel, double *dmax, Cluster** clusters, int *inchain, int i);
double minimaxlink(double *dmax, int n,
		   Cluster *G, Cluster *tG, int iG, Cluster *H, int iH, int nGH);
int minimaxpoint(double *dmax, int n, Cluster *G, int iG, int nG);
double completelink(double *dmax, int n,Cluster *G, int iG, Cluster *H, int iH);

void hier(double *d, int *ndim, int *verbose, int *merge, double *height, int *order, int *protos) {
  int i, j, k, ii, imerge, jmerge, reverse;
  int n = *ndim;
  double *dd = malloc((n*(n-1)/2)*sizeof(double));
  // dd contains the inter-cluster distances as a lower triangular matrix
  
  double *dmax = malloc((n*n)*sizeof(double));
  // dmax[n*i + j] = max_{l\in C_j}d(i,l) where C_j is the jth cluster

  int *un_protos = malloc(n*sizeof(int));
  int *un_merge = malloc((n*2)*sizeof(int));
  // versions of the outputs "protos" and "merge"... except unordered by height
  // these will be ordered based on the nearest-neighbor chain order
  // at the end, we will convert from un_protos to protos and likewise for merge


  // each cluster will be stored as a linked list between clusters[i] and tails[i]
  Cluster** clusters = malloc(n*sizeof(Cluster *));
  Cluster** tails = malloc(n*sizeof(Cluster *));
  
  int clustLabel[n];
  int clustSize[n];

  int nnchain[n];
  int nn,end = -1; // index of end of nnchain
  int inchain[n];// indicates whether element i is in current chain

  double old;
  int nochange=0;

  // initialize n singleton clusters
  for(i = 0; i < n; i++) {
    clusters[i] = malloc(sizeof(Cluster));
    tails[i] = clusters[i];
    clusters[i]->i = i;
    clusters[i]->next = NULL;
    clustLabel[i] = -(i+1);// leaves are negative (as in R's hclust)
    clustSize[i] = 1;
    nnchain[i] = -1;
    inchain[i] = 0;
  }
  for(i = 0; i < n; i++)
    dmax[n*i + i] = 0;
  for(j = 0; j < n-1; j++) {
    for(i = j+1; i < n; i++) {
      ii = lt(i,j,n);
      dmax[n*i + j] = d[ii];
      dmax[n*j + i] = d[ii];
    }
  }
  
  for(ii = 0; ii < n*(n-1)/2; ii++)
    dd[ii] = d[ii];
  
  for(k = 0; k < n-1; k++) {
    // check if chain is empty
    if(end == -1) {
      // start a new chain
      for(i = 0; i < n; i++)
	if(clustLabel[i] != 0)
	  break;
      nnchain[0] = i;
      if(*verbose)
	Rprintf("\nStarting a new chain at %d.\n",i+1);
      end = 0;
    }
    
    // grow nearest neighbor chain until a RNN pair found:
    while(TRUE) {
      nn = findNN(dd,n,clustLabel,dmax,clusters,inchain,nnchain[end]);
      if(end > 0) {
	if(nn == nnchain[end-1]) {
	  // reached a RNN pair
	  break;
	}
	inchain[nnchain[end-1]] = 1;
	// we purposely delay this update to allow findNN to 
	// select a cluster visited a step earlier
	// (which happens for RNN pairs)
      }
      nnchain[++end] = nn;
    }
    if(*verbose) {
      Rprintf(" NN-Chain: ");
      for(i = 0; i <= end; i++)
	Rprintf("%d ",nnchain[i]+1);
      Rprintf("\n");
    }
    if(nnchain[end] < nnchain[end-1]) {
      imerge = nnchain[end-1];
      jmerge = nnchain[end];
    }
    else {
      imerge = nnchain[end];
      jmerge = nnchain[end-1];
    }
    
    // remove RNN pair from chain
    inchain[nnchain[end]] = 0;
    inchain[nnchain[end-1]] = 0;
    
    if(end > 1)
      inchain[nnchain[end-2]] = 0;// again, so that a RNN can be detected
    if(end > 2)
      inchain[nnchain[end-3]] = 0;// again, so that a RNN can be detected
    
    nnchain[end] = -1;
    nnchain[end-1] = -1;
    end -= 2;
    
    // create merged cluster from this RNN pair:
    ii = lt(imerge,jmerge,n);
    height[k] = dd[ii];
    if(*verbose)
      Rprintf("  Merged reciprocal nearest neighbor pair at height %g\n",height[k]);
    reverse = 0;
    
    if(clustLabel[imerge] > clustLabel[jmerge])
      reverse = 1;// put smaller cluster on left
    
    if(clustLabel[imerge] < 0 && clustLabel[jmerge] < 0)
      reverse = 1;  // since imerge > jmerge
    
    if(reverse) {
      un_merge[2*k] = clustLabel[jmerge];
      un_merge[2*k+1] = clustLabel[imerge];
    }
    else {
      un_merge[2*k] = clustLabel[imerge];
      un_merge[2*k+1] = clustLabel[jmerge];
    }
      
    // update the imerge cluster:
    clustLabel[imerge] = k + 1;
    clustSize[imerge] += clustSize[jmerge];
      
    if(reverse)	{
      // put jmerge's elements to left of imerge's
      tails[jmerge]->next = clusters[imerge];
      clusters[imerge] = clusters[jmerge];
    }
    else {
      // put imerge's elements to left of jmerge's
      tails[imerge]->next = clusters[jmerge];
      tails[imerge] = tails[jmerge];
    }
      
    // at jmerge, we no longer have a cluster:
    clustLabel[jmerge] = 0;
    clustSize[jmerge] = 0;
    clusters[jmerge] = NULL;
    tails[jmerge] = NULL;
    
    /// update dmax:
    for(i = 0; i < n; i++) {
      // point i
      if(dmax[n*i+imerge] < dmax[n*i+jmerge]) {
	dmax[n*i+imerge] = dmax[n*i+jmerge];//i.e. max{dmax(i,imerge),dmax(i,jmerge)}
      }
    }
    
    // get the minimax prototype for this newly formed cluster
    un_protos[k] = minimaxpoint(dmax,n,clusters[imerge],imerge,clustSize[imerge]) + 1;
    
    /// update dd:
    
    // imerge is now a new cluster, so update its distances:
    for(j = 0; j < imerge; j++)
      if(clustLabel[j] != 0) { // still an active cluster
	ii = lt(imerge,j,n);
	old = dd[ii];
	dd[ii] = minimaxlink(dmax,n,clusters[j],tails[j],j,clusters[imerge],imerge,
			     clustSize[j]+clustSize[imerge]);
	if(dd[ii]==old)
	  nochange++;
      }
      
      ii = lt(imerge+1,imerge,n);
      for(i = imerge + 1; i < n; i++) {
	if(clustLabel[i] != 0) {
	  old = dd[ii];
	  dd[ii] = minimaxlink(dmax,n,clusters[i],tails[i],i,clusters[imerge],imerge,
			       clustSize[i] + clustSize[imerge]);
	  if(dd[ii]==old)
	    nochange++;
	}
	ii++;
      }
  }

  //// List merges and protos in order of increasing height:
  int o[n-1];
  for(i = 0; i < n-1; i++)
    o[i] = i;
    
  // sort heights and "o = order(height)" (in R speak)
  rsort_with_index(height,o,n-1);

  // if there are ties, want indices ordered (to match R's convention).
  int count;
  
  for(i = 0; i < n-1; i++) {
    count = 0;
    if((i+count+1)<(n-1))
      while(height[i+count+1]==height[i])
	count++;
    // heights are constant from i to i+count
    if(count > 0) {
      R_isort(&o[i],count+1);
    }
    i += count;
  }

  for(i = 0; i < n-1; i++) {
    protos[i] = un_protos[o[i]];
    merge[2*i] = un_merge[2*o[i]];
    merge[2*i+1] = un_merge[2*o[i]+1];
  }

  // shuffling merge rows around messes up the positive indices in merge:
  int ranks[n-1];
  for(i = 0; i < n-1; i++)
    ranks[o[i]] = i;

  int mtemp;
  for(i = 0; i < n-1; i++) {
    for(j = 0; j < 2; j++) {
      if(merge[2*i+j] > 0)
	merge[2*i+j] = ranks[merge[2*i+j]-1] + 1;
      //the -1 and +1 are to match R's indexing
    }
    if(merge[2*i] > 0 && merge[2*i+1] > 0)
      if(merge[2*i] > merge[2*i+1]) {
	// hclust has positive rows in increasing order:
	mtemp = merge[2*i];
	merge[2*i] = merge[2*i+1];
	merge[2*i+1] = mtemp;
      }
  }
  // get order by following "merge"

  // using ranks for a different purpose... ranks[k] will be s.t. clusters[ranks[k]] contains
  // the cluster created at step k.
  for(k=0; k < n-1; k++)
    ranks[k] = 0;

  Cluster *cur = clusters[imerge];
  for(i = 0; i < n; i++) {
    clusters[cur->i] = cur;
    tails[cur->i] = clusters[cur->i];
    cur = cur->next;
  }
  for(i = 0; i < n; i++)
    clusters[i]->next = NULL;

  for(k = 0; k < n-1; k++) {
    if(merge[2*k] < 0)
      imerge = -merge[2*k] - 1;
    else
      imerge = ranks[merge[2*k] - 1];
    if(merge[2*k+1] < 0)
      jmerge = -merge[2*k+1] - 1;
    else
      jmerge = ranks[merge[2*k+1] - 1];
    
    tails[imerge]->next = clusters[jmerge];
    tails[imerge] = tails[jmerge];
    clusters[jmerge] = NULL;
    tails[jmerge] = NULL;
    ranks[k] = imerge;
  }
  
  cur = clusters[imerge];
  for(i = 0; i < n; i++) {
    order[i] = cur->i+1;
    cur = cur->next;
  }

  cur = clusters[imerge];
  Cluster *curnext;
  for(i = 0; i < n; i++) {
    curnext = cur->next;
    free(cur);
    cur = curnext;
  }
  free(dd);
  free(dmax);
  free(un_protos);
  free(un_merge);
  free(clusters);
  free(tails);
}

// returns the nearest neighbor cluster of cluster i that is not 
// already in the chain.
int findNN(double *dd, int n, int *clustLabel, double *dmax, Cluster **clusters, int *inchain, int i)
{
  int j,ii;
  double mind = BIGGEST;
  double dcomp, mincomplete = 0;
  int nn;
  for(j = 0; j < i; j++)
    {
      if(clustLabel[j] == 0 || inchain[j]==1)
	continue;
      ii = lt(i,j,n);
      
      if(dd[ii] < mind)
	{
	  mind = dd[ii];
	  nn = j;
	  mincomplete = 0;// reset mincomplete
	}
      else if(dd[ii]==mind)
	{
	  if(mincomplete==0)
	    {
	      // this is the first duplicate
	      mincomplete = completelink(dmax, n, clusters[nn], nn, clusters[i], i);
	    }
	  dcomp = completelink(dmax, n, clusters[j], j, clusters[i], i);
	  if(dcomp < mincomplete)
	    {
	      mincomplete = dcomp;
	      nn = j;
	    }
	}
    }
    
  for(j = i+1; j < n; j++)
    {
      if(clustLabel[j] == 0 || inchain[j]==1)
	continue;
      ii = lt(j,i,n);
      if(dd[ii] < mind)
	{
	  mind = dd[ii];
	  nn = j;
	  mincomplete = 0;// reset mincomplete
	}
      else if(dd[ii]==mind)
	{
	  if(mincomplete==0)
	    {
	      // this is the first duplicate
	      mincomplete = completelink(dmax, n, clusters[nn], nn, clusters[i], i);
	    }
	  dcomp = completelink(dmax, n, clusters[j], j, clusters[i], i);
	  if(dcomp < mincomplete)
	    {
	      mincomplete = dcomp;
	      nn = j;
	    }
	}
    }

  return nn;
}


// Returns the minimax distance
double minimaxlink(double *dmax, int n,
		   Cluster *G, Cluster *tG, int iG, Cluster *H, int iH,int nGH)
{
  //printf("Inside dmax: (%d,%d)\n",iG,iH);
  //printCluster(G);
  //printCluster(H);

  // temporarily combine clusters
  tG->next = H;
  int i;
  double dmm;
  double dmaxGH[nGH];

  Cluster *cur1;
  cur1 = G;
  for(i = 0; i < nGH; i++)
    {
      if(dmax[n*cur1->i + iG] > dmax[n*cur1->i + iH])
	dmaxGH[i] = dmax[n*cur1->i + iG];
      else
	dmaxGH[i] = dmax[n*cur1->i + iH];
      cur1 = cur1->next;
    }

  dmm = BIGGEST;
  for(i = 0; i < nGH; i++)
    {
      if(dmaxGH[i] < dmm)
	dmm = dmaxGH[i];
    }
 
  // uncombine the clusters
  tG->next = NULL;
  return dmm;
}

// Finds the minimax point of the cluster G
int minimaxpoint(double *dmax, int n, Cluster *G, int iG, int nG)
{
  int mm;
  Cluster *cur1 = G;
  double dmm = BIGGEST;
  int i;
  for(i = 0; i < nG; i++)
    {
      //      printf("dmax[%d,%d]=%.2g\n",cur1->i,iG,dmax[n*cur1->i+iG]);
      if(dmax[n*cur1->i+iG] < dmm)
	{
	  dmm = dmax[n*cur1->i+iG];
	  mm = cur1->i;
	} 
      cur1 = cur1->next;
    }
  return mm;
} 

// complete linkage d(G,H)
double completelink(double *dmax, int n, Cluster *G, int iG, Cluster *H, int iH)
{
  //printf("\nComplete linkage:\n");
  //printCluster(G);
  //printCluster(H);

  double dmm = 0;
  Cluster *cur = G;
  while(cur != NULL)
    {
      // dmax(g,H)
      if(dmax[n*cur->i + iH] > dmm)
	dmm = dmax[n*cur->i + iH];
      cur = cur->next;
    }
  cur = H;
  while(cur != NULL)
    {
      // dmax(h,G)
      if(dmax[n*cur->i + iG] > dmm)
	dmm = dmax[n*cur->i + iG];
      cur = cur->next;
    }
  //  printf("%.8g\n",dmm);
  return dmm;
}

// prints a lower triangular matrix
void printLT(double *d,int n,int* clustLabel)
{
  int i,j;
  
  for(j = 0; j < n; j++)
    Rprintf("\t%d",j);
  Rprintf("\n");
  for(i = 1; i < n; i++)
    {
      Rprintf("%d\t",i);
      if(clustLabel[i]==0)
	{
	  for(j = 0; j < i; j++)
	    Rprintf("*\t");
	  Rprintf("\n");
	  continue;
	}
      for(j = 0; j < i; j++)
	{
	  if(clustLabel[j]==0)
	    Rprintf("*\t");
	  else
	    Rprintf("%.2g\t",d[lt(i,j,n)]);
	}
      Rprintf("\n");
    }
}

void printCluster(Cluster * G)
{
  while(G!=NULL)
    {
      Rprintf("%d\t",G->i);
      G = G->next;
    }
  Rprintf("\n");
}

void printMatrix(double *dmax,int n,int *clustLabel)
{
  int i, j;
  for(i = 0; i < n; i++)
    {
      for(j = 0; j < n; j++)
	{
	  if(clustLabel[j]==0)
	    Rprintf("*\t");
	  else
	    Rprintf("%.2g\t",dmax[n*i + j]);
	}
      Rprintf("\n");
    }
}

/*
int main(int argc, char** argv)
{
  return 1;
}
*/
