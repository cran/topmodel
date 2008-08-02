#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <R_ext/Utils.h>

#define	ZERO		0.0000001

void topidx(double *map2,
	    int *rows,
            int *cols,
            double *ew_res,
            double *ns_res,
	    double *topidxmap)
{

  double	**atb,**a, **map;
  int	natb = 0; 	//natb = number of NULL values
  int	i,j,k,snatb;
  int	iter,nroute,nslp;
  double	sum,route[9],tanB[9],dx,dx1,dx2,sumtb,C;
  int	nsink = 0;
  int	ncells = (*rows)*(*cols);

  dx     = *ew_res;
  dx1    = 1 / dx;
  dx2    = 1 / (1.414 * dx);
  snatb  = natb;

  /* memory allocation */

  atb = (double **) R_alloc(*rows, sizeof(double *));
  a = (double **) R_alloc(*rows, sizeof(double *));
  map = (double **) R_alloc(*rows, sizeof(double *));

  for(i=0; i<*rows; i++){
    atb[i] = (double *) R_alloc(*cols, sizeof(double));
    a[i] = (double *) R_alloc(*cols, sizeof(double));
    map[i] = (double *) R_alloc(*cols, sizeof(double));
  }

  /* copy input to map */

  for(i=0;i<*rows;i++){
    for(j=0;j<*cols;j++){
      map[i][j] = map2[j+(*cols)*i];
    }
  }

  /* Initialisation */

  for(i=0;i<*rows;i++){
    for(j=0;j<*cols;j++){
      a[i][j]=(*ew_res)*(*ns_res);
      if(map[i][j] < -9000){	// indication of null value
	natb++;				
	atb[i][j] = -9999;  	// null value
      }else{
	atb[i][j] = -10.0;	// no topidx calculated yet
      }
    }
  }

  /* calculation */

  Rprintf("ncells = %i\n", ncells);
  Rprintf("natb = %i\n", natb);

  Rprintf("Iterations:         ");

  for(iter=1;natb<ncells;iter++){

    R_CheckUserInterrupt();
    Rprintf("\b\b\b\b\b\b\b\b%8i",iter);

    for(i=0;i<*rows;i++){
      for(j=0;j<*cols;j++){
	/* skip null values */
	if(map[i][j] < -9000)
	  continue;		// go to next j in for-loop

				/* skip squares already done */
	if(atb[i][j] == -9999 || atb[i][j]>=ZERO)
	  continue;

	/* check the 8 possible flow directions for
	 * upslope elements without an atb value
	 */
	if(i>0){
	  if(j>0 &&
	     (map[i-1][j-1] < -9000    ||	// upper left is NULL (necessary?)
	      map[i-1][j-1]>map[i][j]) &&	// upper left is upslope
	     atb[i-1][j-1] != -9999    &&	// upper left atb is not NULL
	     atb[i-1][j-1]<ZERO)			// upper left atb < 0 (-10)
	    continue;

	  if((map[i-1][j] < -9000    ||		// idem for upper neighbour
	      map[i-1][j]>map[i][j]) &&
	     atb[i-1][j] != -9999    &&
	     atb[i-1][j]<ZERO)
	    continue;

	  if(j+1 < *cols &&			// idem for upper right map
	     (map[i-1][j+1] < -9000      ||		// also check if not on right limit
	      map[i-1][j+1] > map[i][j]) &&
	     atb[i-1][j+1] != -9999      &&
	     atb[i-1][j+1]<ZERO)
	    continue;
	}
	if(j>0 &&				// for left neighbour
	   (map[i][j-1] < -9000    ||
	    map[i][j-1]>map[i][j]) &&
	   atb[i][j-1] != -9999    &&
	   atb[i][j-1]<ZERO)
	  continue;
	if(j+1<*cols &&				// for right neighbour
	   (map[i][j+1] < -9000    ||
	    map[i][j+1]>map[i][j]) &&
	   atb[i][j+1] != -9999    &&
	   atb[i][j+1]<ZERO)
	  continue;
	if(i+1<*rows){
	  if(j>0 &&
	     (map[i+1][j-1] < -9000    ||	// lower left 
	      map[i+1][j-1]>map[i][j]) &&
	     atb[i+1][j-1] != -9999    &&
	     atb[i+1][j-1]<ZERO)
	    continue;
	  if((map[i+1][j] < -9000    ||		// lower neighbour
	      map[i+1][j]>map[i][j]) &&
	     atb[i+1][j] != -9999    &&
	     atb[i+1][j]<ZERO)
	    continue;
	  if(j+1<*cols &&
	     (map[i+1][j+1] < -9000    ||	// lower right
	      map[i+1][j+1]>map[i][j]) &&
	     atb[i+1][j+1] != -9999    &&
	     atb[i+1][j+1]<ZERO)
	    continue;
	}
	/* find the outflow directions and calculate 
	 * the sum of weights
	 */
	sum=0.0;
	for(k=0;k<9;k++)
	  route[k]=0.0;
	nroute=0;
	if(i>0){
	  if(j>0 &&
	     map[i-1][j-1] >= -9000 &&			// upper left is not NULL
	     map[i][j]-map[i-1][j-1]>ZERO){		// upper left higher (samen!)
	    	tanB[0]=(map[i][j]-map[i-1][j-1])*dx2;	// then calculate tanBdx2
	    	route[0]=0.354*dx*tanB[0];		// = dy/4 !?
	    	sum+=route[0];
	    nroute++;
	  }
	  if(map[i-1][j] >= -9000 &&			// upper
	     map[i][j]-map[i-1][j]>ZERO){
	    	tanB[1]=(map[i][j]-map[i-1][j])*dx1;
	    	route[1]=0.5*dx*tanB[1];			// = dy/2 !?
	    	sum+=route[1];
	    nroute++;
	  }
	  if(j+1<*cols &&				// upper right
	     map[i-1][j+1] >= -9000 &&
	     map[i][j]-map[i-1][j+1]>ZERO){
	    	tanB[2]=(map[i][j]-map[i-1][j+1])*dx2;
	    	route[2]=0.354*dx*tanB[2];
	    	sum+=route[2];
	    nroute++;
	  }
	}
	if(j>0 &&					// left
	   map[i][j-1] >= -9000 &&
	   map[i][j]-map[i][j-1]>ZERO){
	  	tanB[3]=(map[i][j]-map[i][j-1])*dx1;
	  	route[3]=0.5*dx*tanB[3];
	  	sum+=route[3];
	  nroute++;
	}
	if(j+1<*cols){					// right
	  if(map[i][j+1] >= -9000 &&
	     map[i][j]-map[i][j+1]>ZERO){
	    	tanB[5]=(map[i][j]-map[i][j+1])*dx1;
	    	route[5]=0.5*dx*tanB[5];
	    	sum+=route[5];
	    nroute++;
	  }
	}
	if(i+1<*rows){					// lower left
	  if(j>0 &&
	     map[i+1][j-1] >= -9000 &&
	     map[i][j]-map[i+1][j-1]>ZERO){
	    	tanB[6] = (map[i][j] - map[i+1][j-1]) * dx2;
	    	route[6]=0.354*dx*tanB[6];
	    	sum+=route[6];
	    nroute++;
	  }
	  if(map[i+1][j] >= -9000 &&			// lower
	     map[i][j]-map[i+1][j]>ZERO){
	    	tanB[7]=(map[i][j]-map[i+1][j])*dx1;
	    	route[7]=0.5*dx*tanB[7];
	    	sum+=route[7];
	    nroute++;
	  }
	  if(j+1<*cols &&				// lower right
	     map[i+1][j+1] >= -9000 &&
	     map[i][j]-map[i+1][j+1]>ZERO){
	    	tanB[8]=(map[i][j]-map[i+1][j+1])*dx2;
	    	route[8]=0.354*dx*tanB[8];
	    	sum+=route[8];
	    nroute++;
	  }
	}

	if(!nroute){
	  /* Rprintf("Sink or boundary node at %d, %d\n",i,j); */
	  nsink++;
	  sumtb=0.0;
	  nslp=0;		// number of slopes
	  if(i>0){
	    if(j>0 &&
	       map[i-1][j-1] >= -9000){
	      	sumtb += (map[i-1][j-1]-map[i][j])*dx2;
	      	nslp++;
	    }
	    if(map[i-1][j] >= -9000){
	      	sumtb += (map[i-1][j]-map[i][j])*dx1;
	      	nslp++;
	    }
	    if(j+1<*cols &&
	       map[i-1][j+1] >= -9000){
	      	sumtb += (map[i-1][j+1]-map[i][j])*dx2;
	      	nslp++;
	    }
	  }

	  if(j>0 &&
	     map[i][j-1] >= -9000){
	    	sumtb+=(map[i][j-1]-map[i][j])*dx1;
	    	nslp++;
	  }
	  if(j+1<*cols &&
	     map[i][j+1] >= -9000){
	    	sumtb+=(map[i][j+1]-map[i][j])*dx1;
	    	nslp++;
	  }
	  if(i+1<*rows){
	    if(j>0 &&
	       map[i+1][j-1] >= -9000){
	      	sumtb+=(map[i+1][j-1]-map[i][j])*dx2;
	      	nslp++;
	    }
	    if(map[i+1][j] >= -9000){
	      sumtb+=(map[i+1][j]-map[i][j])*dx1;
	      nslp++;
	    }
	    if(j+1<*cols &&
	       map[i+1][j+1] >= -9000){
	      	sumtb+=(map[i+1][j+1]-map[i][j])*dx2;
	      	nslp++;
	    }
	  }

	  sumtb/=nslp;
	  if(sumtb>ZERO){
		atb[i][j]=log(a[i][j]/(2*dx*sumtb));
	  }else{
		atb[i][j] = -9999;
	  }
	  natb++;
	  continue;
	}
				
	C=a[i][j]/sum;
	atb[i][j]=log(C); 
/*	atb[i][j]=C; */
	natb++;

	if(i>0){
	  if(j>0){
	    a[i-1][j-1]+=C*route[0];
	  }
	  a[i-1][j]+=C*route[1];
	  if(j+1<*cols){
	    a[i-1][j+1]+=C*route[2];
	  }
	}
	if(j>0){
	  a[i][j-1]+=C*route[3];
	}
	if(j+1<*cols){
	  a[i][j+1]+=C*route[5];
	}
	if(i+1<*rows){
	  if(j>0){
	    a[i+1][j-1]+=C*route[6];
	  }
	  a[i+1][j]+=C*route[7];
	  if(j+1<*cols){
	    a[i+1][j+1]+=C*route[8];
	  }
	}
      }
    }
  }
  Rprintf("\nNumber of sinks or boundaries: %i\n",nsink);

  /* copy atb to output */

  for(i=0;i<*rows;i++){
    for(j=0;j<*cols;j++){
      topidxmap[j+(*cols)*i] = atb[i][j];
    }
  }

  return;
}

