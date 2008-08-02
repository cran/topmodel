#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <R_ext/Utils.h>

#define	ZERO		0.0000001

void topidx2(double *map2,
	        int *rows,
            int *cols,
            double *ew_res,
            double *ns_res,
	        double *topidxmap)
{

  double	**atb,**a, **map;
  int	natb = 0; 	//natb = number of NULL values (+ done atb values)
  int	i,j,snatb;
  int	iter,nroute,nslp;
  double	sum,route[9],tanB[9],dx,dx1,dx2,sumtb,C;
  int	nsink = 0;
  int	ncells = (*rows)*(*cols);

/* extra for Jim's implementation */
  int   x1,y1,stop,k;
  double dxval,routefac;

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
      if(map[i][j] < -9000){		// indication of null value
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

	// skip null values
	if(map[i][j] < -9000)
	  continue;

	// skip squares already done 
	if(atb[i][j] == -9999 || atb[i][j]>=ZERO)
	  continue;

	/* check the 8 possible flow directions for
	 * upslope elements without an atb value
	 */

	stop = 0;
	for(x1=-1;x1<2;x1++){
	   for(y1=-1;y1<2;y1++){
		// stay inside map:
		if(i+x1 >= 0 && j+y1 >= 0 && i+x1 < *rows && j+y1 < *cols){
		   if(x1 != 0 || y1 != 0){
		      if((map[i+x1][j+y1] > map[i][j])  &&
                         atb[i+x1][j+y1] != -9999      && 
                         atb[i+x1][j+y1] < ZERO)
                       stop = 1;
                   }
                }
           }
        }

        if(stop == 1) continue;

	/* find the outflow directions and calculate 
	 * the sum of weights
	 */

	sum=0.0;
	for(k=0;k<9;k++){
	   route[k]= 0.0;
	   tanB[k] = 0.0;
	   }
	nroute=0;
	k = 0;
	sumtb = 0.0;

	for(x1=-1;x1<2;x1++){
	   for(y1=-1;y1<2;y1++){
		// stay inside map:
		if(i+x1 >= 0 && j+y1 >= 0 && i+x1 < *rows && j+y1 < *cols){
		   if((x1 != 0 || y1 != 0) && map[i+x1][j+y1] >= -9000){
		      if(x1 == 0 || y1 == 0){
			dxval = dx1;
			routefac = 0.5;
		      }
		      else {
			dxval = dx2;
			routefac = 0.354;
		      }
		      if(map[i][j]-map[i+x1][j+y1]>ZERO){
			 tanB[k]=(map[i][j]-map[i+x1][j+y1])*dxval;
	    		 route[k]=routefac*dx*tanB[k];
	    		 sum+=route[k];
			 //sumtb+=tanB[k];
	    		 nroute++;
		      } 
		   }
		}
		k = k+1;
	   }
	}

/* up to here conversion to Jim's version */

	if(!nroute){
	  Rprintf("Sink or boundary node at %i, %i, map value %f\n",i,j,map[i][j]);
	  nsink++;
	  sumtb=0.0;
	  nslp=0;		// number of slopes
	  if(i>0){
	    if(j>0 && map[i-1][j-1] >= -9000){
	       sumtb += (map[i-1][j-1] - map[i][j])*dx2;
	       nslp++;
	    }
	    if(map[i-1][j] >= -9000){
	       sumtb += (map[i-1][j] - map[i][j])*dx1;
	       nslp++;
	    }
	    if(j+1 < *cols && map[i-1][j+1] >= -9000){
	       sumtb += (map[i-1][j+1] - map[i][j])*dx2;
	       nslp++;
	    }
	  }

	  if(j>0 && map[i][j-1] >= -9000){
	    sumtb+=(map[i][j-1] - map[i][j])*dx1;
	    nslp++;
	  }
	  if(j+1<*cols && map[i][j+1] >= -9000){
	    sumtb+=(map[i][j+1]
		    -map[i][j])*dx1;
	    nslp++;
	  }
	  if(i+1<*rows){
	    if(j>0 && map[i+1][j-1] >= -9000){
	      sumtb+=(map[i+1][j-1] - map[i][j])*dx2;
	      nslp++;
	    }
	    if(map[i+1][j] >= -9000){
	      sumtb+=(map[i+1][j] - map[i][j])*dx1;
	      nslp++;
	    }
	    if(j+1<*cols && map[i+1][j+1] >= -9000){
	      sumtb+=(map[i+1][j+1] - map[i][j])*dx2;
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

  /* atb kopieren naar output */

  for(i=0;i<*rows;i++){
    for(j=0;j<*cols;j++){
      topidxmap[j+(*cols)*i] = atb[i][j];
    }
  }

  /* TODO release memory */


  return;
}

