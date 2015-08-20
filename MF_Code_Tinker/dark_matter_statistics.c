#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "header.h"


void output_halo_mass_function()
{
  int k,nr=15;
  double x,dlogm,m,mvir,cdelta,mmin=powf(10,11.833333333333334),mmax=powf(10,14.166666666666659),delta_vir;
  FILE *fp;
  char aa[100];

  fprintf(stderr,"\n\nCALCULATING HALO MASS FUNCTION.\n");
  fprintf(stderr,    "-------------------------------\n\n");

  sprintf(aa,"%s.dndM",Task.root_filename);
  fp = fopen(aa,"w");

  dlogm = (log(mmax) - log(mmin))/(nr-1);

  for(k=0;k<nr;++k)
    {
      m = exp(k*dlogm)*mmin;
      x = halo_mass_function(m);
      fprintf(fp,"%e %e\n",m,x);
    }
  fclose(fp);
}
