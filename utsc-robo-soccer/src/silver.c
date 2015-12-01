/******************************************************************************
*   This C/C++ function solves a linear system of equations 
*    
*   ΣkAi,kxk=yi
*
*   for the unknown vector xk, with yi the given vector (the data) and Ai,k
*   the given coefficient matrix.
*   This is achieved by algebraically solving the equation for the component xi, i.e.
*
*   xi=(yi -Σk(≠i)Ai,kxk)/Ai,i 
*
*   and then iterating this equation using an initial guess for xk. However, in order to 
*   stabilize the algorithm in case of ill-conditioned matrices Ai,k an additional vector ai.xi
*   is added to the original equation, which thus results in the iteration equation
*
*   xi(n)=(yi +ai.xi(n-1) -Σk(≠i)Ai,kxk(n-1))/(Ai,i+ai)  .
*
*   It was found that convergence of the iteration is achieved even for ill-conditioned cases 
*   if ai is chosen such that Ai,i+ai results in a strictly diagonally dominant matrix.
*   The best way (i.e. providing the best compromise between accuracy of the retrieval and number 
*   of iterations) proved indeed to be to set ai equal to the sum over the absolute values 
*   of the off-diagonal row elements, i.e.
*
*   ai= Σk(≠i)|Ai,k|
*
*
*   Input parameters are:
*    N:  number of vector elements, and dimension of the coefficient array
*    y:  vector of data points (with dimension N)
*    A1: pointer to first element of the coefficient array
*    epsi: desired relative accuracy (should be chosen not smaller than 
*                                     the relative noise in the data y)  
*
*   Output parameters are:
*    x: solution vector (with dimension N)
*    ni: number of iterations
*
*
*   The function is called in the form    inv(N,y,&A[0][0],epsi,x,ni);
*
*    
*
*   Note that the function not only checks whether the solution vector x has converged,
*   but also whether the resulting model data have converged to the original data y (the
*   latter should actually be considered the primary criterion, as there may not be a 
*   unique solution for x considering the given noise level (the solution should be 
*   in particular considered with caution in this sense if convergence was achieved 
*   already after 1 iteration)).
*   In case the algorithm does not converge, it might help to rearrange the order of
*   the indexing for the matrix Ai,k (e.g. to AN-i,k). In any case, the matrix should 
*   be reasonably well behaved (random coefficients are unlikely to lead to a 
*   convergence of the algorithm; also note that I have not tested at all whether
*   the algorithm actually works also in case of mixed positive and negative coefficients 
*   (even though I have formally taken this case into account for the conditioning of the
*   iteration matrix)).
*
*   This is a rather new algorithm I have developed, and I have tested this so far only
*   in context with a radiative transfer inversion problem (where it has worked
*   reasonably well; see http://www.plasmaphysics.org.uk/retrieval.htm ). I would
*   therefore appreciate any feedback regarding its suitability in other contexts 
*   (for contact details see  http://www.plasmaphysics.org.uk/feedback.htm ). 
*   
*
*                                               Thomas Smid, April 2008
*********************************************************************************/

 void inv(int N, double y[], double *A1, double epsi, double x[], int& ni)  {

    int i,k;
    double A[N][N],a[N],xx[N],yy[N],convgc,dev,devmax;

     for (i=0; i<=N-1; i=i+1)  {
       for (k=0; k<=N-1; k=k+1)  {
          A[i][k]=A1[N*i+k];
       }
     }

     for (i=0; i<=N-1; i=i+1)  {
       a[i]=0;
       for (k=0; k<=N-1; k=k+1)  {
       if (k!=i) {a[i]=a[i]+fabs(A[i][k]);}
       }
       if (A[i][i]<=0) {a[i]=1.1*a[i]+fabs(A[i][i]);}
       A[i][i]=A[i][i]+a[i];
       x[i]=y[i]/A[i][i];
     }

     convgc=1; ni=0;

     while (convgc >= epsi) {

       for (i=0; i<=N-1; i=i+1)  {
        xx[i]=x[i];
       }

       for (i=0; i<=N-1; i=i+1)  {
         x[i]=0;
         for (k=0; k<=N-1; k=k+1)  {
           if (k!=i) { x[i]=x[i] +xx[k]*A[i][k]; }
         }
         x[i]=(y[i]+a[i]*xx[i]-x[i])/A[i][i] ;
       }

           
      for (i=0; i<=N-1; i=i+1)  {
       yy[i]=0;
       for (k=0; k<=N-1; k=k+1)  {
         yy[i]=yy[i]+ x[k]*A[i][k];
       }
         yy[i]=yy[i]- x[i]*a[i];
      }



      devmax=0;
      for (i=1; i<=N-1; i=i+1)  {
       if (xx[i]!=0) { dev=fabs(x[i]/xx[i]-1.0); if (dev>devmax) devmax=dev; }
       if (y[i]!=0)  { dev=fabs(yy[i]/y[i]-1.0); if (dev>devmax) devmax=dev; }
      }

      convgc=devmax;
      ni=ni+1;

   }

 }