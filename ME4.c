#include <stdio.h>
#include <math.h>
#include <complex.h> 

enum{
	FL=3
};

void prtRMx(double [FL][FL], char *);
void prtRMx(double matr[FL][FL], char *name)
{	
	int i,j;
	printf("\n%s\n",name);
	for(i=0;i<FL;i++)
	{
      for(j=0;j<FL;j++)
	  {
	     printf("%lf",matr[i][j]);
	  }
	  printf("\n");
    }
}

void prtIMx(complex double [FL][FL], char *);
void prtIMx(complex double matr[FL][FL], char *name)
{
	int i,j;
	printf("\n%s\n",name);
	for(i=0;i<FL;i++)
	{
      for(j=0;j<FL;j++)
	  {
	     printf("%lf+i(%lf)",creal(matr[i][j]),cimag(matr[i][j]));
	  }
	  printf("\n");
    }
}

void prtIVr(complex double [FL], char *);
void prtIVr(complex double vect[FL], char *name)
{	
	int i;
	printf("\n%s\n",name);
	for(i=0;i<FL;i++)
	{
	  printf("%lf+i(%lf)",creal(vect[i]),cimag(vect[i]));
	  printf("\n");
    }
}

typedef double (*dens)(double);
typedef struct 
{
	complex double Psi0[FL];
	double H0[FL][FL], W[FL][FL], H0W[FL][FL], WH0W[FL][FL], H0H0W[FL][FL];
	double tol, a, b;
	dens vS;
		
}Mctx;

typedef struct
{
	complex double Psi[FL];
	double Er, dh, last;
}Exctx;

double vS (double);
double vS (double e)
{	
	return 6.5956e4*exp(-10.54*e);
}


void ME4(Mctx *,Exctx *); 
void ME4(Mctx *mctx, Exctx *exctx)
{
		
	int    i,j;
	double s=0.8, ep=0.,em=0.,e=0., vSep,vSem;
	double complex r0=0.+0.*I, r1=0.+0.*I;
	double l[FL],q,p,z,a,b; 											
	double var;
	double complex A[FL][FL],sqA[FL][FL],unit[FL][FL],exA[FL][FL];
	double complex S1[FL][FL],S2[FL][FL],sqS1[FL][FL],arr[FL];
		
		arr[0]=0.; arr[1]=0.; arr[2]=0.;
		
		unit[0][0]=1.; unit[0][1]=0.; unit[0][2]=0.; 					
		unit[1][0]=0.; unit[1][1]=1.; unit[1][2]=0.;
		unit[2][0]=0.; unit[2][1]=0.; unit[2][2]=1.;		
		//начало цикла по точкам															
		e = mctx->a;	
											
		while(e < mctx->b){												
			ep = e+(1.+1./sqrt(3.))*(exctx->dh/2.);
			em = e+(1.-1./sqrt(3.))*(exctx->dh/2.);
			
			vSem = mctx->vS(em);
			vSep = mctx->vS(ep);
			
			//вычисление матрицы OMG4  в  точке e_n(кси_n)
			for(i=0;i<FL;i++){						  				
		      for(j=0;j<FL;j++){
			    A[i][j] = 
			      ((-1.)*(mctx->H0[i][j]+0.5*(vSep+vSem)*mctx->W[i][j])*exctx->dh*I
			     +(sqrt(3.)/12.)*(vSep-vSem)*mctx->H0W[i][j]*exctx->dh*exctx->dh)*(-1.*I);
				}
			}										 				
			
			//след OMG4 Tr(OMG04)=0
			z = (double)(A[0][0]+A[1][1]+A[2][2])/3.; 				
					
			for(i=0;i<FL;i++){		 				 					
		      A[i][i] = A[i][i]-z;
			}										 					
			//Вычисление sqOMG04=OMG04^2
			for(i=0;i<FL;i++){					     				
		      for(j=0;j<FL;j++){
			    sqA[i][j] = 
				  A[i][0]*A[0][j]
				 +A[i][1]*A[1][j]
				 +A[i][2]*A[2][j];
				}
			}										 				
			//cлед sqA^2
			p = (double)(sqA[0][0]+sqA[1][1]+sqA[2][2])/2.; 		
			//det(A)
			q = (double)(
			  A[0][0]*A[1][1]*A[2][2]+A[1][0]*A[0][2]*A[2][1]	   
	     	 +A[2][0]*A[0][1]*A[1][2]-A[2][0]*A[1][1]*A[0][2]
			 -A[0][0]*A[1][2]*A[2][1]-A[1][0]*A[0][1]*A[2][2]);  
			
			//fprintf(stderr,"tol = %lf;\tq = %lf;\tp = %lf;\tz = %lf \n",tol,q,p,z);
			
			//корни характер. многочлена
			l[2] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*0.)/3.);			
			l[1] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*1.)/3.);
			l[0] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*2.)/3.);
			
			//fprintf(stderr,"l0 = %lf l1 = %lf l2 = %lf\n\n",l[0],l[1],l[2]);
																					
			//упорядочиваем корни характер. многочлена
			for( i=0; i < FL; i++) {             					
				for( j = FL-1; j > i; j-- ) {     
					if ( l[j-1] > l[j] ) {
						var = l[j]; l[j] = l[j-1]; l[j-1] = var;
					}
				}
			}														
			
			// a1=a , b1=b смотри статью
			a = l[1]-l[0]; 
			b = l[2]-l[0];
			var = sqrt(p/3.);
			//fprintf(stderr,"a = %lf b = %lf var = %lf\n\n",a,b,var);
			//fprintf(stderr,"\n");
			
			//r0 , r1 для exp(..A)
			r0 = (-1.)*(2.*sin((var*a*exctx->dh)/2.)*sin((var*a*exctx->dh)/2.)
				-I*sin(var*a*exctx->dh))/a;								
				 
			r1 = (-1./(a-b))*((2.*sin((var*a*exctx->dh)/2.)*sin((var*a*exctx->dh)/2.)
			    -I*sin(var*a*exctx->dh))/a
			    -(2.*sin((var*b*exctx->dh)/2.)*sin((var*b*exctx->dh)/2.)
			    -I*sin(var*b*exctx->dh))/b);
			//fprintf(stderr,"r1 = %g + i*(%g), r0 = %g + i*(%g)\n\n",creal(r1),cimag(r1),creal(r0),cimag(r0));
			
			//вычисление exp(OMG04)
			for(i=0;i<FL;i++){										
			  for(j=0;j<FL;j++){
			    exA[i][j] = 
			      exp(var*l[0]*I*exctx->dh)*(
				  (1.-l[0]*(r0-l[1]*r1))*unit[i][j]
				 +(1./var)*(r0+l[2]*r1)*A[i][j]
				 +1./(var*var)*r1*sqA[i][j]);
				}
			}	
			//уравнение эволюции нейтрино в точке en
			for(i=0;i<FL;i++){										
			  exctx->Psi[i] = (
			    exA[i][0]*mctx->Psi0[0]
			   +exA[i][1]*mctx->Psi0[1]
			   +exA[i][2]*mctx->Psi0[2])*exp(exctx->dh*z*I);
			}														
			//S1 и S2 понадобятся для расчета величины следующего шага	
			for(i=0;i<FL;i++)
			{										
			  for(j=0;j<FL;j++)
			  {
			    S1[i][j] = 
			      (-sqrt(3.)/12.)*(vSep-vSem)
				 *mctx->H0W[i][j]*exctx->dh;
			  }
			}
			
			for(i=0;i<FL;i++)
			{
			  for(j=0;j<FL;j++)
			  {
			    S2[i][j] = 
			      (I*sqrt(3.)/24.)*(vSep-vSem)
				 *(mctx->H0H0W[i][j]
				 +(1./2.)*(vSep+vSem)*mctx->WH0W[i][j])*exctx->dh;
			  }
			}		
															
			//sqS1
			for(i=0;i<FL;i++)
			{					    				
			  for(j=0;j<FL;j++){
			    sqS1[i][j] =
			      (S1[i][0]*S1[0][j]
				 +S1[i][1]*S1[1][j]
				 +S1[i][2]*S1[2][j])*exctx->dh*exctx->dh;
			  }														
			}														
			arr[0]=0.+I*0.; 
			arr[1]=0.+I*0.;
			arr[2]=0.+I*0.;
			
			for(i=0;i<FL;i++)
			{										
			  for(j=0;j<FL;j++)
			  {
			    arr[i] += 
			      (S1[i][j]+exctx->dh*S2[i][j]
				 +0.5*exctx->dh*sqS1[i][j])*exctx->Psi[j]*exctx->dh;
			  }
			}
			
			exctx->Er = 0.;
			for(i=0;i<FL;i++)
			{
			  exctx->Er += (double)(arr[i]*conj(arr[i]));
			}														
			
			//fprintf(stderr,"e = %lf\n",e);
			
			//изменение шага exctx->dh
			if(exctx->Er >= mctx->tol)
			{											
			  exctx->dh = s*exctx->dh*pow((mctx->tol/exctx->Er),0.3333);														
			}
			else
			{
			  e = e + exctx->dh;
			  mctx->Psi0[0] = exctx->Psi[0];
			  mctx->Psi0[1] = exctx->Psi[1]; 
			  mctx->Psi0[2] = exctx->Psi[2]; 
			}															 
		}
	    exctx->last = e;
	
}

int main ()
{	
	int i,j;
	complex double sqPsi=0.+I*0.;
	double s12=0.308, s13=0.0234,c12 , c13, E=7.05,
		   q1=4.35196e6, q2 =0.030554,var;
	Mctx mctx;
	Exctx exctx;
	
	mctx.H0[0][0] = 0.; mctx.H0[0][1] = 0.;        mctx.H0[0][2] = 0.;      				
	mctx.H0[1][0] = 0.; mctx.H0[1][1] = (q1*q2)/E; mctx.H0[1][2] = 0.;
	mctx.H0[2][0] = 0.; mctx.H0[2][1] = 0.;        mctx.H0[2][2] = (q1*(1./E));
	
	var = 1. - s12*s12;						         				 
	c12 = sqrt(var);
	var = 1. - s13*s13;
	c13 = sqrt(var);
	
	exctx.Psi[0]=0.+I*0.;
    exctx.Psi[1]=0.+I*0.; 
    exctx.Psi[2]=0.+I*0.;
	exctx.dh=1e-4;
	exctx.Er=0.;
	
	mctx.Psi0[0]=c12*c13+I*0.;
    mctx.Psi0[1]=s12*c13+I*0.; 
    mctx.Psi0[2]=s13+I*0.;
    
    mctx.a = 0.1;
    mctx.b = 1.;
	mctx.tol=1e-4;
	
	mctx.W[0][0] = c13*c13*c12*c12;mctx.W[0][1] = c12*s12*c13*c13; mctx.W[0][2] = c12*c13*s13;
	mctx.W[1][0] = mctx.W[0][1];   mctx.W[1][1] = s12*s12*c13*c13; mctx.W[1][2] = s12*c13*s13;
	mctx.W[2][0] = mctx.W[0][2];   mctx.W[2][1] = mctx.W[1][2];    mctx.W[2][2] = s13*s13;
	
	mctx.vS = vS;
		//Вычисление [H0,W]
		for(i=0;i<FL;i++)
		{					    					
		  for(j=0;j<FL;j++)
		  {
		    mctx.H0W[i][j] = 
			  mctx.H0[i][0]*mctx.W[0][j]-mctx.W[i][0]*mctx.H0[0][j]
			 +mctx.H0[i][1]*mctx.W[1][j]-mctx.W[i][1]*mctx.H0[1][j]
		     +mctx.H0[i][2]*mctx.W[2][j]-mctx.W[i][2]*mctx.H0[2][j];
	      }
		}     					
											
		//Вычисление [H0,[H0,W]]
		for(i=0;i<FL;i++)
		{					    					
		  for(j=0;j<FL;j++)
	      {
		    mctx.H0H0W[i][j] =   
			  mctx.H0[i][0]*mctx.H0W[0][j]-mctx.H0W[i][0]*mctx.H0[0][j]
			 +mctx.H0[i][1]*mctx.H0W[1][j]-mctx.H0W[i][1]*mctx.H0[1][j]
			 +mctx.H0[i][2]*mctx.H0W[2][j]-mctx.H0W[i][2]*mctx.H0[2][j];	
			
		  }															
		}
		//Вычисление [W,[H0,W]]
		for(i=0;i<FL;i++)
		{											
	      for(j=0;j<FL;j++)
		  {
		    mctx.WH0W[i][j] =  
		      mctx.W[i][0]*mctx.H0W[0][j]-mctx.H0W[i][0]*mctx.W[0][j]
			 +mctx.W[i][1]*mctx.H0W[1][j]-mctx.H0W[i][1]*mctx.W[1][j]
			 +mctx.W[i][2]*mctx.H0W[2][j]-mctx.H0W[i][2]*mctx.W[2][j];
		  }
		}
															
		ME4(&mctx,&exctx);
		
		for(i=0;i<FL;i++)
		{
		  sqPsi += exctx.Psi[i]*conj(exctx.Psi[i]);
		}
		
		printf("#b,E,1-|Psi|^2,Psi, dh\n");				
		printf("%lf\t%lf\t |%g\t%g|\t %lf\n",
		  exctx.last,E,1.-creal(sqPsi),cimag(sqPsi),exctx.dh);
		
		prtIVr(exctx.Psi,"Psi");	
		prtRMx(mctx.H0,"H0");
													
return 0;	
}
