#include <stdio.h>
#include <math.h>
#include <complex.h> 

enum
{
  FL=3,
  NORM,
  PREV,
  LAST,
  OVER
};

const char *pref = "#"; 

typedef double (*dens)(double);

typedef struct 
{
  double complex Psi0[FL];
  double H0[FL][FL], W[FL][FL], H0W[FL][FL], WH0W[FL][FL], H0H0W[FL][FL];
  double tol, a, b;
  dens vS;
  
}basic_ctx;

typedef struct
{
  double complex Psi[FL];
  double Er, dh, last;
  
}variative_ctx;

double atof (const char *str); 
double vS (double); 
void lambda (double [FL], double, double);
void prt_matr(double [FL][FL], char *);
void prt_Cmatr(double complex [FL][FL], char *);
void prt_Cvect(double complex [FL], char *);
void ME4(const basic_ctx *,variative_ctx *); 

int main (int argc, char* argv[])
{
  int i,j;
  double complex sqPsi;
  double s12, s13, c12, c13, E, PSurv, q1, q2;

  E=atof(argv[1]);
  PSurv=0.;
  sqPsi=0.+I*0.;
  s12=fabs(sqrt(0.308)); s13=fabs(sqrt(0.0234));
  c12=fabs(sqrt(1. - 0.308)); c13=fabs(sqrt(1. - 0.0234));
  q1=4.35196e6; q2 =0.030554;
  
  basic_ctx basic;
  variative_ctx vartve;

  basic.H0[0][0] = 0.; basic.H0[0][1] = 0.;        basic.H0[0][2] = 0.;
  basic.H0[1][0] = 0.; basic.H0[1][1] = (q1*q2)/E; basic.H0[1][2] = 0.;
  basic.H0[2][0] = 0.; basic.H0[2][1] = 0.;        basic.H0[2][2] = q1/E;
  
  vartve.Psi[0]=0.+I*0.;
  vartve.Psi[1]=0.+I*0.; 
  vartve.Psi[2]=0.+I*0.;
  vartve.dh=1e-4;
  vartve.Er=0.;

  basic.Psi0[0]=c12*c13+I*0.;
  basic.Psi0[1]=s12*c13+I*0.; 
  basic.Psi0[2]=s13+I*0.;

  basic.tol = atof(argv[2]);//e-4
  basic.a = atof(argv[3]);//0.1
  basic.b = atof(argv[4]);//1

  basic.W[0][0] = c13*c13*c12*c12; basic.W[0][1] = c12*s12*c13*c13; basic.W[0][2] = c12*c13*s13;
  basic.W[1][0] = basic.W[0][1];   basic.W[1][1] = s12*s12*c13*c13; basic.W[1][2] = s12*c13*s13;
  basic.W[2][0] = basic.W[0][2];   basic.W[2][1] = basic.W[1][2];    basic.W[2][2] = s13*s13;

  basic.vS = vS;
  //Вычисление [H0,W]
  for(i=0;i<FL;i++)
  {
    for(j=0;j<FL;j++)
    {
      basic.H0W[i][j] = 
        basic.H0[i][0]*basic.W[0][j]-basic.W[i][0]*basic.H0[0][j]
       +basic.H0[i][1]*basic.W[1][j]-basic.W[i][1]*basic.H0[1][j]
       +basic.H0[i][2]*basic.W[2][j]-basic.W[i][2]*basic.H0[2][j];
    }
  }
  //Вычисление [H0,[H0,W]]
  for(i=0;i<FL;i++)
  {
    for(j=0;j<FL;j++)
    {
      basic.H0H0W[i][j] =   
        basic.H0[i][0]*basic.H0W[0][j]-basic.H0W[i][0]*basic.H0[0][j]
       +basic.H0[i][1]*basic.H0W[1][j]-basic.H0W[i][1]*basic.H0[1][j]
       +basic.H0[i][2]*basic.H0W[2][j]-basic.H0W[i][2]*basic.H0[2][j];
    }
  }
  //Вычисление [W,[H0,W]]
  for(i=0;i<FL;i++)
  {
    for(j=0;j<FL;j++)
    {
      basic.WH0W[i][j] =  
        basic.W[i][0]*basic.H0W[0][j]-basic.H0W[i][0]*basic.W[0][j]
       +basic.W[i][1]*basic.H0W[1][j]-basic.H0W[i][1]*basic.W[1][j]
       +basic.W[i][2]*basic.H0W[2][j]-basic.H0W[i][2]*basic.W[2][j];
    }
  }
  printf("# model:sun\n# model parameters:\t%lf\t%lf\n",65956.000,-10.540);
  printf("%s a = %e\n%s b = %e\n%s tol = %e\n%s E=%e\n%s dh = %e\n",
    pref,basic.a,pref,basic.b,pref,basic.tol,pref,E,pref,vartve.dh);
  printf("%s s12 = %e\n%s s13 = %e\n%s s12^2 = %e\n%s s13^2 = %e\n",
    pref,s12,pref,s13,pref,0.308,pref,0.0234);
  prt_Cvect(basic.Psi0,"Psi0");	
  prt_matr(basic.H0,"H0");
  prt_matr(basic.W,"W");
  prt_matr(basic.H0W,"H0W");
  prt_matr(basic.H0H0W,"H0H0W");
  prt_matr(basic.WH0W,"WH0W");
  printf("\n");
  printf("#################################################\n");
  printf("##            CALCULATION COMPLETED            ##\n");
  printf("#################################################\n");
  
  ME4(&basic,&vartve);
  
  PSurv = c12*c12*c13*c13*(vartve.Psi[0]*conj(vartve.Psi[0]))
         +s12*s12*c13*c13*(vartve.Psi[1]*conj(vartve.Psi[1]))
         +s13*s13*(vartve.Psi[2]*conj(vartve.Psi[2]));
  
  for(i=0;i<FL;i++)
  {
    sqPsi += vartve.Psi[i]*conj(vartve.Psi[i]);
  }
  prt_Cvect(vartve.Psi,"Psi");
  printf("#ProbSurv b 1-|Psi|^2 dh(last)\n");
  printf("%g\t%g\t%g\t%g\n",PSurv,vartve.last,1.-creal(sqPsi),vartve.dh);


return 0;
}

double vS (double e)
{
  return 65956.*exp(-10.54*e);
}
void lambda (double l[FL], double q, double p)
{
  int i,j;
  double var;
  l[2] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*0.)/3.);
  l[1] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*1.)/3.);
  l[0] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*2.)/3.);
  
  for( i=0; i < FL; i++) 
  {
    for( j = FL-1; j > i; j-- )
    {     
      if( l[j-1] > l[j] ) 
      {
        var = l[j]; l[j] = l[j-1]; l[j-1] = var;
      }
    }
  }   
}
void prt_matr(double matr[FL][FL], char *name)
{
  int i,j;
  printf("%s %s = ",pref,name);
  for(i=0;i<FL;i++)
  {
    printf("\n%s\t",pref);
    for(j=0;j<FL;j++)
    {
      printf(" %e",matr[i][j]);
    }
  }
  printf("\n");
}
void prt_Cmatr(complex double matr[FL][FL], char *name)
{
  int i,j;
  printf("\n#%s",name);
  for(i=0;i<FL;i++)
  {
    printf("\n%s\t",pref);		
    for(j=0;j<FL;j++)
    {
      printf("%lf+i(%lf)",creal(matr[i][j]),cimag(matr[i][j]));
    }
  }
}
void prt_Cvect(complex double vect[FL], char *name)
{
  int i;
  printf("%s %s = \n",pref,name);
  for(i=0;i<FL;i++)
  {
    printf("%s\t%e+i(%e)\n",pref,creal(vect[i]),cimag(vect[i]));
  }
}
void ME4(const basic_ctx *basic, variative_ctx *vartve)
{
  int    i,j, index;
  double s, ep,em,e, vSep,vSem;
  double complex r0, r1;
  double l[FL],q,p,z,a,b;
  double var;
  double complex A[FL][FL],sqA[FL][FL],unit[FL][FL],exA[FL][FL];
  double complex S1[FL][FL],S2[FL][FL],sqS1[FL][FL],arr[FL],Psi[FL];
  
  index=NORM;
  s=0.8;
  ep=0.; em=0.; e=0.;
  r0=0.+0.*I; r1=0.+0.*I;
  
  Psi[0]=basic->Psi0[0];
  Psi[1]=basic->Psi0[1];
  Psi[2]=basic->Psi0[2];
  
  arr[0]=0.; arr[1]=0.; arr[2]=0.;

  unit[0][0]=1.; unit[0][1]=0.; unit[0][2]=0.;
  unit[1][0]=0.; unit[1][1]=1.; unit[1][2]=0.;
  unit[2][0]=0.; unit[2][1]=0.; unit[2][2]=1.;
  //начало цикла по точкам
  e = basic->a;

  while(index != OVER)
  { 
    ep = e+(1.+1./sqrt(3.))*(vartve->dh/2.);
    em = e+(1.-1./sqrt(3.))*(vartve->dh/2.);

    vSem = basic->vS(em);
    vSep = basic->vS(ep);

    //вычисление матрицы OMG4  в  точке e_n(кси_n)
    for(i=0;i<FL;i++)
    {
      for(j=0;j<FL;j++)
      {
        A[i][j] = 
          -(basic->H0[i][j]
         +(1./2.)*(vSep+vSem)*basic->W[i][j])
         -I*(sqrt(3.)/12.)*(vSep-vSem)*basic->H0W[i][j]*vartve->dh;
      }
    }

    //след OMG4 Tr(OMG04)=0
    z = (double)(A[0][0]+A[1][1]+A[2][2])/3.;

    for(i=0;i<FL;i++)
    {
      A[i][i] = A[i][i]-z;
    }
    //Вычисление sqOMG04=OMG04^2
    for(i=0;i<FL;i++)
    {
      for(j=0;j<FL;j++)
      {
        sqA[i][j] = 
          A[i][0]*A[0][j]
         +A[i][1]*A[1][j]
         +A[i][2]*A[2][j];
      }
    }
    //cлед sqA^2
    p = (double)(sqA[0][0]+sqA[1][1]+sqA[2][2])/2.;
    //det(A)
    q = (double)
      (A[0][0]*A[1][1]*A[2][2]+A[1][0]*A[0][2]*A[2][1]	   
      +A[2][0]*A[0][1]*A[1][2]-A[2][0]*A[1][1]*A[0][2]
      -A[0][0]*A[1][2]*A[2][1]-A[1][0]*A[0][1]*A[2][2]);
      
    //корни характер. многочлена//упорядочиваем корни характер. многочлена
    lambda(l,q,p);
    
    // a1=a , b1=b смотри статью
    a = l[1]-l[0]; 
    b = l[2]-l[0];
    var = sqrt(p/3.);

    //r0 , r1 для exp(..A)
    r0 = (-1.)*(2.*sin((var*a*vartve->dh)/2.)*sin((var*a*vartve->dh)/2.)
        -I*sin(var*a*vartve->dh))/a;

    r1 = (1./(b-a))*((2.*sin((var*a*vartve->dh)/2.)*sin((var*a*vartve->dh)/2.)
        -I*sin(var*a*vartve->dh))/a
        -(2.*sin((var*b*vartve->dh)/2.)*sin((var*b*vartve->dh)/2.)
        -I*sin(var*b*vartve->dh))/b);

    //вычисление exp(OMG04)
    for(i=0;i<FL;i++)
    {
      for(j=0;j<FL;j++)
      {
        exA[i][j] = 
          ((1.-l[0]*(r0-l[1]*r1))*unit[i][j]
          +(1./var)*(r0+l[2]*r1)*A[i][j]
          +1./(var*var)*r1*sqA[i][j])*exp(var*l[0]*I*vartve->dh);
      }
    }
    //уравнение эволюции нейтрино в точке en
    for(i=0;i<FL;i++)
    {
      vartve->Psi[i] = 
        (exA[i][0]*Psi[0]
        +exA[i][1]*Psi[1]
        +exA[i][2]*Psi[2])*exp(vartve->dh*z*I);
    }
    //S1 и S2 понадобятся для расчета величины следующего шага	
    for(i=0;i<FL;i++)
    {
      for(j=0;j<FL;j++)
      {
        S1[i][j] = 
          (-sqrt(3.)/12.)*(vSep-vSem)*basic->H0W[i][j]*vartve->dh;
      }
    }

    for(i=0;i<FL;i++)
    {
      for(j=0;j<FL;j++)
      {
        S2[i][j] = 
          (basic->H0H0W[i][j]
         +(1./2.)*(vSep+vSem)*basic->WH0W[i][j])*vartve->dh*(I*sqrt(3.)/24.)*(vSep-vSem);
      }
    }

    //sqS1
    for(i=0;i<FL;i++)
    {
      for(j=0;j<FL;j++)
      {
        sqS1[i][j] =
          (S1[i][0]*S1[0][j]
          +S1[i][1]*S1[1][j]
          +S1[i][2]*S1[2][j])*vartve->dh*vartve->dh;
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
          (S1[i][j]+vartve->dh*S2[i][j]
         +(1./2.)*vartve->dh*sqS1[i][j])*vartve->Psi[j]*vartve->dh;
      }
     }

    vartve->Er = 0.;
    for(i=0;i<FL;i++)
    {
      vartve->Er += (double)(arr[i]*conj(arr[i]));
    }
	
    //изменение шага vartve->dh
    if(vartve->Er >= basic->tol)
    {
      vartve->dh = s*vartve->dh*pow((basic->tol/vartve->Er),1./3.);
    }
    else
    {
      e = e + vartve->dh;
      Psi[0] = vartve->Psi[0];
      Psi[1] = vartve->Psi[1]; 
      Psi[2] = vartve->Psi[2];
      
      if(LAST == index)
      {
        index = OVER;
      }
      if(index == PREV)
      {
        index = LAST;
      }
      if((e+2.*vartve->dh)>1. && index != LAST && index != OVER)
      {
        index = PREV;
        vartve->dh = (1.-e)/2.;
      }
    }
  }
  vartve->last=e;  
}
