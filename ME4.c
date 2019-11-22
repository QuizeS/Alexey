#include <stdio.h>
#include <math.h>
#include <complex.h> 
#include <stdlib.h>
#include <string.h>

enum
{
  FL=3,
  NORM,
  PREV,
  LAST,
  OVER,
  PARAM_NAME=16,
  FNAME=256
};

enum  tag_type 
{
  FLOAT,
  COMPLEX,
  STRING
};
//printf->fprintf
enum
{
  P_A=0,P_B,P_E,P_TOL,
  P_S12,P_S13,P_PSI0,
  P_H,P_OUT,P_MODEL,P_TOTAL
};

const char *pref = "#"; 

typedef double (*dens)(double);

typedef struct
{
  enum tag_type type;
  char name[PARAM_NAME];
  union
  {
    double v;
    complex double c[FL];
    char n[FNAME];
  } val;

} param;

typedef struct 
{
  double complex Psi0[FL];
  double H0[FL][FL], W[FL][FL], H0W[FL][FL], WH0W[FL][FL], H0H0W[FL][FL];
  double tol, a, b;
  dens vS;
  
} basic_ctx;

typedef struct
{
  double complex Psi[FL];
  double Er, dh, last, prev;
  int calls;
  
} variative_ctx;

void initparam (param*);
void prtparam (param*);
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
  param params[P_TOTAL];
  
  initparam(params);
  prtparam(params);
  
  E=atof(argv[1]);
  PSurv=0.;
  sqPsi=0.+I*0.;
  s12=sqrt(0.308); s13=sqrt(0.0234);
  c12=sqrt(1. - 0.308); c13=sqrt(1. - 0.0234);
  q1=4.35196e6; q2 =0.030554;
  
  basic_ctx basic;
  variative_ctx vartve;

  basic.H0[0][0] = 0.; basic.H0[0][1] = 0.;        basic.H0[0][2] = 0.;
  basic.H0[1][0] = 0.; basic.H0[1][1] = (q1*q2)/E; basic.H0[1][2] = 0.;
  basic.H0[2][0] = 0.; basic.H0[2][1] = 0.;        basic.H0[2][2] = q1/E;

  basic.tol = atof(argv[2]);//e-4
  basic.a = atof(argv[3]);//0.1
  basic.b = atof(argv[4]);//1
  
  basic.Psi0[0]=c12*c13+I*0.;
  basic.Psi0[1]=s12*c13+I*0.; 
  basic.Psi0[2]=s13+I*0.;

  basic.W[0][0] = c13*c13*c12*c12; basic.W[0][1] = c12*s12*c13*c13; basic.W[0][2] = c12*c13*s13;
  basic.W[1][0] = basic.W[0][1];   basic.W[1][1] = s12*s12*c13*c13; basic.W[1][2] = s12*c13*s13;
  basic.W[2][0] = basic.W[0][2];   basic.W[2][1] = basic.W[1][2];    basic.W[2][2] = s13*s13;

  basic.vS = vS;

  vartve.Psi[0]=0.+I*0.;
  vartve.Psi[1]=0.+I*0.; 
  vartve.Psi[2]=0.+I*0.;
  vartve.dh=(basic.tol)/10.;
  vartve.Er=0.;
  vartve.calls=0;
  vartve.prev=0.;

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
  printf("%s a = %e\n%s b = %e\n%s tol = %e\n%s E = %e\n%s dh = %e\n",
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
  printf("\n");
  
  ME4(&basic,&vartve);
  
  PSurv = c12*c12*c13*c13*(vartve.Psi[0]*conj(vartve.Psi[0]))
         +s12*s12*c13*c13*(vartve.Psi[1]*conj(vartve.Psi[1]))
         +s13*s13*(vartve.Psi[2]*conj(vartve.Psi[2]));
  
  for(i=0;i<FL;i++)
  {
    sqPsi += vartve.Psi[i]*conj(vartve.Psi[i]);
  }
  printf("# Psi[0]^2 = %lf + I%lf\n",
    creal(vartve.Psi[0]*conj(vartve.Psi[0])),cimag(vartve.Psi[0]*conj(vartve.Psi[0])));
  printf("# Psi[1]^2 = %lf + I%lf\n",
    creal(vartve.Psi[1]*conj(vartve.Psi[1])),cimag(vartve.Psi[1]*conj(vartve.Psi[1])));
  printf("# Psi[2]^2 = %lf + I%lf\n",
    creal((vartve.Psi[2])*conj(vartve.Psi[2])),cimag(vartve.Psi[2]*conj(vartve.Psi[2])));
  
  printf("\n");
  
  printf("# phi0 = %lf \n",atan2(creal(vartve.Psi[0]),cimag(vartve.Psi[0])));
  printf("# phi1 = %lf \n",atan2(creal(vartve.Psi[1]),cimag(vartve.Psi[1])));
  printf("# phi2 = %lf \n",atan2(creal(vartve.Psi[2]),cimag(vartve.Psi[2])));
  
  printf("\n");
  
  prt_Cvect(vartve.Psi,"Psi");
  printf("#ProbSurv b 1-|Psi|^2 dh(last) dh(prev) calls\n");
  printf("%10.8g\t%g\t%g\t%g\t%g\t%d\n",
    PSurv,vartve.last,1.-creal(sqPsi),vartve.dh,vartve.prev,vartve.calls);


return 0;
}

void initparam (param *par)
{
  double c12,c13,s12,s13;
  par[P_A].type = FLOAT; 
  strncpy(par[P_A].name,"a",PARAM_NAME-1);   
  par[P_A].val.v = 0.1;
    
  par[P_E].type = FLOAT;
  strncpy(par[P_E].name,"E",PARAM_NAME-1); 
  par[P_E].val.v = 7.05;
  
  par[P_B].type = FLOAT; 
  strncpy(par[P_B].name,"b",PARAM_NAME-1);    
  par[P_B].val.v = 1.0;  
  
  par[P_TOL].type = FLOAT;
  strncpy(par[P_TOL].name,"tol",PARAM_NAME-1); 
  par[P_TOL].val.v = 1e-4;
  
  par[P_S12].type = FLOAT;        
  strncpy(par[P_S12].name,"s12",PARAM_NAME-1);         
  par[P_S12].val.v = sqrt(0.308); 
  
  par[P_S13].type = FLOAT;
  strncpy(par[P_S13].name,"s13",PARAM_NAME-1); 
  par[P_S13].val.v = sqrt(0.0234);
  
  par[P_PSI0].type = COMPLEX;
  strncpy(par[P_PSI0].name,"psi0",PARAM_NAME-1);
  s12=sqrt(0.308); s13=sqrt(0.0234);
  c12=sqrt(1. - 0.308); c13=sqrt(1. - 0.0234); 
  par[P_PSI0].val.c[0] = c12*c13+I*0.;
  par[P_PSI0].val.c[1] = s12*c13+I*0.; 
  par[P_PSI0].val.c[2] = s13+I*0.;
  
  par[P_H].type = FLOAT;
  strncpy(par[P_H].name,"dh",PARAM_NAME-1); 
  par[P_H].val.v = par[P_TOL].val.v;
  
  par[P_OUT].type = STRING;
  strncpy(par[P_OUT].name,"out",PARAM_NAME-1); 
  strncpy(par[P_OUT].val.n,"SCREEN",FNAME-1);
  
  par[P_MODEL].type = STRING;
  strncpy(par[P_MODEL].name,"model",PARAM_NAME-1); 
  strncpy(par[P_MODEL].val.n,"sun",FNAME-1);
}
void prtparam(param *par)
{ 
  for(int i=0;i<P_TOTAL;i++)
  { 
    if(par[i].type == FLOAT)
      printf("%s %s = %lf\n",pref,par[i].name,par[i].val.v);
    if(par[i].type == COMPLEX)
      printf("%s %s = %lf + I%lf\t%lf + I%lf\t%lf + I%lf\n",
        pref,par[i].name,
        creal(par[i].val.c[0]), cimag(par[i].val.c[0]),
        creal(par[i].val.c[1]), cimag(par[i].val.c[1]),
        creal(par[i].val.c[2]), cimag(par[i].val.c[2]));
        
    if(par[i].type == STRING)
      printf("%s %s = %s\n",pref,par[i].name,par[i].val.n);  
  }
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
      fprintf(stderr,"#!!!Изменение шага dh = %lf!!!\n",vartve->dh);
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
        vartve->prev = vartve->dh;
        vartve->dh = (1.-e)/2.;
        fprintf(stderr,"#!!! if Изменение шага dh = %lf!!!\n",vartve->dh);
      }
    }
    vartve->calls++;
  }
  vartve->last=e;  
}
