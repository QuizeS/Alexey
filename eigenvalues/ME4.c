#include <stdio.h>
#include <math.h>
#include <complex.h> 
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <unistd.h>
#include <poll.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

enum ecodes
{
  SUCCESS = 0,
  ERR_NOCONF,
  ERR_LUA,
  ERR_INCORRECT_VALUE,
  ERR_MODEL_NAME,
  ERR_TOO_LONG,
  ERR_WRONG_VAL
};

enum aux_param
{
  FLAVS         =    3,
  PAR_MAX_WIDTH =    6,
  MAX_LEN       =  255,
  MAX_READ_LEN  = 1024,
  MAX_BUF_SIZE  = 4096
};

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

enum
{
  P_A=0,P_B,P_E,P_TOL,
  P_S12,P_S13,P_PSI0,
  P_H,P_OUT,P_MODEL,P_SCAN_STEPS,P_TOTAL
};

int NUM=100;
double cA=500.,cB=50.;
const char *pref = "#"; 
char prog_name[MAX_LEN + 1];

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
void prtparam (param*,FILE*);
bool get_params(char *,param *,bool );//
bool chkin_poll(void);
bool process_cmd_line(int , char **, param *);
double vS (double);
double pLinear (double);
void lambda (double [FL], double, double);
void lambda2 (double [2*FL-1], double , double );
void prt_matr(double [FL][FL], char *,FILE*);
void prt_Cmatr(double complex [FL][FL], char *,FILE*);
void prt_Cvect(double complex [FL], char *,FILE*);
void ME4(const basic_ctx *,variative_ctx *); 
void eigval(const basic_ctx *,variative_ctx *,param *);

int main (int argc, char* argv[])
{
  int i,j;
  double complex sqPsi;
  double s12, s13, c12, c13, E, PSurv, q1, q2;
  param params[P_TOTAL];
  FILE *out;
  
  {
    char *slash;
    slash = strrchr(argv[0], '/');

    if( NULL == slash )
    {
      strncpy(prog_name, argv[0], MAX_LEN);
    }
    else
    {
      strncpy(prog_name, slash + 1, MAX_LEN);
    }
    prog_name[MAX_LEN] = '\0';
  }
  
  initparam(params);
  
  bool psi0_changed = false;
  
  psi0_changed = process_cmd_line(argc, argv, params);
  
  if(0 == strncmp(params[P_OUT].val.n,"SCREEN",FNAME))
  {
     out = stdout;
  }
  else
  {
    out = fopen(params[P_OUT].val.n,"w+");
  }
  prtparam(params,out);
  
  NUM= (int) params[P_SCAN_STEPS].val.v; 
  E=params[P_E].val.v;
  PSurv=0.;
  sqPsi=0.+I*0.;
  s12=params[P_S12].val.v; s13=params[P_S13].val.v;
  c12=sqrt(1. - s12*s12); c13=sqrt(1. - s13*s13);
  q1=4.35196e6; q2 =0.030554;
  
  basic_ctx basic;
  variative_ctx vartve;

  basic.H0[0][0] = 0.; basic.H0[0][1] = 0.;        basic.H0[0][2] = 0.;
  basic.H0[1][0] = 0.; basic.H0[1][1] = (q1*q2)/E; basic.H0[1][2] = 0.;
  basic.H0[2][0] = 0.; basic.H0[2][1] = 0.;        basic.H0[2][2] = q1/E;

  basic.tol = params[P_TOL].val.v;//e-4
  basic.a = params[P_A].val.v;//0.1
  basic.b = params[P_B].val.v;//1
  
  basic.Psi0[0] = params[P_PSI0].val.c[0]; 
  basic.Psi0[1] = params[P_PSI0].val.c[1]; 
  basic.Psi0[2] = params[P_PSI0].val.c[2];

  basic.W[0][0] = c13*c13*c12*c12; basic.W[0][1] = c12*s12*c13*c13; basic.W[0][2] = c12*c13*s13;
  basic.W[1][0] = basic.W[0][1];   basic.W[1][1] = s12*s12*c13*c13; basic.W[1][2] = s12*c13*s13;
  basic.W[2][0] = basic.W[0][2];   basic.W[2][1] = basic.W[1][2];    basic.W[2][2] = s13*s13;

  //basic.vS = vS;
  basic.vS = pLinear;

  vartve.Psi[0]=0.+I*0.;
  vartve.Psi[1]=0.+I*0.; 
  vartve.Psi[2]=0.+I*0.;
  vartve.dh=params[P_H].val.v;
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
  prt_matr(basic.H0,"H0",out);
  prt_matr(basic.W,"W",out);
  prt_matr(basic.H0W,"H0W",out);
  prt_matr(basic.H0H0W,"H0H0W",out);
  prt_matr(basic.WH0W,"WH0W",out);


  eigval(&basic,&vartve,params);
  
  fprintf(out,"#################################################\n");
  fprintf(out,"##            CALCULATION COMPLETED            ##\n");
  fprintf(out,"#################################################\n");

  ME4(&basic,&vartve);
  
  PSurv = c12*c12*c13*c13*(vartve.Psi[0]*conj(vartve.Psi[0]))
         +s12*s12*c13*c13*(vartve.Psi[1]*conj(vartve.Psi[1]))
         +s13*s13*(vartve.Psi[2]*conj(vartve.Psi[2]));
  
  for(i=0;i<FL;i++)
  {
    sqPsi += vartve.Psi[i]*conj(vartve.Psi[i]);
  }
  fprintf(out,"# Psi[0]^2 = %lf + I%lf\n",
    creal(vartve.Psi[0]*conj(vartve.Psi[0])),cimag(vartve.Psi[0]*conj(vartve.Psi[0])));
  fprintf(out,"# Psi[1]^2 = %lf + I%lf\n",
    creal(vartve.Psi[1]*conj(vartve.Psi[1])),cimag(vartve.Psi[1]*conj(vartve.Psi[1])));
  fprintf(out,"# Psi[2]^2 = %lf + I%lf\n",
    creal((vartve.Psi[2])*conj(vartve.Psi[2])),cimag(vartve.Psi[2]*conj(vartve.Psi[2])));

  fprintf(out,"# phi0 = %lf \n",atan2(creal(vartve.Psi[0]),cimag(vartve.Psi[0])));
  fprintf(out,"# phi1 = %lf \n",atan2(creal(vartve.Psi[1]),cimag(vartve.Psi[1])));
  fprintf(out,"# phi2 = %lf \n",atan2(creal(vartve.Psi[2]),cimag(vartve.Psi[2])));

  prt_Cvect(vartve.Psi,"Psi",out);
  fprintf(out,"#ProbSurv b 1-|Psi|^2 calls E tol dh dh(prev) dh(last)\n");
  fprintf(out,"%14.12g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\n",
    PSurv,vartve.last,1.-creal(sqPsi),vartve.calls,E,params[P_TOL].val.v,
    params[P_H].val.v,vartve.prev,vartve.dh);

  if(stdout != out)
  {
    fclose(out);
  }

return 0;
}
bool process_cmd_line(int argc, char **argv, param *params)
{
  char aux_store[MAX_READ_LEN + 1];
  int64_t m_len;
  FILE *fdata;
  bool psi0_changed = false;

  m_len = (int) snprintf(aux_store, MAX_LEN, "%s.lua", prog_name);
  if(m_len > MAX_LEN)
  {
    fprintf(stderr, "[%s] ERROR: default config file name '%s.lua' is too long.\n",
            prog_name, prog_name);

    return ERR_TOO_LONG;
  }
  fdata = fopen(aux_store, "r");
  if(fdata == NULL)
  {
    fprintf(stderr, "[%s] WARNING: default config file '%s' doesn't exist ",
            prog_name, aux_store);
    fprintf(stderr, "or you don't have read permission for it.\n");
  }
  else
  {
    psi0_changed = get_params(aux_store, params, true);
  }

  if(chkin_poll())
  {
    fgets(aux_store, MAX_READ_LEN, stdin);
    psi0_changed = get_params(aux_store, params, false);
  }

  int k = 1;
  while(k < argc)
  {
    if(strcmp(argv[k], "-c") == 0)
    {
      if(k + 1 < argc)
      {
        strncpy(aux_store, argv[k+1], MAX_READ_LEN);
        aux_store[MAX_READ_LEN] = '\0';

        psi0_changed = get_params(aux_store, params, false);
        k++;
        k++;
        continue;
      }
      else
      {
        fprintf(stderr, "[%s] WARNING: the option '-c' expects string after it.",
                prog_name);
        k++;
        continue;
      }
    }
    else
    {
      fdata = fopen(argv[k], "r");
      if(fdata == NULL)
      {
        fprintf(stderr, "[%s] WARNING: config file '%s' doesn't exist or",
                prog_name,aux_store);
        fprintf(stderr, "you don't have read permission for it.\n");
      }
      else
      {
        fclose(fdata);

        m_len = snprintf(aux_store, MAX_LEN, "%s", argv[k]);
        if(m_len > MAX_LEN)
        {
          fprintf(stderr, "[%s] ERROR: name for config file '%*s' ",
                  prog_name, MAX_LEN, argv[k]);
          fprintf(stderr, "is too long.\n");

          return ERR_TOO_LONG;
        }

        psi0_changed = get_params(aux_store, params, true);
        k++;
        continue;
      }
    }
  }
  return psi0_changed;
}
bool get_params(char *src, param *params, bool file_src)
{
  int lstat;
  char aux_store[MAX_READ_LEN];
  int64_t m_len;
  bool psi = false;

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

  if(file_src)
  {
    lstat = luaL_loadfile(L, src);
    if(lstat != LUA_OK)
    {
      const char *mess = lua_tostring(L, -1);
      fprintf(stderr, "[%s] ERROR: %s\n", prog_name, mess);
      lua_pop(L, 1);
      lua_close(L);

      exit(ERR_LUA);
    }
  }
  else
  {
    lstat = luaL_loadstring(L, src);
    if(lstat != LUA_OK)
    {
      const char *mess = lua_tostring(L, -1);
      fprintf(stderr, "[%s] ERROR: %s\n", prog_name, mess);
      lua_pop(L, 1);
      lua_close(L);

      exit(ERR_LUA);
    }
  }

  lstat = lua_pcall(L, 0, LUA_MULTRET, 0);
  if(lstat != LUA_OK)
  {
    const char *mess = lua_tostring(L, -1);
    fprintf(stderr, "[%s] ERROR: %s\n", prog_name, mess);
    lua_pop(L,1);
    lua_close(L);

    exit(ERR_LUA);
  }

  for(uint8_t j1 = 0; j1 < P_TOTAL; j1++)
  {
    memset(aux_store, '\0', MAX_READ_LEN);
    lstat = lua_getglobal(L, params[j1].name);

    if(lstat == LUA_TNONE || lstat == LUA_TNIL)
    {
      continue;
    }

    char *mm;
    double re, im;
    int isnum;

    switch(lstat)
    {
    case LUA_TNUMBER:
      params[j1].val.v = lua_tonumber(L, -1);
      break;
    case LUA_TSTRING:
      mm = (char *) lua_tolstring(L, -1, (size_t *) &m_len);
      if(mm == NULL)
      {
        /// error, either print it or do something other.
      }
      if(m_len > MAX_LEN)
      {
        fprintf(stderr, "[%s] ERROR: option value '%s' is too long\n",
                prog_name, mm);

        lua_close(L);

        exit(ERR_TOO_LONG);
      }
      else
      {
        m_len = snprintf(params[j1].val.n, MAX_LEN, "%s", mm);
        if(m_len > MAX_LEN)
        {
          fprintf(stderr, "[%s] ERROR: parameter '%s' value '%*s' is too long\n",
                  prog_name, params[j1].name, MAX_LEN, mm);
          lua_close(L);

          exit(ERR_TOO_LONG);
        }
      }
      break;
    case LUA_TTABLE:
      for(uint8_t j2 = 0; j2 < FLAVS; j2++)
      {
        lua_pushinteger(L, j2 + 1);
        lua_gettable(L, -2);

        int tl = lua_type(L, -1);
        if(tl != LUA_TTABLE)
        {
          fprintf(stderr, "[%s] ERROR: parameter '%s' has wrong value. ",
                  prog_name, params[j1].name);
          fprintf(stderr, "Expected type: '%s', obtained: '%s'\n",
                  lua_typename(L, LUA_TTABLE),
                  lua_typename(L, lua_type(L, -1)));

          lua_close(L);

          exit(ERR_WRONG_VAL);
        }

        lua_pushinteger(L, 1);
        lua_gettable(L, -2);
        re = lua_tonumberx(L, -1, &isnum);
        if(!isnum)
        {
          fprintf(stderr, "[%s] ERROR: wrong value for parameter '%s'.\n",
                  prog_name, params[j1].name);

          lua_close(L);

          exit(ERR_WRONG_VAL);
        }
        lua_pop(L,1);

        lua_pushinteger(L, 2);
        lua_gettable(L, -2);
        im = lua_tonumberx(L, -1, &isnum);
        if(!isnum)
        {
          fprintf(stderr, "[%s] ERROR: wrong value for parameter '%s'.\n",
                  prog_name, params[j1].name);

          lua_close(L);

          exit(ERR_WRONG_VAL);
        }
        lua_pop(L,2);

        params[j1].val.c[j2] = re + I*im;
        psi = true;
      }
      break;
    default:
      /// Unsupported
      break;
    }
  }

  lua_close(L);
  return psi;
}
bool chkin_poll(void)
{
  int ret;
  struct pollfd pfd[1] = {0};

  pfd[0].fd = STDIN_FILENO;
  pfd[0].events = POLLIN;
  ret = poll(pfd, 1, 0);

  return (ret > 0);
}
void initparam (param *par)
{
  double c12,c13,s12,s13;
  par[P_A].type = FLOAT; 
  strncpy(par[P_A].name,"a",PARAM_NAME-1);   
  par[P_A].val.v = 0.1;
    
  par[P_E].type = FLOAT;
  strncpy(par[P_E].name,"E",PARAM_NAME-1); 
  par[P_E].val.v = 1.0;
  
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

  par[P_SCAN_STEPS].type = FLOAT;
  strncpy(par[P_SCAN_STEPS].name,"scan_steps",PARAM_NAME-1); 
  par[P_SCAN_STEPS].val.v = 100.; 
  
  par[P_OUT].type = STRING;
  strncpy(par[P_OUT].name,"out",PARAM_NAME-1); 
  strncpy(par[P_OUT].val.n,"SCREEN",FNAME-1);
  
  par[P_MODEL].type = STRING;
  strncpy(par[P_MODEL].name,"model",PARAM_NAME-1); 
  strncpy(par[P_MODEL].val.n,"sun",FNAME-1);
}
void prtparam(param *par, FILE *out)
{ 
  for(int i=0;i<P_TOTAL;i++)
  { 
    if(par[i].type == FLOAT)
      fprintf(out,"%s %s = %lf\n",pref,par[i].name,par[i].val.v);
    if(par[i].type == COMPLEX)
      fprintf(out,"%s %s = %lf + I%lf\t%lf + I%lf\t%lf + I%lf\n",
        pref,par[i].name,
        creal(par[i].val.c[0]), cimag(par[i].val.c[0]),
        creal(par[i].val.c[1]), cimag(par[i].val.c[1]),
        creal(par[i].val.c[2]), cimag(par[i].val.c[2]));
        
    if(par[i].type == STRING)
      fprintf(out,"%s %s = %s\n",pref,par[i].name,par[i].val.n);  
  }
}
double vS (double e)
{
  return 65956.*exp(-10.54*e);
}
double pLinear (double e)
{
  return cA+e*(cB-cA);
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
void prt_matr(double matr[FL][FL], char *name,FILE *out)
{
  int i,j;
  fprintf(out,"%s %s = ",pref,name);
  for(i=0;i<FL;i++)
  {
    fprintf(out,"\n%s\t",pref);
    for(j=0;j<FL;j++)
    {
      fprintf(out," %e",matr[i][j]);
    }
  }
  fprintf(out,"\n");
}
void prt_Cmatr(complex double matr[FL][FL], char *name,FILE *out)
{
  int i,j;
  fprintf(out,"\n#%s",name);
  for(i=0;i<FL;i++)
  {
    fprintf(out,"\n%s\t",pref);		
    for(j=0;j<FL;j++)
    {
      fprintf(out,"%lf+i(%lf)",creal(matr[i][j]),cimag(matr[i][j]));
    }
  }
}
void prt_Cvect(complex double vect[FL], char *name,FILE *out)
{
  int i;
  fprintf(out,"%s %s = \n",pref,name);
  for(i=0;i<FL;i++)
  {
    fprintf(out,"%s\t%e+i(%e)\n",pref,creal(vect[i]),cimag(vect[i]));
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
          +1./(var*var)*r1*sqA[i][j])
          *(cos(var*l[0]*vartve->dh)+I*sin(var*l[0]*vartve->dh));
      }
    }
    //уравнение эволюции нейтрино в точке en
    for(i=0;i<FL;i++)
    {
      vartve->Psi[i] = 
        (exA[i][0]*Psi[0]
        +exA[i][1]*Psi[1]
        +exA[i][2]*Psi[2])*(cos(z*vartve->dh)+I*sin(z*vartve->dh));
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
#ifdef DEBUG
      fprintf(stderr,"#!!!Изменение шага dh = %g\n",vartve->dh);
#endif
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
      if((e+2.*vartve->dh)>(basic->b) && index != LAST && index != OVER)
      {
        index = PREV;
        vartve->prev = vartve->dh;
        vartve->dh = (basic->b-e)/2.;
#ifdef DEBUG
        fprintf(stderr,"#!!! if Изменение шага dh = %lf!!!\n",vartve->dh);
#endif
      }
    }
    vartve->calls++;
  }
  vartve->last=e;  
}
void eigval(const basic_ctx *basic, variative_ctx *vartve, param *pr)
{
    int i,j;
    char aux_store[MAX_READ_LEN + 1];
    double t,p,q; 
    double sqA[FL][FL],A[FL][FL],l[2*FL-1];
    FILE *file;
    snprintf(aux_store, MAX_LEN, "eigen_%d_E%5.3e.dat",NUM,pr[P_E].val.v);
    file = fopen(aux_store,"w+");


    basic_ctx tb;
    memcpy(tb.H0,basic->H0,FL*FL*sizeof(double)); 
    memcpy(tb.W,basic->W,FL*FL*sizeof(double));
    memcpy(tb.H0W,basic->H0W,FL*FL*sizeof(double));
    memcpy(tb.H0H0W,basic->H0H0W,FL*FL*sizeof(double));
    memcpy(tb.WH0W,basic->WH0W,FL*FL*sizeof(double));
    memcpy(tb.Psi0,basic->Psi0,FL*sizeof(complex double));
    tb.tol=basic->tol; tb.a=basic->a; tb.vS=basic->vS;

    t=basic->a+(basic->b-basic->a)/NUM;
    while(t < basic->b+(basic->b-basic->a)/(2.*NUM))
    {
      for(i=0;i<FL;i++)
      {
        for(j=0;j<FL;j++)
        {
          A[i][j] = basic->H0[i][j] + basic->vS(t)*basic->W[i][j];
        }
      }
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

      p = (double)((1./2.)*(sqA[0][0]+sqA[1][1]+sqA[2][2])
            -(5./6.)*(A[0][0]+A[1][1]+A[2][2])*(A[0][0]+A[1][1]+A[2][2]));
      q = (double)(
        (-1.)*(A[0][0]*A[1][1]*A[2][2]+A[1][0]*A[0][2]*A[2][1]     
        +A[2][0]*A[0][1]*A[1][2]-A[2][0]*A[1][1]*A[0][2]
        -A[0][0]*A[1][2]*A[2][1]-A[1][0]*A[0][1]*A[2][2])
        -(13./54.)*(A[0][0]+A[1][1]+A[2][2])
        *(A[0][0]+A[1][1]+A[2][2])
        *(A[0][0]+A[1][1]+A[2][2])
        +(1./6.)*(A[0][0]+A[1][1]+A[2][2])*(sqA[0][0]+sqA[1][1]+sqA[2][2]));
     
      lambda2(l,q,p);
      tb.b=t;
#ifdef DEBUG2
      fprintf(stderr,"%lf\n",t);
#endif
      ME4(&tb,vartve);
      vartve->dh=1e-3;
      fprintf(file,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                      t,l[0],l[1],l[2],l[3],l[4],
                      atan2(creal(vartve->Psi[0]),cimag(vartve->Psi[0])),
                      atan2(creal(vartve->Psi[1]),cimag(vartve->Psi[1])),
                      atan2(creal(vartve->Psi[2]),cimag(vartve->Psi[2]))
                      );

      t = t+(basic->b-basic->a)/NUM;
    }
    fclose(file);
}
void lambda2 (double l[2*FL-1], double q, double p)
{
  int i,j;
  double var;
  l[2] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(-3./p))-(2.*M_PI*0.)/3.);
  l[1] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(-3./p))-(2.*M_PI*1.)/3.);
  l[0] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(-3./p))-(2.*M_PI*2.)/3.);
  
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
  l[4]=l[2]-l[0];
  l[3]=l[1]-l[0];
  l[0]=sqrt(-p/3.)*l[0];
  l[1]=sqrt(-p/3.)*l[1];
  l[2]=sqrt(-p/3.)*l[2];
  l[3]=sqrt(-p/3.)*l[3];
  l[4]=sqrt(-p/3.)*l[4];

}

