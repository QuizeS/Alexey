#include <stdio.h>
#include <math.h>
#include <complex.h> 

double vSun (double , double , double );

enum{
	FLAVS=3
};

double vSun (double g, double eta, double e){
	
double v;
		
		v = g*exp((-1.)*eta*e);
		
return v;
}
int main (){
	
	int i,j,k;
	double  a=4.35196e6, b=0.030554, E=7.05, g=6.5956e4, eta=10.54;	//a=4.35196*10^6 умножить на 10^6
	double  s12=0.308, s13=0.0234;
	double  c12 , c13;
	double  dh=0.01, ep=0.,em=0.,e=0. , e0 = 0.01;
	double  complex z=0.+0.*I, p=0.+0.*I, r0=0.+0.*I, r1=0.+0.*I,YY=0.+I*0.;//z-след ,p - след от квадрата матрицы
	double  l0,l1,l2,q,a1,b1; //корни характеристического многочлена
	double 	var;
	double  H0[3][3],W[3][3],H0W[3][3],H0H0W[3][3],WH0W[3][3],arr[3][3];
	double complex Y[3],Y0[3], OMG4[3][3], OMG04[3][3],sqOMG04[3][3],unit[3][3],expOMG04[3][3];
	
	H0[0][0] = 0.; H0[0][1] = 0.; H0[0][2] = 0.;//матрица H0
	H0[1][0] = 0.; H0[1][1] = (b*a)/E; H0[1][2] = 0.;
	H0[2][0] = 0.; H0[2][1] = 0.; H0[2][2] = (a*(1./E));
	
	var = 1. - s12*s12;						//значение cos-ов углов тэта 12 и 13
	c12 = sqrt(var);
	var = 1. - s13*s13;
	c13 = sqrt(var);

	Y0[0]=c12*c13; Y0[1]=s12*c13; Y0[2]=s13;

	W[0][0] = c13*c13*c12*c12; W[0][1] = c12*s12*c13*c13; W[0][2] = c12*c13*s13;
	W[1][0]= W[0][1]; W[1][1] = s12*s12*c13*c13; W[1][2] = s12*c13*s13;
	W[2][0] = W[0][2]; W[2][1] = W[1][2]; W[2][2] = s13*s13;

    H0W[0][0]=0.; H0W[0][1]=0.; H0W[0][2]=0.;  // [H0,W]
	H0W[1][0]=0.; H0W[1][1]=0.; H0W[1][2]=0.;
	H0W[2][0]=0.; H0W[2][1]=0.; H0W[2][2]=0.;
	
	WH0W[0][0]=0.; WH0W[0][1]=0.; WH0W[0][2]=0.;  // [W,[H0,W]]
	WH0W[1][0]=0.; WH0W[1][1]=0.; WH0W[1][2]=0.;
	WH0W[2][0]=0.; WH0W[2][1]=0.; WH0W[2][2]=0.;
	
	H0H0W[0][0]=0.; H0H0W[0][1]=0.; H0H0W[0][2]=0.;  // [H0,[H0,W]
	H0H0W[1][0]=0.; H0H0W[1][1]=0.; H0H0W[1][2]=0.;
	H0H0W[2][0]=0.; H0H0W[2][1]=0.; H0H0W[2][2]=0.;
	
	unit[0][0]=1.; unit[0][1]=0.; unit[0][2]=0.; //I в алгоритме вычисления exp от матрицы
	unit[1][0]=0.; unit[1][1]=1.; unit[1][2]=0.;
	unit[2][0]=0.; unit[2][1]=0.; unit[2][2]=1.;
	
	sqOMG04[0][0]=0.; sqOMG04[0][1]=0.; sqOMG04[0][2]=0.; // Квадрат A_0 для вычисления Tr(A_0)^2
	sqOMG04[1][0]=0.; sqOMG04[1][1]=0.; sqOMG04[1][2]=0.;
	sqOMG04[2][0]=0.; sqOMG04[2][1]=0.; sqOMG04[2][2]=0.;
	
		for(i=0;i<3;i++){					    //Вычисление [H0,W]
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					arr[i][j] = H0[i][k]*W[k][j]-W[i][k]*H0[k][j];
					H0W[i][j] = H0W[i][j]+arr[i][j];// var
				}
			}
		}     									//
	
	
		for(i=0;i<3;i++){					    //Вычисление [H0,[H0,W]]
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					arr[i][j] = H0[i][k]*H0W[k][j]-H0W[i][k]*H0[k][j];
					H0H0W[i][j] = H0H0W[i][j]+arr[i][j];
				}
			}
		}										//
	
	
		for(i=0;i<3;i++){						//Вычисление [W,[H0,W]]
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					arr[i][j] = W[i][k]*H0W[k][j]-H0W[i][k]*W[k][j];
					WH0W[i][j] = WH0W[i][j]+arr[i][j];
				}
			}
		}										//	
	
		printf("WH0W\n");
		for(i=0;i<3;i++){						//Печать матриц
			for(j=0;j<3;j++){
				printf("%lf  ",WH0W[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		
	    printf("H0H0W\n");	
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				printf("%lf  ",H0H0W[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	
		printf("H0W\n");
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				printf("%lf  ",H0W[i][j]);
			}
			printf("\n");
		}
		
		
												//vSun -> var
		e = e0;	
											
		//while(e<1.){
			
			ep = e+(1.+1./sqrt(3.))*(dh/2.);
			em = e+(1.-1./sqrt(3.))*(dh/2.);;
			
			for(i=0;i<3;i++){						  //вычисление матрицы OMG4  в  точке e_n(кси_n)
				for(j=0;j<3;j++){
					OMG4[i][j] = ((-1.)*(H0[i][j]+0.5*(vSun(g,eta,ep)+vSun(g,eta,em))*W[i][j])*dh*I+(sqrt(3.)/12.)*(vSun(g,eta,ep)-vSun(g,eta,em))*H0W[i][j]*dh*dh)*(-1.*I);
				}
			}
			
			z=(OMG4[0][0]+OMG4[1][1]+OMG4[2][2])/3.; //след OMG4
					
			for(i=0;i<3;i++){		 				 //Tr(OMG04)=0
				for(j=0;j<3;j++){
					OMG04[i][j]=OMG4[i][j]-z*unit[i][j];
				}
			}										 //	
			
			for(i=0;i<3;i++){					     //Вычисление sqOMG04=OMG04^2
				for(j=0;j<3;j++){
					p=0.+0.*I;
					for(k=0;k<3;k++){
						p += OMG04[i][k]*OMG04[k][j];
					}
					sqOMG04[i][j]=p;
				}
			}										 //
			
			p =(sqOMG04[0][0]+sqOMG04[1][1]+sqOMG04[2][2])/2.; //cлед sqOMG04^2
			
			q = OMG04[0][0]*OMG04[1][1]*OMG04[2][2]+OMG04[1][0]*OMG04[0][2]*OMG04[2][1]+OMG04[2][0]*OMG04[0][1]*OMG04[1][2]-OMG04[2][0]*OMG04[1][1]*OMG04[0][2]-OMG04[0][0]*OMG04[1][2]*OMG04[2][1]-OMG04[1][0]*OMG04[0][1]*OMG04[2][2];
			
			printf("q = < %lf > , p = %lf + i*%lf \n",q,creal(p),cimag(p));
			
			l2 = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*0.)/3.);//корни характер.многочлена
			l1 = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*1.)/3.);
			l0 = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*2.)/3.);
			a1 = l1-l0;
			b1 = l2-l0;
			var=sqrt(p/3.);
			printf("<l0 = %lf> <l1 = %lf>  <l2 = %lf>, <qp/3..>=%lf \n",l0,l1,l2,q/(var*var*var));
			
			
			
			r0 = (-1.)*(2.*sin((var*a1*e)/2.)*sin((var*a1*e)/2.)
				 -I*sin(var*a1*e))/a1;//r0 , r1 для exp(..A0)
				 
			r1 = (-1./(a1-b1))*((2.*sin((var*a1*e)/2.)*sin((var*a1*e)/2.)
			     -I*sin(var*a1*e))/a1
			 	 -(2.*sin((var*b1*e)/2.)*sin((var*b1*e)/2.)
				 -I*sin(var*b1*e))/b1
				 );
			
			for(i=0;i<3;i++){
				for(j=0;j<3;j++){
					expOMG04[i][j] = exp(sqrt(p/3.)*l0*I*e)*(
					(1.-l0*(r0-l1*r1))*unit[i][j]
					+sqrt(3./p)*(r0+l2*r1)*OMG04[i][j]
					+(3./p)*r1*sqOMG04[i][j]
					);
				}
			}
			
			for(i=0;i<3;i++){
				for(j=0;j<3;j++){
					Y[i] += expOMG04[i][j]*Y0[j];
				}
				Y[i] *= exp(e*z*I);
			}
			for(i=0;i<3;i++){
					YY += Y[i]*conj(Y[i]);
			}
			printf("<YY* = %lf + i*%lf>\n",creal(YY),cimag(YY));
			
				//e = e+dh;
		//}
			
			
			unit[0][0]=0.; unit[0][1]=0.; unit[0][2]=0.; 
			unit[1][0]=0.; unit[1][1]=0.; unit[1][2]=0.;
			unit[2][0]=0.; unit[2][1]=0.; unit[2][2]=0.;
			
			for(i=0;i<3;i++){					    
				for(j=0;j<3;j++){
					p=0.+I*0.;
					for(k=0;k<3;k++){
						p += expOMG04[i][k]*conj(expOMG04[j][k]);
					}
					unit[i][j]=p;
				}
			}     						
		
		printf("expOMG04expOMG04*\n");
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				printf("%lf + (%lf)*i  ",creal(unit[i][j]),cimag(unit[i][j]));
			}
			printf("\n");
		}
		
	
	
return 0;	
}


/*	
	H0[0][0] = 1.; H0[0][1] = 0.; H0[0][2] = 0.;//матрица H0
	H0[1][0] = 0.; H0[1][1] = 1.; H0[1][2] = 0.;
	H0[2][0] = 0.; H0[2][1] = 0.; H0[2][2] = 1.;

	W[0][0] = 1.; W[0][1] = 0.; W[0][2] = 0.;
	W[1][0]= W[0][1]; W[1][1] = 1.; W[1][2] = 0.;
	W[2][0] = W[0][2]; W[2][1] = W[1][2]; W[2][2] = 1.;
*/
