#include <stdio.h>
#include <math.h>
#include <complex.h> 

enum{
	FL=3
};

double vS (double , double , double );
double vS (double g, double eta, double e){	
  
 double v;
		
		v = g*exp((-1.)*eta*e);
		
 return v;
}


double ME4(double,double,double,double,double,double,
           double,double [FL][FL],double [FL][FL]); 
double ME4(double s12, double s13, double E, double dh, double tol,
           double g,double e0, double H0[FL][FL],double W[FL][FL]){
		
	int     i,j;
	double  q1=4.35196e6, q2=0.030554, eta=10.54;	//a=4.35196*10^6 умножить на 10^6
	double  c12 , c13;
	double  Er=0.,s=0.8, ep=0.,em=0.,e=0., vSep,vSem;
	double complex r0=0.+0.*I, r1=0.+0.*I,sqPsi=0.+I*0.;//z-след ,p - след от квадрата матрицы
	double  l[FL],q,p,z,a,b; 											//корни характеристического многочлена
	double 	var;
	double  H0W[FL][FL],H0H0W[FL][FL],WH0W[FL][FL];
	double complex Psi[FL],Psi0[FL], A[FL][FL], A0[FL][FL],sqA0[FL][FL],unit[FL][FL],exA0[FL][FL];
	double complex S1[FL][FL],S2[FL][FL],sqS1[FL][FL],arr[FL];

	H0[0][0] = 0.; H0[0][1] = 0.; H0[0][2] = 0.;      				//матрица H0
	H0[1][0] = 0.; H0[1][1] = (q1*q2)/E; H0[1][2] = 0.;
	H0[2][0] = 0.; H0[2][1] = 0.; H0[2][2] = (q1*(1./E));
	
	var = 1. - s12*s12;						         				 //значение cos-ов углов тэта 12 и 13
	c12 = sqrt(var);
	var = 1. - s13*s13;
	c13 = sqrt(var);

	arr[0]=0.; arr[1]=0.; arr[2]=0.;
	
	Psi0[0]=c12*c13; Psi0[1]=s12*c13; Psi0[2]=s13;

	W[0][0] = c13*c13*c12*c12; W[0][1] = c12*s12*c13*c13; W[0][2] = c12*c13*s13;
	W[1][0]= W[0][1]; W[1][1] = s12*s12*c13*c13; W[1][2] = s12*c13*s13;
	W[2][0] = W[0][2]; W[2][1] = W[1][2]; W[2][2] = s13*s13;
	
	unit[0][0]=1.; unit[0][1]=0.; unit[0][2]=0.; 					//I в алгоритме вычисления exp от матрицы
	unit[1][0]=0.; unit[1][1]=1.; unit[1][2]=0.;
	unit[2][0]=0.; unit[2][1]=0.; unit[2][2]=1.;
	
		for(i=0;i<FL;i++){					    					//Вычисление [H0,W]
			for(j=0;j<FL;j++){
				H0W[i][j] =  H0[i][0]*W[0][j]-W[i][0]*H0[0][j]
							+H0[i][1]*W[1][j]-W[i][1]*H0[1][j]
							+H0[i][2]*W[2][j]-W[i][2]*H0[2][j];
			}
		}     														//
		
		for(i=0;i<FL;i++){					    					//Вычисление [H0,[H0,W]]
			for(j=0;j<FL;j++){
				H0H0W[i][j] =  H0[i][0]*H0W[0][j]-H0W[i][0]*H0[0][j]
							  +H0[i][1]*H0W[1][j]-H0W[i][1]*H0[1][j]
							  +H0[i][2]*H0W[2][j]-H0W[i][2]*H0[2][j];	
			}
		}															//
	
	
		for(i=0;i<FL;i++){											//Вычисление [W,[H0,W]]
			for(j=0;j<FL;j++){
				WH0W[i][j] =  W[i][0]*H0W[0][j]-H0W[i][0]*W[0][j]
							 +W[i][1]*H0W[1][j]-H0W[i][1]*W[1][j]
							 +W[i][2]*H0W[2][j]-H0W[i][2]*W[2][j];
			}
		}															//	
	
		fprintf(stderr,"Коммутатор [W,[H0,W]]\n");							//Печать матриц
		for(i=0;i<FL;i++){											
			for(j=0;j<FL;j++){
				fprintf(stderr,"%lf  ",WH0W[i][j]);
			}
			fprintf(stderr,"\n");
		}
		fprintf(stderr,"\n");
		
	    fprintf(stderr,"Коммутатор [H0,[H0,W]]\n");	
		for(i=0;i<FL;i++){
			for(j=0;j<FL;j++){
				fprintf(stderr,"%lf  ",H0H0W[i][j]);
			}
			fprintf(stderr,"\n");
		}
		fprintf(stderr,"\n");
	
		fprintf(stderr,"Коммутатор [H0,W] \n");
		for(i=0;i<FL;i++){
			for(j=0;j<FL;j++){
				fprintf(stderr,"%lf  ",H0W[i][j]);
			}
			fprintf(stderr,"\n");
		}
		fprintf(stderr,"H0 \n");
		for(i=0;i<FL;i++){
			for(j=0;j<FL;j++){
				fprintf(stderr,"%lf  ",H0[i][j]);
			}
			fprintf(stderr,"\n");
		}
		fprintf(stderr,"W \n");
		for(i=0;i<FL;i++){
			for(j=0;j<FL;j++){
				fprintf(stderr,"%lf  ",W[i][j]);
			}
			fprintf(stderr,"\n");
		}
		fprintf(stderr,"\n");															//
		
																	
		e = e0;	
											
		while(e<1.){												//начало цикла по точкам
			ep = e+(1.+1./sqrt(3.))*(dh/2.);
			em = e+(1.-1./sqrt(3.))*(dh/2.);;
			
			vSem = vS(g,eta,em);
			vSep = vS(g,eta,ep);
			
			for(i=0;i<FL;i++){						  				//вычисление матрицы OMG4  в  точке e_n(кси_n)
				for(j=0;j<FL;j++){
					A[i][j] = ((-1.)*(H0[i][j]+0.5*(vSep+vSem)*W[i][j])*dh*I
							 +(sqrt(3.)/12.)*(vSep-vSem)*H0W[i][j]*dh*dh)*(-1.*I);
				}
			}										 				//
			
			//fprintf(stderr,"A \n");
			//for(i=0;i<FL;i++){
				//for(j=0;j<FL;j++){
					//fprintf(stderr,"%lf + i(%g)",creal(A[i][j]),cimag(A[i][j]));
				//}
				//fprintf(stderr,"\n");
			//}
			//fprintf(stderr,"\n");
				
			
			z = (double)(A[0][0]+A[1][1]+A[2][2])/3.; 				//след OMG4
					
			for(i=0;i<FL;i++){		 				 				//Tr(OMG04)=0	
				A0[i][i] = A[i][i]-z;
			}										 				//	
			
			for(i=0;i<FL;i++){					     				//Вычисление sqOMG04=OMG04^2
				for(j=0;j<FL;j++){
					sqA0[i][j] = A0[i][0]*A0[0][j]
								+A0[i][1]*A0[1][j]
								+A0[i][2]*A0[2][j];
				}
			}										 				//
			
			p = (double)(sqA0[0][0]+sqA0[1][1]+sqA0[2][2])/2.; 		//cлед sqA0^2
			
			q = (double)(
				 A0[0][0]*A0[1][1]*A0[2][2] 
			   + A0[1][0]*A0[0][2]*A0[2][1]	   //det(A0)
	     	   + A0[2][0]*A0[0][1]*A0[1][2] 
	     	   - A0[2][0]*A0[1][1]*A0[0][2]
			   - A0[0][0]*A0[1][2]*A0[2][1] 
			   - A0[1][0]*A0[0][1]*A0[2][2]);  //
			
			//fprintf(stderr,"tol = %lf;\tq = %lf;\tp = %lf;\tz = %lf \n",tol,q,p,z);
			
			
			l[2] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*0.)/3.);			//корни характер. многочлена
			l[1] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*1.)/3.);
			l[0] = 2.*cos((1./3.)*acos((3.*q/(2.*p))*sqrt(3./p))-(2.*M_PI*2.)/3.);
			
			//fprintf(stderr,"l0 = %lf l1 = %lf l2 = %lf\n\n",l[0],l[1],l[2]);
																					//
																				  
			for( i=0; i < FL; i++) {             					//упорядочиваем корни характер. многочлена
				for( j = FL-1; j > i; j-- ) {     
					if ( l[j-1] > l[j] ) {
						var = l[j]; l[j] = l[j-1]; l[j-1] = var;
					}
				}
			}														//	
			
			
			a = l[1]-l[0]; // a1=a , b1=b смотри статью
			b = l[2]-l[0];
			var = sqrt(p/3.);
			//fprintf(stderr,"a = %lf b = %lf var = %lf\n\n",a,b,var);
			//fprintf(stderr,"\n");
			
			
			r0 = (-1.)*(2.*sin((var*a*dh)/2.)*sin((var*a*dh)/2.)
				-I*sin(var*a*dh))/a;								//r0 , r1 для exp(..A0)
				 
			r1 = (-1./(a-b))*((2.*sin((var*a*dh)/2.)*sin((var*a*dh)/2.)
			    -I*sin(var*a*dh))/a
			    -(2.*sin((var*b*dh)/2.)*sin((var*b*dh)/2.)
			    -I*sin(var*b*dh))/b);
			//fprintf(stderr,"r1 = %g + i*(%g), r0 = %g + i*(%g)\n\n",creal(r1),cimag(r1),creal(r0),cimag(r0));
			
			for(i=0;i<FL;i++){										//вычисление exp(OMG04)
				for(j=0;j<FL;j++){
					 exA0[i][j] = exp(var*l[0]*I*dh)*(
					 (1.-l[0]*(r0-l[1]*r1))*unit[i][j]
				    +(1./var)*(r0+l[2]*r1)*A0[i][j]
				    +1./(var*var)*r1*sqA0[i][j]);
				}
			}	
			////exA0
			//fprintf(stderr,"exA0 \n");
			//for(i=0;i<FL;i++){
				//for(j=0;j<FL;j++){
					//fprintf(stderr,"%lf + i(%g)\t",
					//creal(exA0[i][j]),
					//cimag(exA0[i][j]));
				//}
				//fprintf(stderr,"\n");
			//}
			//fprintf(stderr,"\n");	
															////
			//for(i=0;i<FL;i++){											//Вычисление [W,[H0,W]]
				//for(j=0;j<FL;j++){
					//sqexA0[i][j] =  
					  //exA0[i][0]*conj(exA0[0][j])
		       	   	 //+exA0[i][1]*conj(exA0[1][j])
					 //+exA0[i][2]*conj(exA0[2][j]);
				//}
			//}		
			
				//fprintf(stderr,"sqexA0 \n");
			//for(i=0;i<FL;i++){
				//for(j=0;j<FL;j++){
					//fprintf(stderr,"%lf + i(%g)\t",
					//creal(sqexA0[i][j]),
					//cimag(sqexA0[i][j]));
				//}
				//fprintf(stderr,"\n");
			//}
			//fprintf(stderr,"\n");
			////end
			for(i=0;i<FL;i++){										//уравнение эволюции нейтрино в точке en
					Psi[i] = (
					    exA0[i][0]*Psi0[0]
				       +exA0[i][1]*Psi0[1]
					   +exA0[i][2]*Psi0[2])*exp(dh*z*I);
			}
																	//
			sqPsi = 0.+I*0.;										//Проверка нормировки решения
			for(i=0;i<FL;i++){
					sqPsi += Psi[i]*conj(Psi[i]);
			}
			//fprintf(stderr,"Psi0 = (%lf+i*(%lf),%lf+i*(%lf),%lf+i*(%lf))\n\n",
			        //creal(Psi0[0]),cimag(Psi0[0]),
			        //creal(Psi0[1]),cimag(Psi0[1]),
			        //creal(Psi0[2]),cimag(Psi0[2]));
			
			//fprintf(stderr,"Psi = (%lf+i*(%lf),%lf+i*(%lf),%lf+i*(%lf))\n\n",
					//creal(Psi[0]),cimag(Psi[0]),
					//creal(Psi[1]),cimag(Psi[1]),
					//creal(Psi[2]),cimag(Psi[2]));
			//fprintf(stderr,"Нормировка 1-|Psi|^2 = %g + i*%g\n\n",1. - creal(sqPsi),cimag(sqPsi));
																	//
				
			for(i=0;i<FL;i++){										//S1 и S2 понадобятся для расчета величины следующего шага
				for(j=0;j<FL;j++){
					S1[i][j] = (-sqrt(3.)/12.)*(vSep-vSem)*H0W[i][j];
				}
			}
			
			for(i=0;i<FL;i++){
				for(j=0;j<FL;j++){
					S2[i][j] = (I*sqrt(3.)/24.)*(vSep-vSem)*(H0H0W[i][j]+(1./2.)*(vSep+vSem)*WH0W[i][j]);
				}
			}														//
			
			for(i=0;i<FL;i++){					    				//sqS1
				for(j=0;j<FL;j++){
					sqS1[i][j] = S1[i][0]*S1[0][j]
								+S1[i][1]*S1[1][j]
								+S1[i][2]*S1[2][j];
				}													//	
			}														
			arr[0]=0.+I*0.; 
			arr[1]=0.+I*0.;
			arr[2]=0.+I*0.;
			for(i=0;i<FL;i++){										
				for(j=0;j<FL;j++){
					arr[i] += ((dh*dh)*S1[i][j] 
							  +dh*dh*dh*S2[i][j]
							  +0.5*dh*dh*dh*dh*sqS1[i][j])*Psi[j];
				}
			}
			//fprintf(stderr,"arr = (%lf +i*(%lf),%lf +i*(%lf),%lf +i*(%lf))\n\n",creal(arr[0]),cimag(arr[0]),creal(arr[1]),cimag(arr[1]),creal(arr[2]),cimag(arr[2]));
			Er = 0.;
			for(i=0;i<FL;i++){
				Er += (double)(arr[i]*conj(arr[i]));
			}														
			
			//fprintf(stderr,"e = %lf\n",e);
			
			if(Er >= tol){											//изменение шага dh
				dh = s*dh*pow((tol/Er),0.3333);														
			}else{
				e = e + dh;
				Psi0[0] = Psi[0];
				Psi0[1] = Psi[1]; 
				Psi0[2] = Psi[2]; 
			}
							
			//fprintf(stderr,"e = %lf  dh = %g Er = %g \nEnd\n\n\n",e,dh,Er);

														// конец цикла по точкам 
		}		
		printf("#e,E,1-|Psi|^2,Psi\n");				
		printf("%lf\t%lf\t %g\t %lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf ",
		  e,E,1.-creal(sqPsi),cimag(sqPsi)
			 ,creal(Psi[0]),cimag(Psi[0])
		     ,creal(Psi[1]),cimag(Psi[1])
		     ,creal(Psi[2]),cimag(Psi[2]));
	
return 0;}

int main (){
	
	double  s12=0.308, s13=0.0234, E=7.05,
			dh=1e-4, tol=1e-4, g = 6.5956e4,
			e0 = 0.1;
	double  H0[FL][FL],W[FL][FL];
	
	ME4(s12,s13,E,dh,tol,g,e0,H0,W);
	
	
return 0;	
}

