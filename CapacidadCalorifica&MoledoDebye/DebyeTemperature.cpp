#include <iostream>
#include <cmath>
 const double Ti = 77.0;
 const double TfAl = 291.0;
 const double TfPb = 292.0;
 const double TfZn = 292.0;
 const double TfCu = 290.0;
 const double TfSn = 295.0;
 const double TfBronze = 291.0;
const double MasaMolarAl = 26.981;
const double MasaMolarPb = 207.2;
const double MasaMolarCu = 63.546;
const double MasaMolarSn =  118.71;
const double MasaMolarZn =  65.38 ;
const double MasaMolarBronze = 0.88*MasaMolarCu + 0.12*MasaMolarSn;
const double K_bxAvo=8.314469723;
const double MasaAl=66.5;
const double MasaPb=229.1;
const double MasaCu=63.3;
const double MasaSn=62.8;
const double MasaZn=64.4;
const double MasaBronze=26.3;
const double QAl=10447.5;
const double QPb=6109.3;
const double QCu=4417.8;
const double QSn=2825.8;
const double QBronze=1950.2;
const double QZn=4776;
const double KAl=9*K_bxAvo*MasaAl/MasaMolarAl;
const double KPb=9*K_bxAvo*MasaPb/MasaMolarPb;
const double KCu=9*K_bxAvo*MasaCu/MasaMolarCu;
const double KSn=9*K_bxAvo*MasaSn/MasaMolarSn;
const double KZn=9*K_bxAvo*MasaZn/MasaMolarZn;
const double KBronze=9*K_bxAvo*MasaBronze/MasaMolarBronze;
const int N1=10000;
using namespace std;

double f1(double a)
{ return pow(a,3)/(exp(a)-1);
}

double f2(double a)
{ return pow(a,4)*exp(a)/((exp(a)-1)*(exp(a)-1));
}

double IntegralPorSimpsonf1( double b, double a)
{
  
  // int N = rint( (b-a)*N1/100);
       int N=N1;
  //double h = 0.00001;
  /*
  int i; double xi,suma,integral;

suma=0;
 xi=a;
 for(i=0;xi<=b+h;i++)
    {
    xi=a+i*h;
    if((i==0)||(xi>b))
    suma+=f1(xi);
    else if(i%2==0)
    suma+=2*f1(xi);
    else
suma+=4*f1(xi);
}
integral=h/3*suma;
return integral;
  */

 double h=(b-a)/N;
int i; double xi,suma,integral;

suma=0;
for(i=0;i<=N;i++)
    {
    xi=a+i*h;
    if((i==0)||(i==N))
    suma+=f1(xi);
    else if(i%2==0)
    suma+=2*f1(xi);
    else
suma+=4*f1(xi);
}
integral=h/3*suma;
return integral;

}

double IntegralPorSimpsonf2( double b, double a)
{
  
  // int N = rint( (b-a)*N1/100);
       int N=N1;
  //double h = 0.00001;
  /*
  int i; double xi,suma,integral;

suma=0;
 xi=a;
 for(i=0;xi<=b+h;i++)
    {
    xi=a+i*h;
    if((i==0)||(xi>b))
    suma+=f1(xi);
    else if(i%2==0)
    suma+=2*f1(xi);
    else
suma+=4*f1(xi);
}
integral=h/3*suma;
return integral;
  */

 double h=(b-a)/N;
int i; double xi,suma,integral;

suma=0;
for(i=0;i<=N;i++)
    {
    xi=a+i*h;
    if((i==0)||(i==N))
    suma+=f2(xi);
    else if(i%2==0)
    suma+=2*f2(xi);
    else
suma+=4*f2(xi);
}
integral=h/3*suma;
return integral;

}


int main(void)
{
  double a=0.000000000001, dt=0.01, FTf,FTi,Td,holdnew,holdold;
 
   holdold=1;
  
for(Td=1;Td<450;Td+=dt){
 
   holdnew = (KAl/pow(Td,3))*(pow(TfAl,4)*IntegralPorSimpsonf1(Td/TfAl,a) - pow(Ti,4)*IntegralPorSimpsonf1(Td/Ti,a))-QAl;
 

   if(holdold*holdnew<0){
   cout<<Td<<endl;
   FTf =  (KAl/pow(Td,3))*pow(TfAl,3)*IntegralPorSimpsonf2(Td/TfAl,a);
      FTi =  (KAl/pow(Td,3))*pow(Ti,3)*IntegralPorSimpsonf2(Td/Ti,a);
      cout<<FTf/MasaAl<<" "<<FTi/MasaAl<<endl;
      dt=1;
   }
 holdold=holdnew;


}

return 0;

}
