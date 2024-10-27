#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <math.h>
typedef std::numeric_limits< double > dbl;

using namespace std;

const double G=6.67*pow(10,-11);        // Constante Gravitationnelle 
const double Mt=5.9*pow(10,24);         // Masse de la terre
const double Rt=6371000;                // Rayon de la terre
const double xt=0, yt=0;                // Position fixe du centre de la terre
const double Tmax=1*3600;               // Temps maximum d'intégration
const double roh0=1.22;                 // Densité volumique de l'atmosphère terrestre à h = 0
const double H= 8100;                   // Longueure de la colonne d'atmosphère
const double Cd=1.7,Cl=0.001;           // Coefficient de trainé & portance
const double sigma=5.971*pow(10,-8);    // Constante de Stephan
const double Tempmax=25000;             // Température maximal à la surface de la météoroïde
const double Q=8*pow(10,6);             // Chaleur spécifique d'ablation de l'astéroïde
const double Ch=0.1;                    // Coefficient de transfert de chaleur de l'astéroïde
const double rohm=3500;                 // Masse volumique de l'asteroïde
const double tau=8.3333*pow(10,-8);     // Coefficient de luminosité
const double S=pow(10,8);               // Rigidité



// Approximation de la terre immobile /
#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279502884197169
#endif

double L,m,phi,mi,ri,I,Mdt;




///////////////////// ATMOSPHERE //////////////////
double atmos(double r){
    return roh0*exp(-(r-Rt)/H);
}




///////////////////// RAYON METEOROIDE - sphère //////////////////////

double ra(double m){
return pow((m*3/(rohm*4*M_PI)),1.0/3);
}




///////////////////// SURFACE EFFICACE SPHERE ///////////////
double surf(double r){return M_PI*pow(r,2);}





////////////////// CALCUL DE TEMPS ///////////////////////////////
int heure(double t){return int (t)/3600;}
int min(double t){return int(t)%3600/60;}
int sec(double t){return int(t)%3600%60%60;}
double milisec(double t){return t-int(t);}




////////////////// POLAIRE VERS CARTESIENNE ////////////////////////
double Xcart(double r,double teta){
return r*cos(teta);
}





double Ycart(double r,double teta){
return r*sin(teta);
}






//////////////////// CARTESIENNE VERS POLAIRE //////////////////////:
double r(double x1, double x0, double y1, double y0){
return sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
}








double teta(double x,double y){
return atan2(y,x);
}






//////////////// FORCE GRAVITATIONNELLE ///////////////////////////
double acceleration(double x1, double x0,double y1,double y0,double m){
// x1 : position fusee
// x0 : position astre
// m : masse astre
return -(G*m)*(x1-x0)/(pow(r(x1,x0,y1,y0),3));
}




////////////// FORCE DE TRAINEE //////////////////////////::
double frot(double x,double y,double vx,double vy,double R, double m){
return Cd*atmos(r(x,xt,y,yt))*surf(R)*pow(r(vx,0,vy,0),2.0)/(2*m);
}






////////////// FORCE DE PORTANCE //////////////////////////::
double port(double x,double y,double vx,double vy,double R,double m){
return Cl*atmos(r(x,xt,y,yt))*surf(R)*pow(r(vx,0,vy,0),1)/(2*m);
}





////////////// FORCE DE TRAINEE FRAGMENTATION //////////////
double FRAG(double x,double y,double vx,double vy,double R, double m){  
return Cd*atmos(r(x,xt,y,yt))*surf(R)*pow(r(vx,0,vy,0),2.0)/(2*m);
}






///////////// FORCE DE PORTANCE FRAGMENTATION //////////////////////////
double portfrag(double x,double y,double vx,double vy,double R, double m){
return Cl*atmos(r(x,xt,y,yt))*surf(R)*pow(r(vx,0,vy,0),1)/(2*m);
}








//FONCTIONS DERIVEES POUR RUNGE KUTTA
void deriv(int n, double t, double z[],double dz[]){
        dz[0]=z[2];         
        dz[1]=z[3];


        if(Cd*atmos(r(z[0],0,z[1],0))*pow(r(z[2],0,z[3],0),2)/2<S){

        dz[2]=acceleration(z[0],xt,z[1],yt,Mt)+cos(atan2(z[3],z[2])+M_PI)*frot(z[0],z[1],z[2],z[3],z[5],z[4])+cos(atan2(z[3],z[2])+M_PI/2)*port(z[0],z[1],z[2],z[3],z[5],z[4]); // Ax
        dz[3]=acceleration(z[1],yt,z[0],xt,Mt)+sin(atan2(z[3],z[2])+M_PI)*frot(z[0],z[1],z[2],z[3],z[5],z[4])+sin(atan2(z[3],z[2])+M_PI/2)*port(z[0],z[1],z[2],z[3],z[5],z[4]); // Ay
        }else{

        dz[2]=acceleration(z[0],xt,z[1],yt,Mt)+cos(atan2(z[3],z[2])+M_PI)*FRAG(z[0],z[1],z[2],z[3],z[5],z[4])+cos(atan2(z[3],z[2])+M_PI/2)*portfrag(z[0],z[1],z[2],z[3],z[5],z[4]); // Ax
        dz[3]=acceleration(z[1],yt,z[0],xt,Mt)+sin(atan2(z[3],z[2])+M_PI)*FRAG(z[0],z[1],z[2],z[3],z[5],z[4])+sin(atan2(z[3],z[2])+M_PI/2)*portfrag(z[0],z[1],z[2],z[3],z[5],z[4]); // Ay

        }


        dz[4]=-min(Ch*atmos(r(z[0],xt,z[1],yt))*surf(ra(z[4]))*pow(r(z[2],0,z[3],0),3)/(2*Q),surf(ra(z[4]))*sigma*pow(Tempmax,4)/Q);
        
       if(Cd*atmos(r(z[0],0,z[1],0))*pow(r(z[2],0,z[3],0),2.0)/2>S){
            dz[6]=Cd*atmos(r(z[0],xt,z[1],yt))*pow(r(z[2],0,z[3],0),2.0)/(2*rohm*z[5]);
            dz[5]=z[6];
        }

}





///////////// INTEGRATEUR SIMPLECTIQUE ?? //////////////////////////




// RUNGE KUTTA///////////////////
void rk4(int n, double x,double y[],double dx,void deriv(int, double, double[], double []))
/*-----------------------------------------
 sous programme de resolution d'equations
 differentielles du premier ordre par
 la methode de Runge-Kutta d'ordre 4
 x = abscisse
 y = valeurs des fonctions
 dx = pas
 n = nombre d'equations differentielles
 deriv = variable contenant le nom du
 sous-programme qui calcule les derivees
 ----------------------------------------*/
{
int i ;
double ddx ;
/* d1, d2, d3, d4 = estimations des derivees
   yp = estimations intermediaires des fonctions */
double d1[n], d2[n], d3[n], d4[n], yp[n];

ddx = dx/2;                /* demi-pas */

deriv(n,x,y,d1) ;          /* 1ere estimation */          

for( i = 0; i< n; i++){ yp[i] = y[i] + d1[i]*ddx ; }
deriv(n,x+ddx,yp,d2) ;     /* 2eme estimat. (1/2 pas) */

for( i = 0; i < n; i++){ yp[i] = y[i] + d2[i]*ddx ; }
deriv(n,x+ddx,yp,d3) ; /* 3eme estimat. (1/2 pas) */

for( i = 0; i< n; i++){ yp[i] = y[i] + d3[i]*dx ;}
deriv(n,x+dx,yp,d4) ;      /* 4eme estimat. (1 pas) */
/* estimation de y pour le pas suivant en utilisant
  une moyenne pond�r�e des d�riv�es en remarquant
  que : 1/6 + 1/3 + 1/3 + 1/6 = 1 */
for( i = 0; i < n ; i++)
 { y[i] = y[i] + dx*( d1[i] + 2*d2[i] + 2*d3[i] + d4[i] )/6 ; }
}















int main(){




int n=7;double h=100000,vi=50000,t=0,dt=0.001,Ec,Ep,Et,Et0,Ec0;
double z[n],dz[n];
L=M_PI/4;               // élévation
phi=M_PI/4;                // angle d'incidence
mi=5.04*pow(10,7);       // masse initiale
ri=ra(mi);              // rayon initial





////////////////////// Initialisation des positions et vitesses ///////////////////////////////
z[0]=(Rt+h)*cos(L) ;                                // position x météoroïde
z[1]=(Rt+h)*sin(L) ;                                // position y météoroïde
//vi=sqrt(G*Mt/sqrt(z[0]*z[0]+z[1]*z[1]));          // vitesse pour trajectoire circulaire
z[2]= vi*cos(M_PI/2-L+phi);                         // vitesse x météoroïde
z[3]= -vi*sin(M_PI/2-L+phi);                        // vitesse y météoroïde
z[4]=mi;                                            // masse météoroïde
z[5]=ra(z[4]);                                      // Rayon initial
z[6]=0;                                             // Vitesse de variation de rayon initiale
cout<<"Position initiales :\nx = "<<z[0]/(Rt+h)<<"\ny = "<<z[1]/(Rt+h)<<endl;
cout<<"Vitesses initiales :\nVx = "<<z[2]<<"\nVy = "<<z[3]<<endl;






///////////////// APPEL FICHIER TRAJECTOIRE & ENERGIE /////////////////////////////
ofstream traj("astev50.traj");
ofstream ener("astev50.ener");
ifstream in; 





/////////////////////// CALCUL ENERGIE INITIALE ////////////////////////////
        Ec=(z[2]*z[2]+z[3]*z[3])*z[4]/2;
        Ec0=Ec;
        Ep=-G*z[4]*Mt/(sqrt(z[0]*z[0]+z[1]*z[1]));
        Et0=Ec+Ep;
double Eci,hi,deltaEc;
int bol=0;




///////// Boucle principale //////////////



while((z[0]*z[0]+z[1]*z[1])>=Rt*Rt && t<=Tmax){
            ////////// Calcul d'énergie avant RK4 //////////
          Eci=Ec;
          hi=r(z[0],xt,z[1],yt)-Rt;

        traj<<::fixed << std::setprecision(14)<<t<<" "<<(r(z[0],xt,z[1],yt)-Rt)/1000<<" "<<r(z[2],0,z[3],0)<<" "<<z[4]<<" "<<z[5]<<" "<<24.3-2.5*log10(I)<<endl;
        rk4(n,t,z,dt,deriv);

          Ec=(z[2]*z[2]+z[3]*z[3])*z[4]/2;
          Ep=-G*Mt*z[4]/(sqrt(z[0]*z[0]+z[1]*z[1]));
          Et=Ec+Ep;
          deltaEc=(Ec-Eci)/(r(z[0],xt,z[1],yt)-Rt-hi);


       ener<<::fixed << std::setprecision(14)<<t<<" "<<(r(z[0],xt,z[1],yt)-Rt)/1000<<" "<<abs(Ec-Ec0)/(4.185*pow(10,15))<<" "<<deltaEc/(4.185*pow(10,15))<<" "<<Et/(4.185*pow(10,15))<<endl;
        
        L=atan2(z[1],z[0]);
        phi=-atan2(z[3],z[2])-M_PI/2+L;
        t=t+dt;
        if(r(z[0],0,z[1],0)>30000){
            I=tau*pow(r(z[2],0,z[3],0),3)*Ch*atmos(r(z[0],xt,z[1],yt))*surf(ra(z[4]))*pow(r(z[2],0,z[3],0),3)/(2*Q);
        }
        else {       
            I=tau*pow(r(z[2],0,z[3],0),3)*surf(ra(z[4]))*sigma*pow(Tempmax,4)/Q;
        }
        if(Cd*atmos(r(z[0],0,z[1],0))*pow(r(z[2],0,z[3],0),2.0)/2<S && bol==0){z[5]=ra(z[4]);}
        if(Cd*atmos(r(z[0],0,z[1],0))*pow(r(z[2],0,z[3],0),2.0)/2>S){bol=bol+1;}
}












printf("\n\n\n\n\n//////////////////////////////////////////////////////////////\n\n\n");
if((z[0]*z[0]+z[1]*z[1])<=Rt*Rt){
    printf("Impact sur la Terre !\n");
        cout<<heure(t)<<" heure(s) "<<min(t)<<" minute(s) "<<sec(t)+milisec(t)<<" secondes."<<endl;
    }
else if(t>=Tmax){
    printf("Temps maximum atteint.\n");
        cout<<heure(t)<<" heure(s) "<<min(t)<<" minute(s) "<<sec(t)+milisec(t)<<" secondes."<<endl;
}
printf("\n\n//////////////////////////////////////////////////////////////\n\n\n\n\n\n\n");

return 0;
}

