#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>           //M_PI ne marche pas sans math.h et define, dans cmath          
                            //je n'arrive pas a utiliser M_PI
using namespace std;        //donc j'ai utiliser cette alternative     

vector<vector<char>> graph;

void print_error(string message) {
    cout << message;
    cout << endl;
    exit(0);
}

void ask_signal(string message, string& signal_t) {  //signal input

    cin >> signal_t;
    
    if (signal_t !=("SAWTOOTH") && signal_t !=("SQUARE") && signal_t !=("TRIANGLE")) {
        print_error(message);
    }

}

void ask_nbn(int& nbN, string message) {   //nb de termes input

    cin >> nbN;
    
    if (nbN <= 0) {
        print_error(message);
    }
}

void ask_time_interval(double& tmin, double& tmax, string message, string message2) {
                                        //intervalle de temps input
    cin >> tmin;
    cin >> tmax;
    
    if (tmin > tmax) {
        print_error(message);
    }
   
    if ((tmin < 0 || tmin > 1) || (tmax < 0 || tmax > 1)) {
        print_error(message2);
    }
}

void ask_amplitude(double& signalmin, double& signalmax, string message) {
                                                //amplitude input
    cin >> signalmin;
    cin >> signalmax;

    if (signalmin > signalmax) {
        print_error(message);
    }
}

void ask_nbl(int& nbL, int& nbC, string message, string message2) {
                                //nbL input

    cin >> nbL;
    if (nbL <= 2) {
        print_error(message);
    }

    if (nbL % 2 == 0) {
        print_error(message2);
    }
    nbC = 2 * nbL - 1;

}

void approxsawtooth (int nbN, double t, double& s_t){       //calcul approx sawtooth
    double sum(0);
    
    for (int k(1); k <= nbN ; ++k ){
        
        sum += ((pow(-1,k)/k)*sin(2*(M_PI*k)*(t-0.5)));

    }
    
    s_t= (-2/M_PI)*sum;
}

void approxsquare (int nbN,double t,double& s_t){           //calcul approx square
    double sum(0);

    for (int k(1);k<=nbN;++k ){
        
        sum += ((sin(2*M_PI*(2*k-1)*t))/(2*k-1));
    }
    
    s_t= (4/M_PI)*sum;
}

void approxtriangle (int nbN,double t,double& s_t){         //calcul approx triangle
    double sum(0);

    for (int k(1);k<=nbN;++k ){
        
        sum+= (((pow(-1,k))/(pow(2*k-1,2)))*sin(2*M_PI*(2*k-1)*(t-0.25)));
    }
    
    s_t= (-8/pow(M_PI,2))*sum;
}

void thsawtooth (double t,double& s_t){             //calcul theorique sawtooth
    double valth(0);
    valth=2*(t-0.5);
    
    s_t=valth;
}

void thsquare(double t,double& s_t){                //calcul theorique square
    double valth(0);
    if(t>0 && t< 0.5){
        valth=1;
    }
    if(t>0.5 && t<1){
        valth=-1;
    }
    if(t==0 || t==0.5 || t==1){
        valth=0;
    }

    s_t=valth;
}

void thtriangle(double t, double& s_t,double nbL,double nbC){  //calcul theorique tri
    double valth(0);
    if (t<0.5){
        valth=4*(t-0.25);
    }
    if (t>0.5){
        valth=-4*(t-0.75);
    }
    if(t==0.5){
        valth=1;
    }
    
    s_t=valth;
}


void calculthsaw(string signal,double tmin,double tmax,double max,
                  double min,double nbC,double nbL) {
    int j(0);                               //calcul graphe theorique sawtooth
    double delta_t; 
    double delta_s;
    delta_t=(tmax-tmin)/(nbC-1);
    delta_s=(max-min)/(nbL-1);
    
    graph = vector<vector<char>>(nbL, vector<char>(nbC, ' '));

    for (int j(1) ; j<=nbC-2 ; ++j) {
        double result;
        double t_th;
        t_th=tmin+j*delta_t;
        
        thsawtooth(t_th,result);
        
        double v = ((result - min) / delta_s) + 0.5;
        int z = nbL - 1 - (int)v;     //sans le (int)v (si je met seulement v) 
        
        if (z>=0 && z<nbL){    //mes graphiques seront decale de 2 unite vers le haut!
        graph[z][j]='+';
        }  
    }
}

void calculthsqr(string signal,double tmin,double tmax,double max,
                  double min,double nbC,double nbL) {
    int j(0);                               //calcul graphique theorique square
    double delta_t; 
    double delta_s;
    delta_t=(tmax-tmin)/(nbC-1);
    delta_s=(max-min)/(nbL-1);
    
    graph = vector<vector<char>>(nbL, vector<char>(nbC, ' '));

    for (int j(0) ; j<=nbC-1 ; ++j) {
        double result;
        double t_th;
        t_th=tmin+j*delta_t;
        
        thsquare(t_th,result);
        
        double v = ((result - min) / delta_s) + 0.5;
        int z = nbL - 1 - (int)v;      
        
        if (z>=0 && z<nbL){           
        graph[z][j]='+';
        }  
    }   
}

void calculthtri(string signal,double tmin,double tmax,double max,
                  double min,double nbC,double nbL) {
                                    //calcul graphique theorique triangle
    int j(0);
    double delta_t; 
    double delta_s;
    delta_t=(tmax-tmin)/(nbC-1);
    delta_s=(max-min)/(nbL-1);
    graph = vector<vector<char>>(nbL, vector<char>(nbC, ' '));

    for (int j(0) ; j<=nbC-1 ; ++j) {
        
       double result;
       double t_th;
        
       t_th=tmin+j*delta_t;
        
       double nbL1 =nbL;
       double nbC1 =nbC;
        
       thtriangle(t_th,result,nbL1,nbC1);
        
       double v = ((result - min) / delta_s) + 0.5;
        
       int z = nbL - 1 - int(v);                           
       
       if (z>=0 && z<=nbL){         
        graph[z][j]='+';
       }       
    }
}
      
void calculapprox(string signal,double tmin,double tmax,double max,
                  double min,double nbC,double nbL,double nbN) {
                                    //calcul graphe approx de tous les signaux
    int j(0);
    double delta_t,delta_s;
    delta_t=(tmax-tmin)/(nbC-1);
    delta_s=(max-min)/(nbL-1);
    for (int j(0) ; j<=nbC-1 ; ++j) {
        double t;
        double s_t;
        t=tmin+j*delta_t;
        if (signal=="SAWTOOTH"){
            approxsawtooth(nbN,t,s_t);         
        }
        else if (signal=="SQUARE"){
            approxsquare(nbN,t,s_t);
        }
        else if (signal=="TRIANGLE"){
            approxtriangle(nbN,t,s_t);
        } 
        double v = ((s_t - min) / delta_s) + 0.5;
        int i = nbL - 1 - (int)v; 
        if (i >= nbL) {      
            i = nbL - 1;
        }
        else if (i >= 0 && i < nbL) {
            graph[i][j] = '*';
        }
    }
}


void afficheGraph(double max,double min, int nbL) {
    cout << string(graph[0].size(), '-') << endl;       //fonction affichage graphique

    double delta_s = (max - min) / (nbL - 1);
    int y_zero_row = nbL - 1 - int((-min / delta_s) + 0.5);

    for (int i(0); i < graph.size(); ++i) {
        for (int j(0); j < graph[0].size(); ++j) {
            if (i == y_zero_row && graph[i][j] == ' ') {
                cout << '.'; 
            } 
            else {
                cout << graph[i][j];
            }
        }
        cout << endl;
    }
    cout << string(graph[0].size(), '-'); 
}

double MaxDicho(int nbN, double debut, double fin, const string& signal) {
    double gauche = debut;                      //recherche du max par dichotomie
    double droite = fin;
    const double EPSIL_DICHO(1e-9);
    double milieu, s_t;

    while (droite - gauche > EPSIL_DICHO) {
        milieu = (gauche + droite) / 2;
        double s_t_milieu, s_t_droite;

        if (signal == "SAWTOOTH") {
            approxsawtooth(nbN, milieu, s_t_milieu);
            approxsawtooth(nbN, milieu + EPSIL_DICHO, s_t_droite);
        } 
        else if (signal == "SQUARE") {
            approxsquare(nbN, milieu, s_t_milieu);
            approxsquare(nbN, milieu + EPSIL_DICHO, s_t_droite);
        } 
        else if (signal == "TRIANGLE") {
            approxtriangle(nbN, milieu, s_t_milieu);
            approxtriangle(nbN, milieu + EPSIL_DICHO, s_t_droite);
        }

        if (s_t_milieu < s_t_droite)
            gauche = milieu; 
        else
            droite = milieu;
    }

    if (signal == "SAWTOOTH")
        approxsawtooth(nbN, (gauche + droite) / 2, s_t);
    else if (signal == "SQUARE")
        approxsquare(nbN, (gauche + droite) / 2, s_t);
    else if (signal == "TRIANGLE")
        approxtriangle(nbN, (gauche + droite) / 2, s_t);

    return s_t;
}

int main() {
 const string BAD_SIGNAL("Error: signal type must be SAWTOOTH, SQUARE or TRIANGLE");
 const string NBN_TOO_SMALL("Error: the number of sine terms must be greater than 0");
 const string NBL_TOO_SMALL("Error: the number of lines must be greater than 2");
 const string TIME_MIN_MAX("Error: time max is not bigger than min");
 const string SIGNAL_MIN_MAX("Error: signal max is not bigger than min");
 const string WRONG_TIME_VAL("Error: both time values must belong to [0., 1.]");
 const string NBL_NOT_ODD("Error: the number of lines must be odd");  
 
 string signal;
 int nbN,nbL,nbC;
 double tmin, tmax,min,max,SMaxVal;
 
 ask_signal(BAD_SIGNAL, signal);
 ask_nbn(nbN, NBN_TOO_SMALL);
 ask_time_interval(tmin, tmax, TIME_MIN_MAX, WRONG_TIME_VAL);
 ask_amplitude(min, max, SIGNAL_MIN_MAX);
 ask_nbl(nbL, nbC, NBL_TOO_SMALL, NBL_NOT_ODD);
 
 if(signal=="SAWTOOTH") calculthsaw(signal,tmin,tmax,max,min,nbC,nbL);
 if(signal=="SQUARE") calculthsqr(signal,tmin,tmax,max,min,nbC,nbL);
 if(signal=="TRIANGLE") calculthtri(signal,tmin,tmax,max,min,nbC,nbL);
 
 calculapprox (signal,tmin,tmax,max,min,nbC,nbL,nbN);
 afficheGraph(max,min,nbL);
 
 if (signal == "SAWTOOTH") 
    SMaxVal=MaxDicho(nbN, 1 - 1.0 / (2 * nbN + 1), 1, signal); 

 else if (signal == "SQUARE") 
    SMaxVal = MaxDicho(nbN, 0, 1.0 / (2 * nbN + 1), signal); 
 
 else if (signal == "TRIANGLE") 
    SMaxVal= MaxDicho (nbN,0.5-1.0/(2*(2*nbN+1)),0.5+1.0/(2*(2*nbN+1)),signal);
 
 cout << setprecision(8) << fixed << endl << SMaxVal << endl;
 
 return 0;
}