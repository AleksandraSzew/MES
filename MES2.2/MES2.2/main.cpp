#include <iostream>
#include <math.h>
#include "Element.h"
#include <iomanip>
using namespace std;




int main()
{double rmax,	
		dr, //krok siatki	
		dtau, //
		c, //efektyenr cieplo wlasciwe
		ro, //gêstoœæ metalu
		k, // wspó³czynnik przewodzenia ciep³a
		alfa, // wspó³czynnik konwekcyjnej wymiany ciep³a
		T0, //- temperatura pocz¹tkowa
 		Talfa; // temperatura otoczenia
		double r=0.0;
		double T1; //temperatura pomocnicza do F
		int iteracja=1;
		double taumax; //koncowy czas
	double pom= 0.0; //zmienna pomocnicza do sumowania
	
	int liczbaWezlow;
	int liczbaElementow;
	// T W O R Z E N I E  W E Z L O W  I  E L E M E N T O W
	cout<<"Wprowadz liczbe wezlow: ";
	cin>>liczbaWezlow;
	liczbaElementow=liczbaWezlow-1;
	cout<<"*********************************************************"<<endl<<endl;
	 
	//U Z U P E L N I E N I E  D A N Y C H
	cout<<"Liczba wezlow: "<<liczbaWezlow<<endl;
	cout<<"Liczba elementow: "<<liczbaElementow<<endl;
	cout<<"rmax: ";
	cin>>rmax;
	dr=rmax/liczbaElementow;

	cout<<"delta tau ";
	cin>>dtau;
	cout<<"c ";
	cin>>c;
	cout<<"ro ";
	cin>>ro;
	cout<<"k ";
	cin>>k;
	cout<<"alfa ";
	cin>>alfa;
	cout<<"t0 ";
	cin>>T0;
	cout<<"talfa ";
	cin>>Talfa;
	cout<<"Tau max: ";
	cin>>taumax;
	

	cout<<endl<<"*********************************************************"<<endl;
	// D A N E 
	cout<<" D  A  N  E  :  "<<endl<<endl;
	cout<<"rmax:               "<<rmax<<" m "<<endl;
	cout<<"Krok siatki:        "<<dr<<" m "<<endl;
	cout<<"delta tau:          "<<dtau<<" s "<<endl;
	cout<<"Wlasnosc materialu: "<<c<<" J/(kg*K) "<<endl;
	cout<<"Gestosc metalu:     "<<ro<<"  kg/m3 "<<endl;
	cout<<"alfa:               "<<alfa<< "  W/m2K "<<endl;
	cout<<"Temp. poczatkowa:   "<<T0<<" oC "<<endl;
	cout<<"Temp. otoczenia:    "<<Talfa<<" oC "<<endl;
	cout<<"*********************************************************"<<endl<<endl;

	Element* tablicaEl = new Element[liczbaElementow];

	//cout<<" D A N E  W E Z L O W  I  E L E M E N T O W : "<<endl;
	for(int i=0; i<liczbaElementow; i++)
	{
		
		tablicaEl[i].wyliczMacierzLokalna(k, dr, dtau, c, ro, r, alfa,rmax, T0,T0, Talfa,liczbaElementow,i);
		r=r+dr;
		
		//cout<<endl<<endl;
	}

	cout<<"*************************************************************************************"<<endl<<endl;

	//MACIERZ GLOBALNA
	cout<<"MACIERZ K: "<<endl;
	double **K= new double*[liczbaWezlow]; 
	double *F= new double[liczbaWezlow];
	for(int i=0;i<liczbaWezlow;i++)
	{
		K[i]= new double[liczbaWezlow];

	}
	for(int i=0;i<liczbaWezlow;i++)
	{
		for(int j=0;j<liczbaWezlow;j++)
	{
		K[i][j]= 0;
		F[i]=0;
		}
	}
	for(int i=0;i<liczbaElementow;i++){
		
		K[i][i]=K[i][i]+tablicaEl[i].K[0][0];
		K[i][i+1]=K[i][i+1]+tablicaEl[i].K[0][1];
		K[i+1][i]=K[i+1][i]+tablicaEl[i].K[1][0];
		K[i+1][i+1]=K[i+1][i+1]+tablicaEl[i].K[1][1];
		F[i] +=tablicaEl[i].F0;
		F[i+1] +=tablicaEl[i].F1;
		}
		
	cout<<"ITERACJA "<<iteracja<<endl;
	iteracja++;
	for(int i=0;i<liczbaWezlow;i++)
	{
		for(int j=0;j<liczbaWezlow;j++)
	{	cout.width(8);
		cout<<setw(8)<<K[i][j]<<", ";
		}
		cout<<"              "<<F[i];
		cout<<endl;
	} 
	cout<<"*************************************************************************************"<<endl<<endl;

	//M A C I E R Z  T E M P E R A T U R

	double * wektorT= new double[liczbaElementow];

	//W Y L I C Z E N I E  t  G A U S S E M	
	for(int i=0; i<liczbaWezlow;i++)
	{
		wektorT[i]=0;
	}

	double * b = new double[liczbaWezlow];
	for(int i =0; i<liczbaWezlow; i++)
	b[i] = -F[i];

	double m = 0;
	double apt;
	
	for (int j = 0; j < liczbaWezlow; j++)
	{
		for (int i = j + 1; i< liczbaWezlow; i++)
		{
			m = K[i][j] / K[j][j];
			K[i][j] -= m*K[j][j]; 
			for (int k = j + 1; k < liczbaWezlow; k++) 
			{
				K[i][k] -= (m*K[j][k]);		
			}
			b[i] -= m*b[j];
		}
	}

	for (int i = 0; i < liczbaWezlow; i++)
	{
		for (int j = i + 1; j < liczbaWezlow; j++)
		{
			m = K[i][j] / K[j][j];
			K[i][j] -= m*K[j][j];
			for (int k = j + 1; k < liczbaWezlow; k++)
			{
				K[i][k] -= m*K[j][k];
			}
			b[i] -= m*b[j];
		}
	}

	for (int i = 0; i < liczbaWezlow; i++)
	{
		for (int j = 0; j < liczbaWezlow; j++)
		{
			if (i == j)
			{
				wektorT[i] = b[i] / K[i][j];
				
			}
		}
	}

	//W Y P I S A N I E  W Y N I K O W  
	double tau=dtau+dtau;
	r=0.0;
	cout<<"Wyniki: "<<endl;
	for(int i=0; i<=liczbaElementow;i++)
	{
		cout<<"t"<<i+1<<": "<<wektorT[i];
		cout<<endl;
	}
	//K O L E J N E  I T E R A C J E 
	while(tau<=taumax){
		//wyzerowanie
		for(int i=0;i<liczbaWezlow;i++)
	{
		for(int j=0;j<liczbaWezlow;j++)
	{
		K[i][j]= 0;
		F[i]=0;
		}
	}
		
		
		cout<<endl<<endl<<"ITERACJA "<<iteracja<<endl<<endl;
		iteracja++;
		for(int i=0; i<liczbaElementow; i++){
			//zamiana temperatur
			
			T0=wektorT[i];
			T1=wektorT[i+1];
		
		
		tablicaEl[i].wyliczMacierzLokalna(k, dr, dtau, c, ro, r, alfa,rmax, T0,T1, Talfa, liczbaElementow,i);
		r=r+dr;
		
		cout<<endl<<endl;
		}
		r=0.0;

		//wyliczenie macierzy globalnej
		for(int i=0;i<liczbaElementow;i++){
		
		K[i][i]=K[i][i]+tablicaEl[i].K[0][0];
		K[i][i+1]=K[i][i+1]+tablicaEl[i].K[0][1];
		K[i+1][i]=K[i+1][i]+tablicaEl[i].K[1][0];
		K[i+1][i+1]=K[i+1][i+1]+tablicaEl[i].K[1][1];
		F[i] +=tablicaEl[i].F0;
		F[i+1] +=tablicaEl[i].F1;
		}

		//wypisanie
		for(int i=0;i<liczbaWezlow;i++){
		for(int j=0;j<liczbaWezlow;j++){	
		cout.width(8);
		cout<<setw(8)<<K[i][j]<<", ";
		}
		cout<<"              "<<F[i];
		cout<<endl;
	} 
	cout<<"*************************************************************************************"<<endl<<endl;
	//W Y L I C Z E N I E  t  G A U S S E M	
	

	double * b = new double[liczbaWezlow];
	for(int i =0; i<liczbaWezlow; i++)
	b[i] = -F[i];

	double m = 0;
	double apt;
	
	for (int j = 0; j < liczbaWezlow; j++)
	{
		for (int i = j + 1; i< liczbaWezlow; i++)
		{
			m = K[i][j] / K[j][j];
			K[i][j] -= m*K[j][j]; 
			for (int k = j + 1; k < liczbaWezlow; k++) 
			{
				K[i][k] -= (m*K[j][k]);		
			}
			b[i] -= m*b[j];
		}
	}

	for (int i = 0; i < liczbaWezlow; i++)
	{
		for (int j = i + 1; j < liczbaWezlow; j++)
		{
			m = K[i][j] / K[j][j];
			K[i][j] -= m*K[j][j];
			for (int k = j + 1; k < liczbaWezlow; k++)
			{
				K[i][k] -= m*K[j][k];
			}
			b[i] -= m*b[j];
		}
	}

	for (int i = 0; i < liczbaWezlow; i++)
	{
		for (int j = 0; j < liczbaWezlow; j++)
		{
			if (i == j)
			{
				wektorT[i] = b[i] / K[i][j];
				
			}
		}
	}

	//W Y P I S A N I E  W Y N I K O W  
	cout<<"Wyniki: "<<endl;
	for(int i=0; i<=liczbaElementow;i++)
	{
		cout<<"t"<<i+1<<": "<<wektorT[i];
		cout<<endl;
	}
	tau +=dtau;
	}
	

system("PAUSE");
}