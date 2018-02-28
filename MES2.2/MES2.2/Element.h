#include <iostream>

using namespace std;

class Element{
	private:
	
	
	double E1; //=-0.5773502692
	double E2; //=0.5773502692
	double rp,//zinterpolowany promien punktu
			wp; //waga punktu
	int np; //zalezne od metody calkowania
public:
	double K[2][2]; //macierz lokalna
	double F0, F1;
	double N[2][2];
	int ID1, ID2;

	Element()
	{

	}
	Element(int ID1,int ID2)
	{
		this->ID1=ID1;
		this->ID2=ID2;
		
	}
	void wprowadzID()
	{	
		cout<<"WPROWADZ WARTOSCI: "<<endl;
		cout<<"ID1: "; cin>>ID1;
		cout<<"ID2: "; cin>>ID2;
		
		
	}
	
	void wyliczMacierzLokalna(double k, double dr, double tau, double c, double ro,double r, double alfa,double rmax, double T0,double T1, double Talfa, int liczbaElementow, int iteracja)
	{
		wp=1;
		np=2;
		//funkcje ksztaltu
		E1=-0.5773502692;
		E2=0.5773502692;
		for(int i=0; i<=2; i++)
		{
			for(int j=0; j<=2;j++)
			N[i][j]=0;
		}



		//zerowanie macierzy 
		for(int i=0;i<=2;i++)
		{
			for(int j=0;j<=2;j++)
			{
				K[i][j]=0;
				F0=0; F1=0;
			}
		}
		//**********************************************************************************
		//OBLICZANIE MACIERZY
		
		
				 for (int ip = 0; ip < np;ip++) {
					 double Alfa=0.0;
					 if(ip==np-1){
						Alfa=alfa;
						}

					 N[0][0]=0.5*(1-E1); //0.78...
					 N[0][1]=0.5*(1-E2); //0.21...
	
					 N[1][0]=0.5*(1+E1); //0.21...
					 N[1][1]=0.5*(1+E2); //0.78...

					rp= N[0][ip]*r + N[1][ip]*(r+dr); 
					
					K[0][0] +=((k/dr)*(rp*wp))+((c*ro*dr)/tau)*N[0][ip]*N[0][ip]*rp*wp;
					if(iteracja==liczbaElementow-1){
					K[1][1] +=((k/dr)*(rp*wp))+((c*ro*dr)/tau)*N[1][ip]*N[1][ip]*rp*wp+2*Alfa*rmax;
					}else{
						K[1][1] +=((k/dr)*(rp*wp))+((c*ro*dr)/tau)*N[1][ip]*N[1][ip]*rp*wp/*+2*Alfa*rmax*/;
					}
					K[0][1] +=-1*((k/dr)*(rp*wp))+((c*ro*dr)/tau)*N[0][ip]*N[1][ip]*rp*wp;
					K[1][0] =K[0][1];
					
					
					F0 +=-1*(c*ro*dr/tau)*((N[0][ip]*T0+N[1][ip]*T1)*N[0][ip]*rp*wp);
					if(iteracja==liczbaElementow-1){
					F1 +=-1*(c*ro*dr/tau)*(N[0][ip]*T0+N[1][ip]*T1)*N[1][ip]*rp*wp -(2*Alfa*rmax* Talfa);
					}else{
						F1 += -1*(c*ro*dr/tau)*(N[0][ip]*T0+N[1][ip]*T1)*N[1][ip]*rp*wp;
					}
			
					
				 }

		//**********************************************************************************
				 //wyswietlenie lokalnych
				 
		//cout<<endl<<"MACIERZ LOKALNA: "<<endl;
		//for(int i=0;i<2;i++)
		//{	
		//	
		//cout<<"|";
		//	for(int j=0;j<2;j++)
		//	{
		//		cout<<K[i][j]<<",  ";
		//	}
		//cout<<"|"<<endl;
		//}
		//cout<<"F"<<endl;
		//cout<<F0<<","<<F1<<endl;
	}
};