// inicia los parametros que vamos a estimar, y calcula el primer chi^2
void Inicie(void){
  cout<<"Parametros iniciales"<<endl;
if(esferas==0){
  if((id_stack==0)||(id_stack==1)){
    param[0]=2.8;		sigma[0]=0.01;		//	dc
    param[1]=1.25;		sigma[1]=0.001;	//	rs
    param[2]=1.05;		sigma[2]=0.01;	//	rv
    param[3]=2.1;		sigma[3]=0.05;	//	alpha
    param[4]=3.6;		sigma[4]=0.04;	//	beta
  }
  if(id_stack==0)
    param[0]=2.8;		sigma[0]=0.15;		//	dc
  if((id_stack==2)||(id_stack==3)){
    param[0]=0.94;		sigma[0]=0.05;		//	dc
    param[1]=1.07;		sigma[1]=0.0015;	//	rs
    param[2]=1.1;		sigma[2]=0.015;	//	rv
    param[3]=4.9;		sigma[3]=0.04;	//	alpha
    param[4]=6.6;		sigma[4]=0.04;	//	beta
  }
  if((id_stack==4)||(id_stack==5)){
    param[0]=0.854;		sigma[0]=0.0004;		//	dc
    param[1]=1.12;		sigma[1]=0.002;	//	rs
    param[2]=1.025;		sigma[2]=0.002;	//	rv
    param[3]=5.6;		sigma[3]=0.03;	//	alpha
    param[4]=7.8;		sigma[4]=0.02;	//	beta
  }
  if((id_stack==6)||(id_stack==7)){
    param[0]=0.87;		sigma[0]=0.004;		//	dc
    param[1]=2.26;		sigma[1]=0.1;	//	rs
    param[2]=0.94;		sigma[2]=0.002;	//	rv
    param[3]=1.6;		sigma[3]=0.1;	//	alpha
    param[4]=8.4;		sigma[4]=0.1;	//	beta
  }
}
else{
double fac=0.5;
  if((id_stack==0)){
    param[0]=0.89;		sigma[0]=0.04*fac;	//	dc
    param[1]=0.25;		sigma[1]=0.03*fac;	//	A
    param[2]=1.4;		sigma[2]=0.06*fac;	//	rv
    param[3]=3.9;		sigma[3]=0.6*fac;	//	alpha
    param[4]=4.8;		sigma[4]=0.5*fac;	//	beta
  }
  if((id_stack==1)){
    param[0]=0.865;		sigma[0]=0.01*fac;	//	dc
    param[1]=0.23;		sigma[1]=0.02*fac;	//	A
    param[2]=1.34;		sigma[2]=0.04*fac;	//	rv
    param[3]=4.8;		sigma[3]=0.6*fac;	//	alpha
    param[4]=5.8;		sigma[4]=0.6*fac;	//	beta
  }
  if((id_stack==2)){
    param[0]=0.87;		sigma[0]=0.03*fac;	//	dc
    param[1]=0.2;		sigma[1]=0.02*fac;	//	A
    param[2]=1.25;		sigma[2]=0.04*fac;	//	rv
    param[3]=5.6;		sigma[3]=0.8*fac;	//	alpha
    param[4]=6.8;		sigma[4]=0.7*fac;	//	beta
  }
  if((id_stack==3)){
    param[0]=0.87;		sigma[0]=0.02*fac;	//	dc
    param[1]=0.16;		sigma[1]=0.02*fac;	//	A
    param[2]=1.19;		sigma[2]=0.03*fac;	//	rv
    param[3]=6.5;		sigma[3]=0.9*fac;	//	alpha
    param[4]=7.8;		sigma[4]=0.9*fac;	//	beta
  }
  if((id_stack==4)){
    param[0]=0.86;		sigma[0]=0.03*fac;	//	dc
    param[1]=0.12;		sigma[1]=0.02*fac;	//	A
    param[2]=1.15;		sigma[2]=0.03*fac;	//	rv
    param[3]=8.0;		sigma[3]=1.4*fac;	//	alpha
    param[4]=9.2;		sigma[4]=1.4*fac;	//	beta
  }
  if((id_stack==5)){
    param[0]=0.855;		sigma[0]=0.03*fac;	//	dc
    param[1]=0.04;		sigma[1]=0.05*fac;	//	A
    param[2]=1.15;		sigma[2]=0.03*fac;	//	rv
    param[3]=9.5;		sigma[3]=1.8*fac;	//	alpha
    param[4]=10.5;		sigma[4]=1.8*fac;	//	beta
  }
  if((id_stack==6)){
    param[0]=0.85;		sigma[0]=0.02*fac;	//	dc
    param[1]=-0.1;		sigma[1]=0.5*fac;	//	A
    param[2]=1.06;		sigma[2]=0.04*fac;	//	rv
    param[3]=8.;		sigma[3]=2.5*fac;	//	alpha
    param[4]=11.;		sigma[4]=2.*fac;	//	beta
  }
  if((id_stack==7)){
    param[0]=0.845;		sigma[0]=0.02*fac;	//	dc
    param[1]=-0.4;		sigma[1]=0.5*fac;	//	A
    param[2]=1.02;		sigma[2]=0.04*fac;	//	rv
    param[3]=13.;		sigma[3]=3.8*fac;	//	alpha
    param[4]=15.;		sigma[4]=3.5*fac;	//	beta
  }
  if((id_stack==8)){
    param[0]=0.835;		sigma[0]=0.02*fac;	//	dc
    param[1]=-35.12;		sigma[1]=5.5*fac;	//	A
    param[2]=0.91;		sigma[2]=0.04*fac;	//	rv
    param[3]=40.;		sigma[3]=9.8*fac;	//	alpha
    param[4]=40.0;		sigma[4]=9.6*fac;	//	beta
  }
}

#if deltac_constante>0
    param[0]=0.8673581;		sigma[0]=0.0;	//	dc
#endif

  Chi();				//	compara con observaciones
  ChiOld=ChiNew;
  ChiMax=ChiNew;
  ChiMin=ChiNew;
}
