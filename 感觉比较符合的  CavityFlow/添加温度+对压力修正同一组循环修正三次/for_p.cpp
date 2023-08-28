#include<iostream>
# include<cmath>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>
using namespace std;

double p[96][96],u[95][96],ru[95][96],ru2[95][96],v[96][95],rv[96][95],rv2[96][95];
double p_[96][96],p1[96][96],pg[96][96]; 
double ru_[95][96],rv_[96][95];
double T[95][95],T1[95][95],uT[95][95],vT[95][95];
int i,j;
double u_,un,v_,vn;

double rho = 1;
double niu = 1; 
double A,B;
float a,b,c,d;
double dx,dy,dt;
double C,D,E;//二阶偏微分温度 

double Lx,Ly;



int main()
{
	
	dx = 1;

	dy = 1;

	dt = 0.001;
	rho = 1;
	niu = 10;
	
	Lx = dx*94;
	Ly = dy*94;
	
	
	//温度初始化	
for(i=0;i<=94;i++)
{
	for(j=0;j<=94;j++)
	{
		
	T[i][j] = 0;	
	uT[i][j] = u[i][j]*T[i][j];
	vT[i][j] = v[i][j]*T[i][j];
	
	T[i][94] = 1;
	uT[i][94] = u[i][94]*T[i][94];
	vT[i][94] = v[i][94]*T[i][94];
		
	}
	
		}	
	
	
	
	for(i=0;i<=94;i++)
	{
		for(j=0;j<=95;j++)
		
	{
		u[i][j] = 0;
		ru[i][j] = u[i][j]*rho;
		ru2[i][j] = ru[i][j]*u[i][j];
		u[i][95] = 1;
		ru[i][95] = u[i][95]*rho;
		ru2[i][95] = ru[i][95]*u[i][95];
		
	}
		
	}
	
	for(i=0;i<=95;i++)
	{
		for(j=0;j<=94;j++)
		
	{
		v[i][j] = 0;
		rv[i][j] = v[i][j]*rho;
		rv2[i][j] = rv[i][j]*v[i][j];
		
	}
		
	}
	
	for(i=0;i<=95;i++)
	{
		for(j=0;j<=95;j++)
		{
			p[i][j] = 0;
			
		}
		
		
	}
	
	for(i=0;i<=95;i++)
	{
		for(j=0;j<=95;j++)
		{
			p_[i][j] = 0;
			
		}
		
		
	}
	
	
	
	
	
	for(int t=0;t<1300000;t++)
	{
	
	cout << t<<endl;
	
		for(i=1;i<94;i++)
	{
		for(j=1;j<95;j++)
		{
		
		v_ = 1/2*(v[i][j]+v[i+1][j]);
		 vn = 1/2*(v[i][j-1]+v[i+1][j-1]);
		A = -1*(((ru2[i+1][j]-ru2[i-1][j])/(2*dx))+((ru[i][j+1]*v_-ru[i][j-1]*vn)/(2*dy)))+niu*((u[i+1][j]-2*u[i][j]+u[i-1][j])/(dx*dx)+(u[i][j+1]-2*u[i][j]+u[i][j-1])/(dy*dy));	
		ru_[i][j] = ru[i][j] + A*dt-(dt/dx)*(p[i+1][j]-p[i][j]);	//ru在n+1时刻的数值 
			
		}
				
	}
	
	//更新变量 Ux 
	
	for(i=1;i<94;i++)
	{
		for(j=1;j<95;j++)
		{
			ru[i][j] = ru_[i][j];
			u[i][j] = ru[i][j]/rho;
			ru2[i][j] = ru[i][j]*u[i][j];
			// v值还没有更新 
			
		}
		
		
	 } 
	
	//左右边界 
	for(j=1;j<95;j++)
	{
		ru_[0][j] = 0; 
		//ru_[0][j] = 2*ru_[1][j]-ru_[2][j]; //边界条件外推 
		ru[0][j] = ru_[0][j];
		u[0][j] = ru[0][j]/rho; 
		ru2[0][j] = ru[0][j]*u[0][j];	
	
		ru_[94][j] = 0;
		//ru_[94][j] = 2*ru_[93][j]-ru_[92][j]; //边界条件外推
		ru[94][j] = ru_[94][j];
		u[94][j] = ru[94][j]/rho;
		ru2[94][j] = ru[94][j]*u[94][j]; 
	
	
	
	
	
	 } 
	
	
	//上下边界
	for(i=0;i<=94;i++) 
	{
		ru_[i][0] = 0;
	//	ru_[i][0] = 2*ru_[i][1]-ru_[i][2]; 
		ru[i][0] = ru_[i][0];
		u[i][0] = ru[i][0]/rho; 
		ru2[i][0] = ru[i][0]*u[i][0];
		
		
		
		ru_[i][95] = 1;		
		ru[i][95] =ru_[i][95];
		u[i][95] = ru[i][95]/rho; 
		ru2[i][95] = ru[i][95]*u[i][95];
	}
	
	//========================================================================================//
	//计算 Uy 
	
	for(i=1;i<95;i++)
	{
		for(j=1;j<94;j++)
		{
		u_ = 1/2*(u[i][j+1]+u[i][j]); 
	    un = 1/2*(u[i-1][j+1]+u[i-1][j]);
	    
	    B = -1*((rv[i+1][j]*u_-rv[i-1][j]*un)/(2*dx)+(rv2[i][j+1]-rv2[i][j-1])/(2*dy))+niu*((v[i+1][j]-2*v[i][j]+v[i-1][j])/(dx*dx)+(v[i][j+1]-2*v[i][j]+v[i][j-1])/(dy*dy));
		
		rv_[i][j] = rv[i][j]+B*dt-(dt/dy)*(p[i][j+1]-p[i][j]);
			
		}
		
		
	} 
	
	//更新变量 Uy
	
	 for(i=1;i<95;i++)
	 {
	 	for(j=1;j<94;j++)
	 	{
	 		rv[i][j] = rv_[i][j];
			v[i][j] = rv[i][j]/rho;
			rv2[i][j] = rv[i][j]*v[i][j];
	 		
	 		
		 }
	 	
	 	
	  } 
	
	//左右边界 
	
	for(j=1;j<94;j++)
	{
		rv_[0][j] = 0;
		//rv_[0][j] = 2*rv_[1][j]-rv_[2][j];
		rv[0][j] = rv_[0][j];
		v[0][j] = rv[0][j]/rho;
		rv2[0][j] = rv[0][j]*v[0][j];
		
		
		rv_[95][j] = 0;
	//	rv_[95][j] = 2*rv_[94][j]-rv_[93][j];
		rv[95][j] =rv_[95][j];
		v[95][j] = rv[95][j]/rho; 
		rv2[95][j] = rv[95][j]*v[95][j];
	 } 
	
	//上下边界
	for(i=0;i<=95;i++)
	{
	
	rv_[i][0] = 0;
	//rv_[i][0] = 2*rv_[i][1]-rv_[i][2];
	rv[i][0] = rv_[i][0];
	v[i][0] = rv[i][0]/rho;
	rv2[i][0] = rv[i][0]*v[i][0];
	
	rv_[i][94] = 0; 
	//rv_[i][94] = 2*rv_[i][93]-rv_[i][92];
	rv[i][94] = rv_[i][94];
	v[i][94] = rv[i][94]/rho;
	rv2[i][94] = rv[i][94]*v[i][94];
	
	
	
   }
//=============================================================================//
//压力修正(内部流场点) 
for(int l=0;l<3;l++)  //压力修正三次 
{


for(i=0;i<=95;i++)
{
	
	for(j=0;j<=95;j++)
	{
		
		pg[i][j] = p[i][j];
		
	 } 
	
 } 









for(i=1;i<95;i++)
{
for(j=1;j<95;j++)
{

a = 2*((dt)/(dx*dx)+(dt)/(dy*dy));
b = -1*(dt/(dx*dx));
c = -1*(dt/(dy*dy));
d = (ru[i][j]-ru[i-1][j])+(rv[i][j]-rv[i][j-1]);


if(isnan(d))
    {
        
        d = 0;       //在无穷小时会出现 nan; 
    }
 else
 {
 	d = d;
 	
   }  
    
/*    
if(isnan(d))
    {
        cout << d << endl;
       
    }
    */
/*
else
{
	cout << "yes" << endl;
	
}

*/


p[i][j] = 1e-4*(-1/a)*(b*p[i+1][j]+b*p[i-1][j]+c*p[i][j+1]+c*p[i][j-1]+d); //压力 p 的误差会扩散，这里就乘上1e-4修正，自己试出来的;

}
}

for(i=1;i<95;i++)
{
for(j=1;j<95;j++)
{

p1[i][j] = pg[i][j]+0.035*p[i][j]; //再乘上0.035的松弛因子系数; 
}
}

//边界上的压力

//左右边界

for(j=1;j<95;j++)
{
	p1[0][j] = p1[1][j];
	p1[95][j] = p1[94][j];
 } 

//上下边界

for(i=0;i<=95;i++)
{
	p1[i][0] = p1[i][1];
	p1[i][95] = p1[i][94];
	
 } 


//更新 p 值

for(i=0;i<=95;i++)
{
	
	for(j=0;j<=95;j++)
	{
		
		p[i][j] = p1[i][j];
		
	 } 
	
 } 

//利用更新的 p值重新计算 u、v; 

	for(i=1;i<94;i++)
	{
		for(j=1;j<95;j++)
		{
		
		//v_ = 1/2*(v[i][j]+v[i+1][j]);
		 //vn = 1/2*(v[i][j-1]+v[i+1][j-1]);
		//A = -1*(((ru2[i+1][j]-ru2[i-1][j])/(2*dx))+((ru[i][j+1]*v_-ru[i][j-1]*vn)/(2*dy)))+niu*((u[i+1][j]-2*u[i][j]+u[i-1][j])/(dx*dx)+(u[i][j+1]-2*u[i][j]+u[i][j-1])/(dy*dy));	
		ru_[i][j] = ru[i][j] -(dt/dx)*(p[i+1][j]-p[i][j]);	//ru在n+1时刻的数值 
			
		}
				
	}
	
	//更新变量 Ux 
	
	for(i=1;i<94;i++)
	{
		for(j=1;j<95;j++)
		{
			ru[i][j] = ru_[i][j];
			u[i][j] = ru[i][j]/rho;
			ru2[i][j] = ru[i][j]*u[i][j];
			// v值还没有更新 
			
		}
		
		
	 } 
	
	//左右边界 
	for(j=1;j<95;j++)
	{
		ru_[0][j] = 0; 
		//ru_[0][j] = 2*ru_[1][j]-ru_[2][j]; //边界条件外推 
		ru[0][j] = ru_[0][j];
		u[0][j] = ru[0][j]/rho; 
		ru2[0][j] = ru[0][j]*u[0][j];	
	
		ru_[94][j] = 0;
	//	ru_[94][j] = 2*ru_[93][j]-ru_[92][j]; //边界条件外推
		ru[94][j] = ru_[94][j];
		u[94][j] = ru[94][j]/rho;
		ru2[94][j] = ru[94][j]*u[94][j]; 
	
	
	
	
	
	 } 
	
	
	//上下边界
	for(i=0;i<=94;i++) 
	{
		ru_[i][0] = 0;
	//	ru_[i][0] = 2*ru_[i][1]-ru_[i][2]; 
		ru[i][0] = ru_[i][0];
		u[i][0] = ru[i][0]/rho; 
		ru2[i][0] = ru[i][0]*u[i][0];
		
		
		
		ru_[i][95] = 1;		
		ru[i][95] =ru_[i][95];
		u[i][95] = ru[i][95]/rho; 
		ru2[i][95] = ru[i][95]*u[i][95];
	}
	
	//========================================================================================//
	//计算 Uy 
	
	for(i=1;i<95;i++)
	{
		for(j=1;j<94;j++)
		{
	//	u_ = 1/2*(u[i][j+1]+u[i][j]); 
	  //  un = 1/2*(u[i-1][j+1]+u[i-1][j]);
	    
	    //B = -1*((rv[i+1][j]*u_-rv[i-1][j]*un)/(2*dx)+(rv2[i][j+1]-rv2[i][j-1])/(2*dy))+niu*((v[i+1][j]-2*v[i][j]+v[i-1][j])/(dx*dx)+(v[i][j+1]-2*v[i][j]+v[i][j-1])/(dy*dy));
		
		rv_[i][j] = rv[i][j]-(dt/dy)*(p[i][j+1]-p[i][j]);
			
		}
		
		
	} 
	
	//更新变量 Uy
	
	 for(i=1;i<95;i++)
	 {
	 	for(j=1;j<94;j++)
	 	{
	 		rv[i][j] = rv_[i][j];
			v[i][j] = rv[i][j]/rho;
			rv2[i][j] = rv[i][j]*v[i][j];
	 		
	 		
		 }
	 	
	 	
	  } 
	
	//左右边界 
	
	for(j=1;j<94;j++)
	{
		rv_[0][j] = 0;
		//rv_[0][j] = 2*rv_[1][j]-rv_[2][j];
		rv[0][j] = rv_[0][j];
		v[0][j] = rv[0][j]/rho;
		rv2[0][j] = rv[0][j]*v[0][j];
		
		
		rv_[95][j] = 0;
		//rv_[95][j] = 2*rv_[94][j]-rv_[93][j];
		rv[95][j] =rv_[95][j];
		v[95][j] = rv[95][j]/rho; 
		rv2[95][j] = rv[95][j]*v[95][j];
	 } 
	
	//上下边界
	for(i=0;i<=95;i++)
	{
	
	rv_[i][0] = 0;
	//rv_[i][0] = 2*rv_[i][1]-rv_[i][2];
	rv[i][0] = rv_[i][0];
	v[i][0] = rv[i][0]/rho;
	rv2[i][0] = rv[i][0]*v[i][0];
	
	rv_[i][94] = 0; 
	//rv_[i][94] = 2*rv_[i][93]-rv_[i][92];
	rv[i][94] = rv_[i][94];
	v[i][94] = rv[i][94]/rho;
	rv2[i][94] = rv[i][94]*v[i][94];
	
	
	
   }
}

	//对uT、vT更改
for(i=0;i<=94;i++)
{
	for(j=0;j<=94;j++)
	{
			
	uT[i][j] = u[i][j]*T[i][j];
	vT[i][j] = v[i][j]*T[i][j];
	
	T[i][94] = 1;
	uT[i][94] = u[i][94]*T[i][94];
	vT[i][94] = v[i][94]*T[i][94];
		
	}
	
		}
	
	
	for(i=1;i<94;i++)
{
	for(j=1;j<94;j++)
	{
		
		C =   (uT[i+1][j]-uT[i-1][j])/(2*dx);
	    D =   (vT[i][j+1]-vT[i][j-1])/(2*dy);
	    E = 0.1*(T[i+1][j]-2*T[i][j]+T[i-1][j])/(dx*dx)+0.1*(T[i][j+1]-2*T[i][j]+T[i][j-1])/(dy*dy);	
		T1[i][j] = dt*(-1*(C+D)+E)+T[i][j];
		uT[i][j] = u[i][j]*T1[i][j];
		vT[i][j] = v[i][j]*T1[i][j];
 	}
	
	
	
}

//左右边界

for(j=1;j<94;j++)
{
	T1[0][j] = 2*T1[1][j]-T1[2][j]; 
	uT[0][j] = u[0][j]*T1[0][j];
	vT[0][j] = v[0][j]*T1[0][j];
	
	T1[94][j] = 2*T1[93][j]-T1[92][j];
	uT[94][j] = u[94][j]*T1[94][j];
	vT[94][j] = v[94][j]*T1[94][j];
	
 } 
 
//上下边界

for(i=0;i<=94;i++)
{
	
	T1[i][0] = 2*T1[i][1]-T1[i][2];
	uT[i][0] = u[i][0]*T1[i][0];
	vT[i][0] = v[i][0]*T1[i][0];
	
	T1[i][94] = 1; 
	uT[i][94] = u[i][94]*T1[i][94];
	vT[i][94] = v[i][94]*T1[i][94];
	
 } 

	//更新宏观量 T
	
	
	for(i =0;i<=94;i++)
	{
		for(j=0;j<=94;j++)
		{
			
				T[i][j] = T1[i][j];
		}
		
		
	
	 } 
	 
	
	
	
	
	
	
	
	
	
	
	
	
	
}
	
	
	ostringstream name;
	  name<<"cavity_"<<6<<".dat";
	  ofstream out(name.str().c_str());
	  out<< "Title= \"LBM Lid Driven Flow\"\n" << "VARIABLES=\"X\",\"Y\",\"u\",\"v\",\"p\",\"T\"\n" << "ZONE T=\"BOX\",I=" << 95 << ",J=" << 95 << ",F=POINT" << endl;
	  for(int j=0;j<=94;j++)
	     for(int i=0;i<=94;i++)
	     {
	     	
			 
			 out<<double(i*dx)/Lx<<" "<<double(j*dx)/Ly<<" "<<u[i][j]<<" "<<v[i][j]<<" "<<p[i][j]<<" "<<T[i][j]<<endl;
		       
		
		 }
	
	
	
	
	
	return 0; 
 } 
  
