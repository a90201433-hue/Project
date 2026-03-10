#include "MaderSolver.h"
#include "TransformValues.h"
#include <cmath>

Field omega;
Field pressure;
Field temperature;

Field q1,q2,q3,q4;

Field u_tilde;
Field v_tilde;

Field dMass;
Field dEnergy;
Field dMomentumX;
Field dMomentumY;
Field dOmega;

int NX, NY;

extern int Nx;
extern int Ny;
extern int fict;

extern double gamma_gas;
extern double VISC;
extern double Z_freq;
extern double E_act;
extern double Rgas;
extern double MINWT;
extern double GASW;
extern double MINGRHO;
extern double M_molar;


void InitializeMaderSolver(int Nx,int Ny)
{
    NX=Nx+2*fict;
    NY=Ny+2*fict;

    omega.resize(NX,std::vector<State>(NY));
    pressure.resize(NX,std::vector<State>(NY));
    temperature.resize(NX,std::vector<State>(NY));

    q1.resize(NX,std::vector<State>(NY));
    q2.resize(NX,std::vector<State>(NY));
    q3.resize(NX,std::vector<State>(NY));
    q4.resize(NX,std::vector<State>(NY));

    u_tilde.resize(NX,std::vector<State>(NY));
    v_tilde.resize(NX,std::vector<State>(NY));

    dMass.resize(NX,std::vector<State>(NY));
    dEnergy.resize(NX,std::vector<State>(NY));
    dMomentumX.resize(NX,std::vector<State>(NY));
    dMomentumY.resize(NX,std::vector<State>(NY));
    dOmega.resize(NX,std::vector<State>(NY));
}


void ComputeEquationOfState(Field& U)
{
    for(int i=fict;i<Nx+fict-1;i++)
    for(int j=fict;j<Ny+fict-1;j++)
    {
        double rho=std::max(U[i][j][0],MINGRHO);

        double u=U[i][j][1]/rho;
        double v=U[i][j][2]/rho;

        double E=U[i][j][3]/rho;

        double I=E-0.5*(u*u+v*v);

        double P=(gamma_gas-1.0)*rho*I;

        if(P<1e-10) P=1e-10;

        pressure[i][j][0]=P;

        double R=Rgas/M_molar;

        temperature[i][j][0]=P/(rho*R);
    }
}


void UpdateChemicalReaction(double dt)
{
    for(int i=fict;i<Nx+fict-1;i++)
    for(int j=fict;j<Ny+fict-1;j++)
    {
        double T=temperature[i][j][0];

        if(T>MINWT && omega[i][j][0]>GASW)
        {
            double exponent=-E_act/(Rgas*T);

            if(exponent<-700) exponent=-700;

            double rate=Z_freq*omega[i][j][0]*exp(exponent);

            omega[i][j][0]-=dt*rate;

            if(omega[i][j][0]<GASW)
                omega[i][j][0]=0;
        }
    }
}


// ИСПРАВЛЕНИЕ 1: знак Q1 и Q2 изменён на правильный.
// По документу (стр. 8):
//   Q1(i,j) = VISC · ρ · (V(i,j) − V(i,j+1))  если V(i,j) ≥ V(i,j+1)
//   Q2(i,j) = VISC · ρ · (U(i,j) − U(i+1,j))  если U(i,j) ≥ U(i+1,j)
// Вязкость всегда неотрицательна — она диссипирует энергию при сжатии.
void ComputeArtificialViscosity(Field& U)
{
    for(int i=fict;i<Nx+fict-1;i++)
    for(int j=fict;j<Ny+fict-1;j++)
    {
        double rho=U[i][j][0];

        double v=U[i][j][2]/rho;
        double v_up=U[i][j+1][2]/U[i][j+1][0];

        // Сжатие в Z-направлении: v >= v_up
        if(v_up-v<0)
            q1[i][j][0]=VISC*rho*(v-v_up);   // ✓ было: (v_up-v) → отрицательное
        else
            q1[i][j][0]=0;

        double u=U[i][j][1]/rho;
        double u_r=U[i+1][j][1]/U[i+1][j][0];

        // Сжатие в R-направлении: u >= u_r
        if(u_r-u<0)
            q2[i][j][0]=VISC*rho*(u-u_r);    // ✓ было: (u_r-u) → отрицательное
        else
            q2[i][j][0]=0;

        q3[i][j][0]=q1[i][j-1][0];
        q4[i][j][0]=q2[i-1][j][0];
    }
}


void UpdateVelocities(Field& U,double dt,double dx,double dy)
{
    for(int i=fict;i<Nx+fict-1;i++)
    for(int j=fict;j<Ny+fict-1;j++)
    {
        double rho=U[i][j][0];

        double u=U[i][j][1]/rho;
        double v=U[i][j][2]/rho;

        double P1=pressure[i][j-1][0];
        double P3=pressure[i][j+1][0];

        double P2=pressure[i-1][j][0];
        double P4=pressure[i+1][j][0];

        v_tilde[i][j][0]=
            v-dt/(rho*dy)*(P3-P1+q3[i][j][0]-q1[i][j][0]);

        u_tilde[i][j][0]=
            u-dt/(rho*dx)*(P4-P2+q4[i][j][0]-q2[i][j][0]);
    }
}


void UpdateInternalEnergy(Field& U,double dt,double dx,double dy)
{
    for(int i=fict;i<Nx+fict-1;i++)
    for(int j=fict;j<Ny+fict-1;j++)
    {
        double rho=U[i][j][0];

        double u=U[i][j][1]/rho;
        double v=U[i][j][2]/rho;

        double uL=U[i-1][j][1]/U[i-1][j][0];
        double uR=U[i+1][j][1]/U[i+1][j][0];

        double vD=U[i][j-1][2]/U[i][j-1][0];
        double vU=U[i][j+1][2]/U[i][j+1][0];

        double U1=uL+u_tilde[i-1][j][0];
        double U2=uR+u_tilde[i+1][j][0];

        double V1=vD+v_tilde[i][j-1][0];
        double V2=vU+v_tilde[i][j+1][0];

        double T3=u+u_tilde[i][j][0];
        double T1=v+v_tilde[i][j][0];

        double P=pressure[i][j][0];

        double I=(U[i][j][3]/rho)-0.5*(u*u+v*v);

        double I_new=
            I-(dt/(4*rho))*(
                (P/dx)*(U2-U1)
              +(q4[i][j][0]/dx)*(U2-T3)
              +(q2[i][j][0]/dx)*(T3-U1)
              +(P/dy)*(V2-V1)
              +(q3[i][j][0]/dy)*(V2-T1)
              +(q1[i][j][0]/dy)*(T1-V1)
            );

        double E=
            I_new
            +0.5*(u_tilde[i][j][0]*u_tilde[i][j][0]
            +v_tilde[i][j][0]*v_tilde[i][j][0]);

        U[i][j][3]=rho*E;
    }
}


// ИСПРАВЛЕНИЕ 2: добавлен перенос по Y (Z-направление).
// Ранее TransportMass обрабатывала только грани в R-направлении (по i).
// По документу (стр. 18-20) необходимо обойти все 4 грани — и по R, и по Z.
void TransportMass(Field& U,double dt,double dx,double dy)
{
    for(int i=0;i<NX;i++)
    for(int j=0;j<NY;j++)
    {
        dMass[i][j][0]=0;
        dEnergy[i][j][0]=0;
        dMomentumX[i][j][0]=0;
        dMomentumY[i][j][0]=0;
        dOmega[i][j][0]=0;
    }

    // --- Перенос по R (X) --- между (i-1,j) и (i,j)
    for(int i=fict+1;i<Nx+fict-1;i++)
    for(int j=fict;j<Ny+fict-1;j++)
    {
        double u_face=
        0.5*(u_tilde[i-1][j][0]+u_tilde[i][j][0]);

        double alpha=u_face*dt/dx;

        if(alpha>1) alpha=1;
        if(alpha<-1) alpha=-1;

        int donor=(alpha>=0)?i-1:i;
        int acc=(alpha>=0)?i:i-1;

        double rho=U[donor][j][0];

        double dm=rho*fabs(alpha);

        dMass[acc][j][0]+=dm;
        dMass[donor][j][0]-=dm;

        double momx=U[donor][j][1];
        double momy=U[donor][j][2];
        double rhoE=U[donor][j][3];

        dMomentumX[acc][j][0]+=momx*fabs(alpha);
        dMomentumX[donor][j][0]-=momx*fabs(alpha);

        dMomentumY[acc][j][0]+=momy*fabs(alpha);
        dMomentumY[donor][j][0]-=momy*fabs(alpha);

        dEnergy[acc][j][0]+=rhoE*fabs(alpha);
        dEnergy[donor][j][0]-=rhoE*fabs(alpha);

        double w=omega[donor][j][0];

        dOmega[acc][j][0]+=w*dm;
        dOmega[donor][j][0]-=w*dm;
    }

    // --- Перенос по Z (Y) --- между (i,j-1) и (i,j)
    for(int i=fict;i<Nx+fict-1;i++)
    for(int j=fict+1;j<Ny+fict-1;j++)
    {
        double v_face=
        0.5*(v_tilde[i][j-1][0]+v_tilde[i][j][0]);

        double beta=v_face*dt/dy;

        if(beta>1) beta=1;
        if(beta<-1) beta=-1;

        int donor=(beta>=0)?j-1:j;
        int acc=(beta>=0)?j:j-1;

        double rho=U[i][donor][0];

        double dm=rho*fabs(beta);

        dMass[i][acc][0]+=dm;
        dMass[i][donor][0]-=dm;

        double momx=U[i][donor][1];
        double momy=U[i][donor][2];
        double rhoE=U[i][donor][3];

        dMomentumX[i][acc][0]+=momx*fabs(beta);
        dMomentumX[i][donor][0]-=momx*fabs(beta);

        dMomentumY[i][acc][0]+=momy*fabs(beta);
        dMomentumY[i][donor][0]-=momy*fabs(beta);

        dEnergy[i][acc][0]+=rhoE*fabs(beta);
        dEnergy[i][donor][0]-=rhoE*fabs(beta);

        double w=omega[i][donor][0];

        dOmega[i][acc][0]+=w*dm;
        dOmega[i][donor][0]-=w*dm;
    }
}


void ApplyTransport(Field& U)
{
    for(int i=fict;i<Nx+fict-1;i++)
    for(int j=fict;j<Ny+fict-1;j++)
    {
        double rho_old=U[i][j][0];

        double rho_new=rho_old+dMass[i][j][0];

        if(rho_new<MINGRHO) rho_new=MINGRHO;

        double momx=
        U[i][j][1]+dMomentumX[i][j][0];

        double momy=
        U[i][j][2]+dMomentumY[i][j][0];

        double rhoE=
        U[i][j][3]+dEnergy[i][j][0];

        double u=momx/rho_new;
        double v=momy/rho_new;

        double E=rhoE/rho_new;

        double I=E-0.5*(u*u+v*v);

        U[i][j][0]=rho_new;
        U[i][j][1]=momx;
        U[i][j][2]=momy;
        U[i][j][3]=rho_new*(I+0.5*(u*u+v*v));

        if(rho_new>MINGRHO)
        {
            omega[i][j][0]=
            (omega[i][j][0]*rho_old+dOmega[i][j][0])/rho_new;
        }
    }
}


void MaderTimeStep(
    Field& W,
    Field& W_new,
    double dt,
    double dx,
    double dy)
{
    Field U;

    ConvertWtoU(W,U,0);

    ComputeEquationOfState(U);

    UpdateChemicalReaction(dt);

    ComputeArtificialViscosity(U);

    UpdateVelocities(U,dt,dx,dy);

    UpdateInternalEnergy(U,dt,dx,dy);

    TransportMass(U,dt,dx,dy);

    ApplyTransport(U);

    ComputeEquationOfState(U);

    ConvertUtoW(W_new,U,0);
}
