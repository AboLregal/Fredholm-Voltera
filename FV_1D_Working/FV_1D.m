% clc;clear;close all

function f=FV_1D(c)

%% Hell_2_Project

global N H1 H2 Y const1 const2 lamda1 lamda2
% global N
% N=10;


%% WK

for k=0:N
    for n=0:2:N
        if n==0 || n==N
            wz(n+1,k+1)=(2/N)*cos(n*k*pi/N)/(1-(n*n));
        else
            wz(n+1,k+1)=(4/N)*cos(n*k*pi/N)/(1-(n*n));
        end
    end
end

wk=sum(wz);

%% T's Hell

for i=0:N
    Si=cos(i*pi/N);
    for omega=1:N+1
        if omega==1
            T_L(omega,i+1)=1;
        elseif omega==2
            T_L(omega,i+1)=Si;
        else
            T_L(omega,i+1)=((2*Si)*T_L(omega-1,i+1))-T_L(omega-2,i+1);
        end
    end
end

for k=0:N
    Sk=(1+cos(k*pi/N));
    for omega=1:N+1
        if omega==1
            T_2(omega,k+1)=1;
        elseif omega==2
            T_2(omega,k+1)=Sk/2;
        else
            T_2(omega,k+1)=(Sk*T_2(omega-1,k+1))-T_2(omega-2,k+1);
        end
    end
end

for i=0:N
    for k=0:N
        Si=cos(i*pi/N);
        Sk=(1+cos(k*pi/N));
        for omega=1:N+1
            if omega==1
                T_1(omega,k+1,i+1)=1;
            elseif omega==2
                T_1(omega,k+1,i+1)=(Si*Sk)/2;
            else
                T_1(omega,k+1,i+1)=(Si*Sk*T_1(omega-1,k+1,i+1))-T_1(omega-2,k+1,i+1);
            end
        end
    end
    
end

%% Right Hand Calculations

% c=rand(1,N+1);
for i=1:N+1
    T_Z=T_1(:,:,i);
    Si=cos((i-1)*pi/N);
    y=Y(Si);
    Lamda1=Si*lamda1/2;
    Lamda2=lamda2/2;
    for k=0:N
        Sk=cos(k*pi/N);
        C1(1,k+1)=const1(Si,Sk);
        C2(1,k+1)=const2(Si,Sk);
    end
    for k=1:N+1
        T1x=T_Z(:,k);
        T2x=T_2(:,k);
        I1=Lamda1*C1(k)*H1(c,T1x);
        I2=Lamda2*C2(k)*H2(c,T2x);
        I(k)=wk(k)*(I1+I2);
    end
    Ix=sum([0.5*I(1) I(2:end-1) 0.5*I(end)]);
    R_H(1,i)=y+Ix;
end

%% L_H
for j=1:N+1
    T_L_J=T_L(:,j);
    R_L(1,j)=c*T_L_J;
end

%% Solving Target

f=R_L-R_H;



