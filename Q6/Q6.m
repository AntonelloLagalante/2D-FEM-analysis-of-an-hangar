%% Question 6

close all
clear all
clc

load('C:\Users\anton\OneDrive - Politecnico di Milano\Università\Magistrale\Dynamics of mechanical system\Yearwork\YEARWORK_V2\Q1\hangar_mkr')

ndof=72;
ndoc=6;
ntot=ndof+ndoc;

MFF=M(1:ndof,1:ndof);
CFF=R(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);

MFC=M(1:ndof,ndof+1:ntot);
CFC=R(1:ndof,ndof+1:ntot);
KFC=K(1:ndof,ndof+1:ntot);

MCF=M(ndof+1:ntot,1:ndof);
CCF=R(ndof+1:ntot,1:ndof);
KCF=K(ndof+1:ntot,1:ndof);

MCC=M(ndof+1:ntot,ndof+1:ntot);
CCC=R(ndof+1:ntot,ndof+1:ntot);
KCC=K(ndof+1:ntot,ndof+1:ntot);

dof_xB=idb(18,1);
dof_yB=idb(18,2);

dof_xC=idb(9,1);
dof_yC=idb(9,2);

dof_xA=idb(5,1);
dof_yA=idb(5,2);

epsilon=0.005;
ms=5;

vect_force=zeros(ndof,1);
vect_force(dof_yA,1)=1;
vect_force(dof_xA,1)=1;

vect_f=0:0.01:20;

for k=1:length(vect_f)
    omega=vect_f(k)*2*pi; 
    A=-omega^2*MFF+i*omega*CFF+KFF;
    vect_f0=vect_force*ms*epsilon*omega^2;
    x0=A\vect_f0;

    xB=x0(dof_xB);
    xBdd=-omega^2*xB;
    yB=x0(dof_yB);
    yBdd=-omega^2*yB;
    xC=x0(dof_xC);
    xCdd=-omega^2*xC;
    yC=x0(dof_yC);
    yCdd=-omega^2*yC;

    mod1(k)=abs(xBdd);
    phase1(k)=angle(xBdd);
    mod2(k)=abs(yBdd);
    phase2(k)=angle(yBdd);
    mod3(k)=abs(xCdd);
    phase3(k)=angle(xCdd);
    mod4(k)=abs(yCdd);
    phase4(k)=angle(yCdd);

end
figure;
subplot 211;plot(vect_f,mod1);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('xBdd')
subplot 212;plot(vect_f,phase1);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')
figure;
subplot 211;plot(vect_f,mod2);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('yBdd')
subplot 212;plot(vect_f,phase2);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')

figure;
subplot 211;plot(vect_f,mod3);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('xCdd')
subplot 212;plot(vect_f,phase3);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')
figure;
subplot 211;plot(vect_f,mod4);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('yCdd')
subplot 212;plot(vect_f,phase4);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')