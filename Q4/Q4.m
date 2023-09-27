%% Question 4

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


% Frequency response function

%row = number of the node (1st coloumn: x, 2nd coloumn: x, 3rd coloumn: theta) 
xC0=zeros(ndoc,1);
doc_x1=idb(1,1);
doc_x26=idb(26,1);
xC0(1)=1;
xC0(4)=1;

i=sqrt(-1);
vect_f=0:0.01:20;

for j=1:length(vect_f)
    ome=2*pi*vect_f(j);
    
    A=-ome^2*MFF+i*ome*CFF+KFF;
    B=-ome^2*MFC+i*ome*CFC+KFC;
    vect_x0=-inv(A)*B*xC0;
    
    C=-ome^2*MCF+i*ome*CCF+KCF;
    D=-ome^2*MCC+i*ome*CCC+KCC;
    
    QC=0;

    vect_R=C*vect_x0+D*xC0-QC; 

    H=vect_R(1);
    M=vect_R(3);

    mod1(j)=abs(H);
    phase1(j)=angle(H);

    mod2(j)=abs(M);
    phase2(j)=angle(M);
end

figure;
subplot 211;plot(vect_f,mod1);grid;xlabel('Freq. [Hz]');ylabel('Amp. [N/m]');title('Horizontal reaction force in the clamp of node 1')
subplot 212;plot(vect_f,phase1*180/pi);grid;xlabel('Freq. [Hz]');ylabel('Phase [Deg]')

figure;
subplot 211;plot(vect_f,mod2);grid;xlabel('Freq. [Hz]');ylabel('Amp. [Nm/m]');title('Clamping moment in the clamp of node 1')
subplot 212;plot(vect_f,phase2*180/pi);grid;xlabel('Freq. [Hz]');ylabel('Phase [Deg]')