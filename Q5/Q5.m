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


% Triangular waveform
T = 0.3; % period
A = 2000; % amplitude
t = 0:1e-3:1.5; % time vector
ty=0:1e-3:T;
y = zeros(1,length(ty));
y = A*(2/T)*abs(mod(ty-T/4,T) - T/2) - A/2;

figure()
plot(ty, y);
xlabel('Time (seconds)');
ylabel('Amplitude');


%Fourier transform
fftout=fft(y);
N=length(y);
df=1/T;
fmax=(N/2-1)*df;
vett_freq=0:df:fmax;
modf(1)=1/N*abs(fftout(1));
modf(2:N/2)=2/N*abs(fftout(2:N/2));
fasf(1:N/2)=angle(fftout(1:N/2));

figure
subplot 211;bar(vett_freq,modf);
subplot 212;plot(vett_freq,fasf);

dof_yD=idb(7,2); %vertical displacement for point D --> node 7

ome0=2*pi/T; %foundamental frequency
output_yD=zeros(1,length(t));
output_yDdd=zeros(1,length(t));

vect_force=zeros(ndof,1);
vect_force(dof_yD,1)=1;

i=sqrt(-1); %the imaginary operator i has to be defined in the script

for k = 1:N/2
   omega(k)=2*pi*vett_freq(k);
   if omega(k)<=(pi*fmax*2)
    A=-omega(k)^2*MFF+i*omega(k)*CFF+KFF;
    x0 = A\(modf(k)*vect_force);
    yD=x0(dof_yD);

    yDdd=-omega(k)^2*yD;
    output_yD=output_yD+abs(yD)*cos(omega(k)*t+angle(yD));
    output_yDdd=output_yDdd+abs(yDdd)*cos(omega(k)*t+angle(yDdd));
   end
end

figure;
plot(t,output_yD,'LineWidth',1)
title("Vertical displacement of point D");
figure;
plot(t,output_yDdd,'LineWidth',1)
title("Vertical acceleration of point D");


%Resonance
f1=2.2286;
f2=6.5068;
f3=17.4237;
f4=20.766;

T1=1/f1
T2=1/f2
T3=1/f3
T4=1/f4