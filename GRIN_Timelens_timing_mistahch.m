clc;
clear all;
%% Time windows
global Tmax Tmin delta_t;
Tmin=-5e-10; Tmax=5e-10; delta_t=1e-14; t=Tmin:delta_t:Tmax;
Ns=floor((Tmax-Tmin)/delta_t+1); df=1/(Tmax-Tmin);  k=-(Ns-1)/2:(Ns-1)/2;  f=fftshift(-k*df);
%% Input dispersion
beta2_DCF=20e-27;                               % GVD of DCF (unit: s^2/m)
L_DCF=5e3;                                      % Length of DCF (unit:m)
%% GRIN time-lens parameters
gama=20e-3;                                     % Nonlinear coefficients (unit: W^-1*m^-1)
beta2_GRIN=-20e-27;                             % GVD of GRIN time lens (unit: s^2/m)
P_pump=1;                                       % Pump peak power (unit: W)
t_pump=40e-12;                                  % Half duration of pump (unit: s)
alfa=sqrt(-4*beta2_GRIN*gama*P_pump/t_pump^2);
Z_T=2*pi/alfa;                                  % Self-imaging period (unit: m)
L_GRIN=Z_T;                                     % Length of GRIN time lens (unit:m)
%% Output dispersion
beta2_SMF=-20e-27;                              % GVD of SMF (unit:s^2/m)
L_SMF=5e3;                                      % Length of SMF (unit:m)
%% Input signal
t0=2e-12;                                       % Âö¿í (unit:ps)
A_input=exp(-t.^2/(2*t0^2));
fA_input=fft(A_input);
figure
plot(t/1e-12,abs(A_input),'r');
xlim([-100,100]);ylim([0,1.2]);
xlabel('Time (ps)');ylabel('Amplitude');
title('Input signal');
%% Propagation through input dispersion
fA_incident=fA_input.*exp(1i*beta2_DCF*L_DCF*(2*pi*f).^2/2);
A_incident=ifft(fA_incident);
figure
plot(t/1e-12,abs(A_incident),'r-');
xlim([-100,100]);ylim([0,1.2]);
xlabel('Time (ps)');ylabel('Amplitude');
title('Temporal waveform at the incident plane of the lens');
%% Propagation through GRIN time lens
Z=L_GRIN;                                           % GRIN time lens length                         
num=4e2;                                            % Simulation step number
dZ=Z/num;                                           % Step length                          
ufft=fft(A_incident);
for mm=2:num+1
    mm
    % Linear propagation dZ/2
    linear_operator=1i*beta2_GRIN*(2*pi*f).^2/2;                  
    Dispersion=exp(linear_operator*dZ/2);     
    A_half=ifft(Dispersion.*ufft);
    % Nonlinear propagation dZ
    I_pump=1-(t-20e-12).^2/(t_pump^2);                                           
    XPM_opertor=1i*2*gama*P_pump*I_pump;                                              
    u1=A_half.*exp(XPM_opertor*dZ).*rectpuls((t-20e-12)/(2*t_pump));
    % linear propagation dZ/2
    ufft=Dispersion.*fft(u1);
    A_output(mm,:)=ifft(ufft);
    I_output(mm,:)=abs(A_output(mm,:)).^2;
    fA_output(mm,:)=ufft;
    IfA_output(mm,:)=abs(fA_output(mm,:)).^2;
end
A_output(1,:)=A_incident;
I_output(1,:)=abs(A_incident).^2;
fA_output(1,:)=fA_incident;
IfA_output(1,:)=abs(fA_incident).^2;
z=0:dZ:Z;
A_outputmax=max(max(abs(A_output)));
figure
grid off
mesh(t/1e-12,z/1e3,abs(A_output)/A_outputmax);
view(0,90);
xlim([-60,60]);ylim([0,Z_T/1e3]);
xlabel('Time (ps)');ylabel('Propagation distance(km)');
colorbar;
%% Temporal waveforms at the input and output planes of the lens
figure
plot(t/1e-12,abs(A_output(mm,:)),'r-');
hold on
plot(t/1e-12,abs(A_incident),'b--');
xlim([-100,100]);ylim([0,1.2]);
xlabel('Time (ps)');ylabel('Amplitude');
legend('after GRIN','before GRIN');
title('Temporal waveforms at the input and output planes of the lens');
%% Propagation through SMF
A_out1=A_output(mm,:);
fA_out1=fA_output(mm,:);
D0=linspace(0,10e3,200);
for nj=1:length(D0);
    nj
    z_SMF=D0(nj);
    fAout=fA_out1.*exp(1i*beta2_SMF*z_SMF*(2*pi*f).^2/2);
    Aout(nj,:)=ifft(fAout);
    Iout(nj,:)=abs(Aout(nj,:));
end
Aoutmax=max(max(abs(Aout)));
figure
grid off
mesh(t/1e-12,D0*beta2_SMF/1e-24,Iout/Aoutmax); 
view(0,90);
xlim([-60,60]);ylim([-150,0]);
xlabel('Time (ps)');ylabel('GDD(ps^2)');
set(gca,'ydir','reverse');
colorbar;
%% Input signal and output images
fA_image=fA_out1.*exp(1i*beta2_SMF*L_SMF*(2*pi*f).^2/2);
A_image=ifft(fA_image);
figure
plot(t/1e-12,abs(A_input),'k.');
hold on
plot(t/1e-12,abs(Aout(100,:)),'r--');
hold on
plot(t/1e-12,abs(Aout(120,:)),'b-');
xlim([-60,60]);ylim([0,1.2]);
xlabel('Time (ps)');ylabel('Amplitude');
legend('Input pulse','phi_out=-100ps^2','phi_out=-120ps^2');
title('Input signal and output image');