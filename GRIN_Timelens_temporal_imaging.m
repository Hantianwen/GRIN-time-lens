clc;
clear all;
%% Time windows
global Tmax Tmin delta_t;
Tmin=-5e-10; Tmax=5e-10; delta_t=1e-13; t=Tmin:delta_t:Tmax;
Ns=floor((Tmax-Tmin)/delta_t+1); df=1/(Tmax-Tmin);  k=-(Ns-1)/2:(Ns-1)/2;  f=fftshift(-k*df);
%% Input dispersion
beta2_DCF=20e-27;                                 % GVD of DCF (unit: s^2/m)
L_DCF=0.5*1e3;                                    % Length of DCF (unit:m)
phi_input=beta2_DCF*L_DCF; 
%% GRIN time lens parameters
gama=20e-3;                                     % Nonlinear coefficients (unit: W^-1*m^-1)
beta2_GRIN=-20e-27;                             % GVD of GRIN time lens (unit: s^2/m)
P_pump=1;                                       % Pump peak power (unit: W)
t_pump=40e-12;                                  % Half duration of pump (unit: s)
alfa=sqrt(-4*beta2_GRIN*gama*P_pump/t_pump^2);
Z_T=2*pi/alfa;                                  % Self-imaging period (unit: m)
L_GRIN=Z_T/4;                                   % Length of GRIN time lens (unit:m)
%% Output dispersion
phi_output=(beta2_GRIN/alfa*sin(alfa*L_GRIN)+phi_input*cos(alfa*L_GRIN))/(-1*cos(alfa*L_GRIN)+alfa/beta2_GRIN*phi_input*sin(alfa*L_GRIN));                                 
%% Temporal magnification factor
M_t=cos(alfa*L_GRIN)-alfa/beta2_GRIN*phi_output*sin(alfa*L_GRIN);
%% Input signal
t0=2e-12;                                       % Input pulse width (unit: s)
A_input=exp(-t.^2/(2*t0^2))+0.5*exp(-(t-5*t0).^2/(2*t0^2));
IA_input=abs(A_input).^2;
fA_input=fft(A_input);
figure
plot(t/1e-12,IA_input,'r');
xlim([-50,50]);ylim([0,1.2]);
xlabel('Time (ps)');ylabel('Amplitude');
title('Input temporal waveform');
%% Propagation through input dispersion
fA_incident=fA_input.*exp(1i*beta2_DCF*L_DCF*(2*pi*f).^2/2);
A_incident=ifft(fA_incident);
%% Propagation through GRIN time lens
Z=L_GRIN;                                           % GRIN time lens length                         
num=4e2;                                            % Simulation step number
dZ=Z/num;                                           % Step length                          
ufft=fft(A_incident);
for mm=2:num+1
    mm
    % Linear propagation dZ/2
    linear_operator=1i*beta2_GRIN*(2*pi*f).^2/2;                  % Dispersion opertor
    Dispersion=exp(linear_operator*dZ/2);     
    A_half=ifft(Dispersion.*ufft);
    % Nonlinear propagation dZ
    I_pump=1-t.^2/(t_pump^2);                                  % Pump pulse            
    XPM_opertor=1i*2*gama*P_pump*I_pump;                                              
    u1=A_half.*exp(XPM_opertor*dZ).*rectpuls((t)/(2*t_pump));  % XPM opertor
    % linear propagation dZ/2
    ufft=Dispersion.*fft(u1);
    A_output(mm,:)=ifft(ufft);
    I_output(mm,:)=abs(A_output(mm,:)).^2;
    fA_output(mm,:)=ufft;
    IfA_out(mm,:)=abs(fA_output(mm,:)).^2;
end
%% Propagation through output dispersion
A_out1=A_output(mm,:);
fA_out1=fA_output(mm,:);
fA_image=fA_out1.*exp(1i*phi_output*(2*pi*f).^2/2);
A_image=ifft(fA_image);
IA_image=abs(A_image).^2;
%% Output image
figure
plot(t/1e-12,IA_image,'r-');
hold on
plot(t/1e-12,IA_input,'b--');
xlim([-50,50]);ylim([0,2.2]);
xlabel('Time (ps)');ylabel('Amplitude')
legend('output','input');
title('Input signal and output image');