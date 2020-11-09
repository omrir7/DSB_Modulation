% DSB-LC - Amplitude modulation - full channel model
% Author : Omri Rafaeli
%
%                   Table of Contents:
%   Chapter:
%------------------------------------------------------------------------
%   1.2.1.1 - Ploting Data and carrier signals
%   1.2.1.2 - Creating the DSB Modulated signal and its Fourier Transform
%   1.2.2.1 - Generating the Noise vectors 
%   1.2.2.2 - Calculating the recieved signals and its fourier transform
%   1.2.3.1 - Generating a Band Pass Filter 
%   1.2.3.2 - Calculating the filtered signal as convolution with the BPF
%   1.2.3.3 - Extracting the envelope from the filtered signal  

close all;

%--------------------------------------------------------------------------------------------------------------------------

% Defining constants
A_1=1;          % first cos Amp
A_2=2;          % second cos Amp
A_c=2;          % carrier Amp
f_1=1;          % first cos freq
f_2=2.5;        % second cos freq
f_c=100;        % carrier freq
kAM=0.2;        % Modulation index
f_s=1000;       % sampling freq
T_s=1/f_s;      % sampling time
t=-3:T_s:3;     % Generating the time vector
t=setdiff(t,0); % take out 0

%--------------------------------------------------------------------------------------------------------------------------
%Modulator
%1.2.1.1 - Ploting Data and carrier signals

v_m_t=A_1*cos(2*pi*f_1*t)+A_2*cos(2*pi*f_2*t);  % data signal - Freq=2.5Hz
v_c_t=A_c*cos(2*pi*f_c*t);                      % carrier signal - Freq=100Hz

%Ploting Vm(t) - The Data Signal
figure;
subplot(2,1,1);
plot(t,v_m_t);
xlabel('Time [sec]','FontWeight','bold','Fontname', 'Arial');
ylabel('v_{m}(t)','FontWeight','bold','Fontname', 'Arial');
title('v_{m}(t) vs. Time');
grid;

%Ploting Vm(t) - The Carrier Signal
hold on;
subplot(2,1,2);
plot(t,v_c_t);
xlabel('Time [sec]','FontWeight','bold','Fontname', 'Arial');
ylabel('v_{c}(t)','FontWeight','bold','Fontname', 'Arial');
title('v_{c}(t) vs. Time');
grid;
hold off;
%--------------------------------------------------------------------------------------------------------------------------
%1.2.1.2 - Creating the DSB Modulated signal and its Fourier Transform

 % (1+kAM*v_m).*v_c - Amplitude Modulation DSB-LC
v_mod_t=DSB_LC_MOD(v_m_t,v_c_t,kAM);

%Ploting the modulated signal
figure;
plot(t,v_mod_t);
xlabel('Time [sec]','FontWeight','bold','Fontname', 'Arial');
ylabel('v_{M}(t)','FontWeight','bold','Fontname', 'Arial');
title('v_{M}(t) vs. Time');
grid;

%Fourier Transform for the modulated signal
v_mod_f=fftshift(fft(v_mod_t)/length(t)); % fourier transform fft shift shifts the transform to be centerd and length t normalizing the transform 

%generating the frequency vector
f=linspace(-f_s/2,f_s/2,length(t)); 

%Ploting the Fourier transform of the modulated signal
figure;
plot(f,v_mod_f);
xlabel('Frequency [Hz]','FontWeight','bold','Fontname', 'Arial');
ylabel('V_{M}(f)','FontWeight','bold','Fontname', 'Arial');
title('V_{M}(f) vs. Frequency');
grid;
%--------------------------------------------------------------------------------------------------------------------------

%1.2.2.1 - Generating the Noise vectors 
N_1=2; % variances
N_2=20;
N_3=200;

% generating a uniformed distributed vector with var=N1/2 and same length as time vector
z_1_t=sqrt(N_1/2)*randn(1,length(t)); 
z_2_t=sqrt(N_2/2)*randn(1,length(t));
z_3_t=sqrt(N_3/2)*randn(1,length(t));
%--------------------------------------------------------------------------------------------------------------------------

%1.2.2.2 - Calculating the recieved signals and its fourier transform

%the recieved signal is being calculated according to a Added gaussian
%noise channel
x_r_1_t=v_mod_t+z_1_t; % calculating the recieved signal 
x_r_2_t=v_mod_t+z_2_t;
x_r_3_t=v_mod_t+z_3_t;

% fourier transform of the recieved signal and normalizing it
X_r_1_f=fftshift(fft(x_r_1_t)/length(t)); 
X_r_2_f=fftshift(fft(x_r_2_t)/length(t));
X_r_3_f=fftshift(fft(x_r_3_t)/length(t));

%Ploting the 3 recieved signals in the frequency plane - x_r_i | 1<=i<=3
figure;
plot(f,X_r_3_f);
hold on;
plot(f,X_r_2_f);
hold on;
plot(f,X_r_1_f);
hold off;

xlabel('Frequency [Hz]','FontWeight','bold','Fontname', 'Arial');
ylabel('X_{r}(f)','FontWeight','bold','Fontname', 'Arial');
title('X_{r}(f) vs. Frequency: For Different Noise Variances');
legend('Noise Variance=200','Noise Variance=20','Noise Variance=2');
grid;
%--------------------------------------------------------------------------------------------------------------------------

%       Demodulator
%1.2.3.1 -  Calculating the filtered signal as convolution with the BPF

%the width of the window will include the entire bandwidth of the original
%data - max(f1,f2)
B=2*f_2;              

%sinc in the freq domain is rect width B , cos is shifting the rect from 0 to fc
%BPF_t=B*sinc(B*t).*cos(2*pi*f_c*t);           
BPF_t=4*pi*sin(t*2*f_2*pi)./(pi*t).*cos(2*pi*t*f_c)/1.07;
BPF_f=fftshift(fft(BPF_t)/length(t));       

%Ploting the Band Pass Filter - not ideal rect but fast ascend and fast descend
figure;
plot(f,BPF_f);
xlabel('Frequency [Hz]','FontWeight','bold','Fontname', 'Arial');
ylabel('BPF(f)','FontWeight','bold','Fontname', 'Arial');
title('Band Pass Filter vs. Frequency');

grid;
%--------------------------------------------------------------------------------------------------------------------------
%1.2.3.2 - Calculating the filtered signal as convolution with the BPF

%Generating the filtered signal by covoluting the recieved signal with the
%BPF (output of a linear time invarient system)

x_1_filtered_t=conv(x_r_1_t,BPF_t,'same')/length(t); 
x_2_filtered_t=conv(x_r_2_t,BPF_t,'same')/length(t);
x_3_filtered_t=conv(x_r_3_t,BPF_t,'same')/length(t);

% Fourier transform of the filtered signals
x_1_filtered_f=fftshift(fft(x_1_filtered_t)/length(t));
x_2_filtered_f=fftshift(fft(x_2_filtered_t)/length(t));
x_3_filtered_f=fftshift(fft(x_3_filtered_t)/length(t));

%ploting the fourier transform of the filtered signals
figure;
plot(f,x_3_filtered_f);
hold on;
plot(f,x_2_filtered_f);
hold on;
plot(f,x_1_filtered_f);
hold off;

xlabel('Frequency [Hz]','FontWeight','bold','Fontname', 'Arial');
ylabel('x_{filtered}(f)','FontWeight','bold','Fontname', 'Arial');
title('x_{filtered}(f) vs. Frequency: For Different Noise Variances');
legend('Noise Variance=200','Noise Variance=20','Noise Variance=2');
grid;
%--------------------------------------------------------------------------------------------------------------------------

%1.2.3.3 - Extracting the envelope from the filtered signal

% extract envelopes from the filtered signal in order to get the data signal
[upper_env_1,lower_env_1]=envelope(x_1_filtered_t); 
[upper_env_2,lower_env_2]=envelope(x_2_filtered_t);
[upper_env_3,lower_env_3]=envelope(x_3_filtered_t);

%ploting the extracted envelope from the recieved filtered signals - can be
%seen that low noise varience -> better reconstruction of the data signal
figure;
plot(t,upper_env_3);
hold on;
plot(t,upper_env_2);
hold on;
plot(t,upper_env_1);
hold on;
plot(t,v_m_t);
hold off;

xlabel('Time [sec]','FontWeight','bold','Fontname', 'Arial');
ylabel('x_{d}(t)','FontWeight','bold','Fontname', 'Arial');
title('x_{d}(t) vs. Time: For Different Noise Variances');
legend('Noise Variance=200','Noise Variance=20','Noise Variance=2','v_{m}(t)');
grid;
%--------------------------------------------------------------------------------------------------------------------------
% DSB-LC modulation function
%inputs: v_m -   Data signal
%        v_c -   Carrier Signal 
%        kAM -   Modulation index
%output: v_mod - Modulated signal

function v_mod = DSB_LC_MOD(v_m,v_c,kAM)
    v_mod = (1+kAM*v_m).*v_c;
end

