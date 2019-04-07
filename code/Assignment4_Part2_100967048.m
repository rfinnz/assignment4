%% Assignment 4 Part 2

% In this part of the assignment, the transient response of the circuit was
% simulated

% Part a) Upon inspection of the circuit, this is a RLC circuit, due to the
% presence of capacitors, resistors and and inductors. DC transient response
% in RLC circuits can be observed when a sudden voltage or current is
% applied to it, which is what we will be doing in part 2

%Part b) In simple terms, RLC circuits experience "resonance" at a certain
%frequency. This frequency point is where the reactive inductance of the
%inductor equals the value of the capacitance reactance of the capacitor
%(xL = xC). At this frequency the gain of the circuit sharply rises, and is
%low at other frequencies. For this reason, this circuit has a very sharp
%bandpass response, and can be used as a filter or an amplifier functioning
%at a certain frequency.

% Definition of variables based on the components present in the circuit
R1 = 1;
G1 = 1/R1;
C = 0.25;
R2 = 2;
G2 = 1/R2;
L = 0.2;
R3 = 10;
G3 = 1/R3;
alpha = 100;
R4 = 0.1;
G4 = 1/R4;
RO = 1000;
GO = 1/RO;
Vin = 1;

% Define Matrices
C_Matrix = [0 0 0 0 0 0 0;
           -C C 0 0 0 0 0;
            0 0 -L 0 0 0 0;
            0 0 0 0 0 0 0;
            0 0 0 0 0 0 0;
            0 0 0 0 0 0 0;
            0 0 0 0 0 0 0;];

G_Matrix = [1 0 0 0 0 0 0;
           -G2 G1+G2 -1 0 0 0 0;
            0 1 0 -1 0 0 0;
            0 0 -1 G3 0 0 0;
            0 0 0 0 -alpha 1 0;
            0 0 0 G3 -1 0 0;
            0 0 0 0 0 -G4 G4+GO];

F_Matrix = [Vin;
             0;
             0;
             0;
             0;
             0;
             0;];

F0_Matrix = [Vin-Vin;       
                0;
                0;
                0;
                0;
                0;
                0;];

%d) we will be simulating the circuit for 1 second using 1000 steps
step = 1000;              


vol_1 = zeros(7, step);
vol_start = zeros(7, 1);
dt = 10^-3;
 %setting up the plot for the first input signal, a step that transitions
 %from 0 to 1 at 30 miliseconds
for i = 1:step

    if i < 30
        vol_1(:,i) = (C_Matrix./dt+G_Matrix)\(F0_Matrix+C_Matrix*vol_start/dt);
        
    elseif i == 30
        vol_1(:,i) = (C_Matrix./dt+G_Matrix)\(F_Matrix+C_Matrix*vol_start/dt);
        
    else
        vol_1(:,i) = (C_Matrix./dt+G_Matrix)\(F_Matrix+C_Matrix*vol_old/dt);
        
    end
    
    vol_old = vol_1(:, i);
    
end

figure(1)
plot(1:step, vol_1(7,:), 'g')
hold on
plot(1:step, vol_1(1,:), 'm')
title('Vin and Vo with a step that transitions 0-1 at 0.03s')
xlabel('Time in miliseconds')
ylabel('Voltage in volts')
grid on


vol_2 = zeros(7, step);
function_F = zeros(7,1);

%setting up the plot for the second input signal, a sin(2*pi*f*t) signal,
%at a frequency of 1/(30) 1/ms.

for i_2 = 1:step

    function_vol = sin(2*pi*(1/0.03)*i_2/step);
    function_F(1,1) = function_vol;
    if i_2 == 1
        vol_2(:,i_2) = (C_Matrix./dt+G_Matrix)\(function_F+C_Matrix*vol_start/dt);
    else
        vol_2(:,i_2) = (C_Matrix./dt+G_Matrix)\(function_F+C_Matrix*vol_old/dt);
    end
    vol_old = vol_2(:, i_2);
        
end

figure(2)
plot(1:step, vol_2(7,:), 'b')
hold on
plot(1:step, vol_2(1,:), 'g')
title('Vin and Vo with input signal function sin(2pift),f = 1/30 1/ms')
xlabel('Time in miliseconds')
ylabel('Voltage in volts')
grid on


%setting up the plot using a guassian pulse with mag=1, std dev = 30ms and
%delay of 60 ms

vol_3 = zeros(7, step);
Gaussian_F = zeros(7,1);

for i_3 = 1:step

    Guassian_vol = exp(-1/2*((i_3/step-0.06)/(0.03))^2);
    Gaussian_F(1,1) = Guassian_vol;
    
    if i_3 == 1
        vol_3(:,i_3) = (C_Matrix./dt+G_Matrix)\(Gaussian_F+C_Matrix*vol_start/dt);
        
    else
        vol_3(:,i_3) = (C_Matrix./dt+G_Matrix)\(Gaussian_F+C_Matrix*vol_old/dt);
        
    end
    
    vol_old = vol_3(:, i_3);
        
end

figure(3)
plot(0:step-1, vol_3(7,:), 'r')
hold on
plot(0:step-1, vol_3(1,:), 'b')
title('Vin and Vo with Guassian pulse (mag =1, std dev= 30ms, delay = 60ms')
xlabel('Time in miliseconds')
ylabel('Voltage in volts')
grid on

% Part d) iv. Now that the simulation is complete, the frequency content of
% the input and output signals will be plotted using the built-in matlab
% functions fft() and fftshift().

freq = (-step/2:step/2-1);               

%Plot of Vin, Vo with first input signal in f-domain

fft_vol1_in = fft(vol_1(1, :));
fft_vol1_out = fft(vol_1(7, :));
ffts_vol1_in = fftshift(fft_vol1_in);
ffts_vol1_out = fftshift(fft_vol1_out);

figure(4)
plot(freq, abs(ffts_vol1_in), 'r')
hold on
plot(freq, abs(ffts_vol1_out), 'b')
title('Vin and Vo in f-domain with a step 0-1 at 30ms')
xlabel('frequency in 1/ms')
ylabel('Voltage in volts')
grid on

%Plot of Vin, Vo with second input signal in f-domain

fft_vol2 = fft(vol_2.');
ffts_vol2 = fftshift(fft_vol2);

figure(5)
plot(freq, abs(ffts_vol2(:, 1)), 'r')
hold on
plot(freq, abs(ffts_vol2(:, 7)), 'b')
title('Vin and Vout in f-domain with function sin(2pift), f = 1/30ms')
xlabel('frequency in 1/ms')
ylabel('Voltage in v')
grid on

%Plot of Vin, Vo with third input signal in f-domain

fft_vol3 = fft(vol_3.');
ffts_vol3 = fftshift(fft_vol3);


figure(6)
plot(freq, abs(ffts_vol3(:, 1)), 'r')
hold on
plot(freq, abs(ffts_vol3(:, 7)), 'b')
title('Vin and Vout in f-domain with Guassian pulse (mag =1, std dev = 30ms, delay = 60ms')
xlabel('frequency in 1/ms')
ylabel('Voltage in volts')
grid on    

%%%%%%%%end of part 2%%%%%%%%%%%%%%