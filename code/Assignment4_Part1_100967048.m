%% Assignment 4 Part 1

% Richard Finney 100967048

%In this part of the assignment, we are repeating the work completed in
%pa9 and reporting on it.

%In part a)the C, G matrices and the F vector was created to describe the
%circuit network

%In pat b), the input voltage was DC swept from -10V to 10V and Vo and V3 was
%plotted. Then for the AC case, Vo was plotted as a function of w, and the
%gain, Vo/V1 was plotted in dB. Then,for the AC case, the gain was plotted
%as a function of random peturbations on C using a normal   distribution
%with stf = 0.5 at w = pi. The gain was then plotted using histograms

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

% Definition of Matrices
C_matrix = [0 0 0 0 0 0 0;
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

% Defining DC and AC voltage matrices, as well as the F matrix

V_DC = zeros(7,1);       
V_AC = zeros(7,1);       
F_Matrix = zeros(7,1);

% DC Sweep Plot
for vol = -10:0.1: 10
    
    F_Matrix(1,1) = vol;
    V_DC = G_Matrix\F_Matrix;                          % DC sweep calculation
    
    figure(1)
    plot(vol, V_DC(7,1), 'g.')
    hold on
    
    plot(vol, V_DC(4,1), 'r.')
    hold on
    title('DC sweep of Vo (green) and V3 (red)')
    xlabel('Vin')
    ylabel('V')
    
end

% AC Sweep and Gain Plot
w = logspace(1,2,500);                  
F_Matrix(1) = 1;

for i = 1:length(w)
    
    V_AC = (G_Matrix+C_matrix*1j*w(i))\F_Matrix;  % calculating the voltage matrix using AC sweep
    figure(2)
    semilogx(w(i), abs(V_AC(7,1)), 'b.')
    hold on
    title('AC sweep of Vo')
    
    dB = 20*log(abs(V_AC(7,1))/F_Matrix(1));    % Calculating the gain
    figure(3)
    plot(i, dB, 'r.')
    hold on
    title('Gain Vo/Vin in dB')
end

% AC case: voltage gain calculation as a function of random perturbations
% on C using a normal distribution with std = .05 and w = pi
pert =  0.25 + 0.05.*randn(1,1000);
w = pi;
Gain = zeros(1000,1);

for n = 1:length(Gain)
    
    C = pert(n);
    C_matrix(2,1) = -C;
    C_matrix(2,2) = C;
    V_AC = (G_Matrix+C_matrix*1j*w)\F_Matrix;   % Voltage calculation using AC sweep
    Gain(n,1) = abs(V_AC(7,1))/F_Matrix(1);    % Gain calculation using AC voltage
end

% Gain histogram
figure(4)
hist(Gain,50);
title('Histogram of Gain with random perturbations')