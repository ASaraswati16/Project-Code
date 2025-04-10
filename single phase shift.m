clear all;
nair=1.0;  %refractive index of air
narc1=2.0; %refractive index of anti-reflection coating 1
narc2=1.37; %refractive index of anti-reflection coating 2
nsubs=1.52; %refractive index of substrate i.e. crown glass
lamdac=535; %sets centre avelength in nm
Ls=350;Le=750; %range in nm
Np=1000;%array points
x=1:Np;
lamda=Ls+x*(Le-Ls)/Np;
lambda1 = 450;  % Design wavelength for left stack (nm)
lambda2 = 620;  % Design wavelength for right stack (nm)
N = 8;    % Number of layers for both right and left matrices
nm = input('enter the value of refractive index of middle layer'); %will either take value of narc1 or narc2
theta_i = input ('enter the value of the phase shift of the middle layer'); %will either take value of pi or pi/2
M_i = [cos(theta_i) , 1i*sin(theta_i)/nm ; 1i*nm*sin(theta_i) , cos(theta_i)]; %matrix for the inserted layer
d=lamdac/4;
%making a for loop
Rpow = zeros(size(lamda));
for i = 1:length(lamda)
    % Calculate phase shifts for each layer at current wavelength
theta = 2 *pi*d./lamda(i); % High index layer
% Transfer matrix for a single layer on left stack
    MatrixL1 = [cos(theta), (1i*sin(theta))/narc2 ; 1i*narc2*sin(theta), cos(theta)];
    MatrixL2 = [cos(theta), (1i*sin(theta))/narc1 ; 1i*narc1*sin(theta), cos(theta)];
    % Total matrix for N double layers for left stack
    M_L = (MatrixL1 * MatrixL2)^N;
    % Transfer matrix for a single layer on right stack
    MatrixR1 = [cos(theta), (1i*sin(theta))/narc1 ; 1i*narc1*sin(theta), cos(theta)];
    MatrixR2 = [cos(theta), (1i*sin(theta))/narc2 ; 1i*narc2*sin(theta), cos(theta)];
    % Total matrix for N double layers for right stack
    M_R = (MatrixR1 * MatrixR2)^N;
% Total characteristic matrix
    M = (M_L * M_i * M_R);
    M11 = M(1,1);
    M12 = M(1,2);
    M21 = M(2,1);
    M22 = M(2,2);
Ramp = ((M11*nair) + (nsubs*nair*M12) - M21-(nsubs*M22))/((M11*nair) + (nsubs*nair*M12) + M21 + (nsubs*M22));
Rpow(i)=(abs(Ramp)).^2;
end
figure;
plot(lamda,(Rpow*100));
xlabel('Wavelength (nm)');
ylabel('Reflectance (%)');

