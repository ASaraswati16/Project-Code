clear all;
no=1;  %refractive index of air
n1=3.7; %GaAs
%no=n1; %put '%' when using the other case
n2=3.1; %AlGaAs
ns=1; %refractive index of substrate taken the same as air
%ns=n1; %put '%' when using the other case
lamda0=100; %sets centre avelength in micrometer(3 THz in microm is 100)
Ls=90;Le=115; %range in micrometer (2.6 THz is 115 mircom to 3.2 THz is 94 microm)
Np=1000;%array points
x=1:Np;
lamda=Ls+x*(Le-Ls)/Np;
N = 13;    % Number of layers for both exteme right and exteme left
nm = n1; 
theta_m1 = 0.5*pi; %this is when Lamda is 2
M_m1 = [cos(theta_m1) , 1i*sin(theta_m1)/nm ; 1i*nm*sin(theta_m1) , cos(theta_m1)]; %defect matrix when Lamda is 2
theta_m2 = 2*pi; %this is when Lamda is 3.5
M_m2 = [cos(theta_m2) , 1i*sin(theta_m2)/nm ; 1i*nm*sin(theta_m2) , cos(theta_m2)]; %defect matrix when Lamda is 3.5
theta_m3 = 1*pi; %this is for when Lamda is 2.5
M_m3 = [cos(theta_m3) , 1i*sin(theta_m3)/nm ; 1i*nm*sin(theta_m3) , cos(theta_m3)]; %defect matrix when Lamda is 2.5
d=lamda0/4;
a= -0.0035; % a = -0.008 when n0=ns=n1 ; a= -0.0035 when n0=ns=1
%a=0; %put '%' when using the other case
%making a for loop
Rpow = zeros(size(lamda));
for i = 1:length(lamda)
    % Calculate phase shifts for each layer at current wavelength
theta = (2 *pi*d./lamda(i)) - (1i*a); % High index layer
% Transfer matrix for a single layer at beginning
    MatrixL1 = [cos(theta), (1i*sin(theta))/n1 ; 1i*n1*sin(theta), cos(theta)];
    MatrixL2 = [cos(theta), (1i*sin(theta))/n2 ; 1i*n2*sin(theta), cos(theta)];
    % Total matrix for 13 double layers for beginning
    M_L1 = (MatrixL1 * MatrixL2)^N;
    % Transfer matrix for a single layer for after the 13 double layers
    MatrixR1 = [cos(theta), (1i*sin(theta))/n1 ; 1i*n1*sin(theta), cos(theta)];
    MatrixR2 = [cos(theta), (1i*sin(theta))/n2 ; 1i*n2*sin(theta), cos(theta)];
    % Total matrix for N double layers for different cases
    M_R1 = (MatrixR1 * MatrixR2);
    M_R2 = (MatrixR1 * MatrixR2)^3;
    %matrix for facet
    M_f = [cos(theta), (1i*sin(theta))/no ; 1i*no*sin(theta), cos(theta)]^50;
% Total characteristic matrix
    M = (M_f * M_L1 * M_m1 * M_R1 * M_m1 * MatrixR2 * M_m2 * MatrixR2 * M_m1 * M_R1 * M_m1 * M_R2 * M_m3 * MatrixR2 * M_m1 * MatrixR2 * M_m1 * M_L1 * M_f);
    M11 = M(1,1);
    M12 = M(1,2);
    M21 = M(2,1);
    M22 = M(2,2);
    %reflectance
    Ramp = ((M11*no) + (ns*no*M12) - M21-(ns*M22))/((M11*no) + (ns*no*M12) + M21 + (ns*M22));
    Rpow(i)=(abs(Ramp)).^2;
    % transmittance
    A = M11*no + M12*ns;
    B = M21*no + M22*ns;
    Tpow(i) = (4*ns*no) / ((abs((no*A)+B)).^2);
end
figure;
subplot(1,2,1);
plot(lamda,(Rpow));
xlabel('Wavelength (um)');
ylabel('Reflectance');
subplot(1,2,2);
plot(lamda,(Tpow*100));
xlabel('wavelength(um)');
ylabel('Transmittance');

