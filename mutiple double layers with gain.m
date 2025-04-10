clear all;
no=1.0;  %refractive index of air
n1=2.1; %refractive index of anti-reflection coating 1
n2=1.37; %refractive index of anti-reflection coating 2
ns=1.52; %refractive index of substrate i.e. crown glass
lamda0=510; %sets centre avelength in nm 
Ls=350;Le=800; %range in nm
Np=1000;%array points
x=1:Np;
lamda=Ls+x*(Le-Ls)/Np;
d=lamda0/4;
a = input('enter the value of a'); %value would be -0.05 or 0.001
nm = input('enter the value of nm');%value would be 1.37
alpha = (a*nm*4./lamda0); % a= (alpha*lamda0)/(4*nm) nm in this case is equal to n2
%adding the layer
N=6;
%making a for loop
Rpow = zeros(size(lamda));
Tpow = zeros(size(lamda));
for i = 1:length(lamda)
% Calculate phase shifts for each layer at current wavelength
    theta = (2 *pi*d./lamda(i)) - (1i*a); % defining theta based on equation 6.4 pg 34
    
 % Transfer matrix for a single period (one double layer of n1 and n2)
    Matrix1 = [cos(theta), (1i*sin(theta))/n1 ; 1i*n1*sin(theta), cos(theta)];
    Matrix2 = [cos(theta), (1i*sin(theta))/n2 ; 1i*n2*sin(theta), cos(theta)];
    % Total matrix for N double layers
    M = (Matrix1 * Matrix2)^N;
    M11 = M(1,1);
    M12 = M(1,2);
    M21 = M(2,1);
    M22 = M(2,2);
    % Reflectivity calculation using M_total
    Ramp = ((M11*no) + (ns*no*M12) - M21-(ns*M22))/((M11*no) + (ns*no*M12) + M21 + (ns*M22));
    Rpow(i)=(abs(Ramp)).^2;
    % Plot reflectance vs. wavelength
    A = M11*no + M12*ns;
    B = M21*no + M22*ns;
    Tpow(i) = (4*ns*no) / ((abs((no*A)+B)).^2);
end
figure;
subplot(1,2,1);
plot(lamda,(Rpow*100));
xlabel('Wavelength (nm)');
ylabel('Reflectance (%)');
subplot(1,2,2);
plot(lamda,(Tpow*100));
xlabel('wavelength(nm)');
ylabel('Transmittance (%)');

