clear all;
nair=1.0;  %refractive index of air
narc1=2.2; % refractive index of anti-reflection coating 1
narc2=1.37; %refractive index of anti-reflection coating 2
nsubs=1.52; %refractive index of substrate i.e. crown glass
lamdac=500; %sets centre avelength in nm
Ls=350;Le=750; %range in nm
Np=1000;%array points
x=1:Np;
lamda=Ls+x*(Le-Ls)/Np;
d=lamdac/4;
%adding the layer
N=6;
%making a for loop
Rpow = zeros(size(lamda));
for i = 1:length(lamda)
% Calculate phase shifts for each layer at current wavelength
theta = 2 *pi*d./lamda(i); % High index layer
    
 % Transfer matrix for a single period (one double layer of n1 and n2)
    Matrix1 = [cos(theta), (1i*sin(theta))/narc1 ; 1i*narc1*sin(theta), cos(theta)];
    Matrix2 = [cos(theta), (1i*sin(theta))/narc2 ; 1i*narc2*sin(theta), cos(theta)];
    % Total characteristic matrix for N double layers
    M = (Matrix1 * Matrix2)^N;
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

