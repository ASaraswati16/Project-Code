clear all;
nair=1.0; %refractive index of air
nsubs=1.52; %refractive index of substrate;usually crown glass
lamdac=550; %central wavelength in nm
Ls=450;Le=700; %wavelength range in nm
Np=100; %number of points to be calculated
x=1:Np; % sets up array of 100 equally spaced points
lamda=Ls+x*(Le-Ls)/Np; %sets up array of wavelengths
narc1=1.37; %refractive index of the first anti-reflecting coating
d1=lamdac/(4*narc1); %calculates ideal thickness of layer at central wavelength
theta=2*pi*narc1*d1./lamda; %calculates array of phase or optical thickness
%matrix terms being defined of the characterstic matrix
M11=cos(theta); 
M22=M11;
M12=j*sin(theta)/narc1;
M21=j*narc1*sin(theta);
A=M11+nsubs*M12; 
B=M21+nsubs*M22;
Ramp=(nair*A-B)./(nair*A+B);
Rpow=(abs(Ramp)).^2;
figure
plot(lamda,Rpow*100);
xlabel('wavelength (nm)')
ylabel('reflectance (%)')