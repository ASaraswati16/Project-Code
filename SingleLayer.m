clear all;
nair=1.0;
nsubs=1.52;
lamdac=550;
Ls=450;Le=700;
Np=100;
x=1:Np;
lamda=Ls+x*(Le-Ls)/Np;
narc1=1.37;
d1=lamdac/(4*narc1);
theta=2*pi*narc1*d1./lamda;
M11=cos(theta);
M22=M11;
M12=j*sin(theta)/narc1;
M21=j*narc1*sin(theta);
%A=M11*n0+ns*n0*M12; % Equation 3.12
A=M11+nsubs*M12; % Equation 3.12

B=M21+nsubs*M22;
%Ramp=(A-B)./(A+B);
Ramp=(nair*A-B)./(nair*A+B);
Rpow=(abs(Ramp)).^2;
figure
plot(lamda,Rpow*100);
xlabel('wavelength nm')
ylabel('reflectance %')