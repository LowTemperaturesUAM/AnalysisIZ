% This function helps to make (log(Ga), Gb) density plot starting from list of data pairs.

% Before running this function, create a vector in the Command Window with
% the following line: Data=[];
% Then copy the two column data ((log(Ga), Gb), from origin) in 'Data'

close all; %Closing existing figures

%%Border of the histogram
Ny=100;
Nx=178; % Nx=118 for size 1; Nx=? for size 2
M=zeros(Ny,Nx);
ProfileGa=zeros(Nx);
ProfileGb=zeros(Ny);

% log(Ga/G0) range [-3, 0]
XXmax=0;
XXmin=-3;
% Gb/G0 range [0, 2.5]
YYmax=0;
YYmin=2.5; 

dXX=(XXmax-XXmin)/Nx;
dYY=(YYmax-YYmin)/Ny;
XX=linspace(XXmin,XXmax,Nx);
YY=linspace(YYmin,YYmax,Ny);

Nz=100;
ZZmin=0;
ZZmax=15;
dZZ=(ZZmax-ZZmin)/Nz;
ZZ=linspace(ZZmin,ZZmax,Nz);
WorkFunction=zeros(Nz);

for i= 1:length(Result(:,1))
    a = floor((Result(i,1)-XXmin)/dXX);
    b = floor((Result(i,2)-YYmin)/dYY);
    if (a>0) & (a<Nx) & (b>0) & (b<Ny)
        M(b,a)=M(b,a)+1;
        ProfileGa(a)=ProfileGa(a)+1;
        ProfileGb(b)=ProfileGb(b)+1;
    end
end 

for i= 1:length(Result(:,1))
    c = floor((Result(i,3)-ZZmin)/dZZ);
    if (c>0) & (c<Nz) 
        WorkFunction(c)=WorkFunction(c)+1;
    end
end 

figure(1)
s2=pcolor(XX,YY,M) ;
s2.EdgeColor = 'none';
axis([-3 0 0 2.5])

% The color scale is normalized to the maximum value in the counts.

cc = max(M,[],'all')
caxis([0 floor(cc)]);
colorbar
xlabel('log(Ga/G0)')
ylabel('Gb/G0')

figure(2)
plot(XX, ProfileGa);
xlabel('log(Ga/G0)')
ylabel('Counts')

figure(3)
plot(YY, ProfileGb);
xlabel('Gb/G0')
ylabel('Counts')

figure(4)
plot(ZZ, WorkFunction);
xlabel('Work function (eV)')
ylabel('Counts')
xlim([0 6])

