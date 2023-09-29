clear all; clc; close all;

%Carpeta Histograma fuera de las otras---------------------------%
%num_folder = '2'
dir_ = ['\'];
label = ['0T_100mV_'];

%Carpeta Histograma en la total---------------------------
% dir_ = '';
% label = '10T_100mV';
%---------------------------------------------------------


Files=dir(fullfile('*.blq'));

Bias=0.1; % En V
count = 0;
total_num = 0;

%  d=0.0007877;  %Good value for Bias=0.2 V
%  d=0.003936;   %Good value for Bias=0.01 V
d=0.00788;  %Good value for Bias=0.1 V
edge=[0:d:4.5];
M=floor(4.5/d);

factor=1/0.91*1000;

for k=1:length(Files)
   file_number=k
   FileNames=Files(k).name;
   data = ReducedblqreaderV6(FileNames);

   for i = 1:M
      Nseparate(i,k)=0;
      Napproach(i,k)=0;
      Nall(i,k)=0;
   end
   
   for i = 1:length(data)
       if mod(i,2)        
         current(1:2048)=(data(i).data(:,2) - data(i).data(1,2))/(7.748091E-5*Bias)*factor;
         N=histcounts(current,edge);        
         N(1)=0;
         
%          if sum(N(3000:3850))>1
         if sum(N(floor(M*2.4/4.5):floor(M*3/4.5)))>1
                  Napproach(:,k)=Napproach(:,k)+N.';
                  Nall(:,k)=Nall(:,k)+N.';
                  count=count+1;
         end
      
       else
         current(1:2048)=(data(i).data(:,2) - data(i).data(2048,2))/(7.748091E-5*Bias)*factor;
         N=histcounts(current,edge);        
         N(1)=0;
         
%          if sum(N(3000:3850))>1
           if sum(N(floor(M*2.4/4.5):floor(M*3/4.5)))>1
                  Nseparate(:,k)=Nseparate(:,k)+N.';
                  Nall(:,k)=Nall(:,k)+N.';
                  count=count+1;
           end
        
       end
   end
   
   total_num = total_num +i;
   amount_of_curves=i
   
end

total_num

conductance=edge(1:M).';
% total_separate = sum(Nseparate,2);
% total_approach = sum(Napproach,2);
% total = sum(Nall,2);
Table = table(conductance,Nseparate,Napproach,Nall);

figure(1)
plot(conductance,Nall)
hold on
plot(conductance,Nseparate)
plot(conductance,Napproach)
xlabel('d (AA)')
ylabel('log(G/G0)')
hold off


% writetable(T_separate,strcat('histogram\histogram_',label,'_separate.txt'),'Delimiter','\t');
% writetable(T_approach,strcat('histogram\histogram_',label,'_approach.txt'),'Delimiter','\t');
% writetable(T_total,strcat('histogram\histogram_',label,'.txt'),'Delimiter','\t');

count