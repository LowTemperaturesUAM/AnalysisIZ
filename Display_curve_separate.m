% This program is used to visualize one by one the IZ curves in a .blq, to
% check the data validity. One might add a criteria for selecting certain
% type of IZ curves.

clear all; clc; close all;

% Change the Bias voltge here if it's different than 0.1V.
% offset for 10^6 data: 0.001418; for 10^5 data: 0.00623
% En V
       Bias=0.05;
       offset=0.0097;

% Change the filename once inside the folder that contains the .blq file
% FileNames='Ag_julio23_F_8';
FileNames='0T_1';

data = ReducedblqreaderV6([FileNames '.blq']);

size=1; %size in Z when making IZ
factor=1/0.92; 

% figure(1) to show the IZ curves z vs G/G0.
fig = figure(1);
fig.WindowKeyPressFcn = @saveEasy;

% xlabel('Z (A)')
% ylabel('log(G/G0)')
%         xlim([-20 2])
%         ylim([-1 5])      

for i = 1:length(data)
    current(1:2048)= data(i).data(:,2)/(7.748091E-5*Bias)*factor;
    %The vector 'current' equals G/G0. 'fact' is left for correcting the data if the I recording is wrong by a factor.  
    % Initialize x vector (reserved for z position in AA, but by now a vector of indices.)
    x = 1:2048;
    y = real(log10(current+offset));
    % y = log(G/G0)
    % offset for 10^6 data: 0.001418; for 10^5 data: 0.00623, modified before.
    
    if (mod(i,2)==0)&(current(2048)<0.05) % We only look at approaching curves
         flag=1;
         j = 2048;
         
         while ((current(j-2)-current(j))<0.5) %Find position index j for the origin d=0
             j = j-1;
             if (j<400)|(current(j)>0.8) %curve not valid if the steps occur too close to the edge.
                 flag=0; 
                 break
             end
         end
         
         x = (x-j)/1024*0.1*size*14*1.2361*10*1.3318;
         % Rescaling of x vector for the z position in unit of AA. Calibration 1.2361 nm/V. Full
         % size=1 -- voltage in Z ranges from -100 mV to 100 mV with 2048 points. 
         % 1.3318:correcting factor for the z scaling
         
         m=0;
         logG=[];
         d=[];
                  
         jj=j-4;
         while (current(jj)>0.65)&(jj>(j-30))
             jj=jj-1;
         end
         if jj>j-25
             flag=0;
         end
         
         if flag

            i
            
            figure(1)
                    
            plot(x,current,'b--o','LineWidth',0.0001)
            xlabel('d (AA)')
            ylabel('G/G0')
            xlim([-9 5])
            ylim([0 4])        

            fig.UserData = strcat(FileNames,'-(',num2str(i),').txt'); 

            uiwait(fig)
         end
         
    end        


end
    
    



 
  

