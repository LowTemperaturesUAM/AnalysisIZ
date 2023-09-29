% This program is used to visualize one by one the IZ curves in a .blq, to
% check the data validity. One might add a criteria for selecting certain
% type of IZ curves.

clear all; clc; close all;

% Change the Bias voltge here if it's different than 0.1V.
% offset for 10^6 data: 0.001418; for 10^5 data: 0.00623
% En V
       Bias=0.1;
       offset=0.0097;

%        Bias=0.05;
%        offset=0.0016;

% Change the filename once inside the folder that contains the .blq file
FileNames='Ag_julio23_F_7';
% FileNames='0T_1';

data = ReducedblqreaderV6([FileNames '.blq']);

size=1; %size in Z when making IZ
factor=1/0.92*1000; 

% figure(1) to show the IZ curves z vs G/G0.
fig = figure(1);
fig.WindowKeyPressFcn = @saveEasy;

% xlabel('Z (A)')
% ylabel('log(G/G0)')
%         xlim([-20 2])
%         ylim([-1 5])      

for i = 1000:length(data)
    current(1:2048)= data(i).data(:,2)/(7.748091E-5*Bias)*factor;
    %The vector 'current' equals G/G0. 'fact' is left for correcting the data if the I recording is wrong by a factor.  
    % Initialize x vector (reserved for z position in AA, but by now a vector of indices.)
    x = 1:2048;
    y = real(log10(current+offset));
    % y = log(G/G0)
    % offset for 10^6 data: 0.001418; for 10^5 data: 0.00623, modified before.
    
    if mod(i,2)&(current(1)<0.05) % We only look at approaching curves
         flag=1;
         j = 1;
         
         while ((current(j+2)-current(j))<0.4) %Find position index j for the origin d=0
             j = j+1;
             if (j>1800)|(current(j)>0.8) %curve not valid if the steps occur too close to the edge.
                 flag=0; 
                 break
             end
         end
         
         x = (j-x)/1024*0.1*size*14*1.2361*10*1.3318;
         % Rescaling of x vector for the z position in unit of AA. Calibration 1.2361 nm/V. Full
         % size=1 -- voltage in Z ranges from -100 mV to 100 mV with 2048 points. 
         % 1.3318:correcting factor for the z scaling
         
         m=0;
         logG=[];
         d=[];
         
         
         jj=j+4;
         while (current(jj)>0.65)&(jj<(j+30))
             jj=jj+1;
         end
         if jj<j+25
             flag=1;
         end
         
         
         if flag & (j>2)%If the curve is correct, select the linear part and record the two coordinate in the two vectors d and logG
                        %Do the linear fit --> value stocked in the vector "fit", find the slor "kk" to
                        %calcutae the work function and find the value of
                        %logGa
             jj = j-2;
             while y(jj)>-2.8 % Find the points (d, logG) along the linear segment
                 m=m+1;
                 logG(m)=y(jj);
                 d(m)=(j-jj)/1024*0.1*size*14*1.2361*10*1.3318; % d(m) in Angstrom, d<0 for contact.
                 jj=jj-1;
                 if jj==0
                     break
                 end
             end
             
             if (m>10)&(y(jj+1)<-2.5)&(y(jj+2)<-2.5) %if the segment is long enough and does not start with a jump
                 
                 linear_fit = polyfit(d,logG,1);
                 fit = polyval(linear_fit,d); % The fitting line
                 kk = linear_fit(1) % slope
                 work_function = kk*kk/0.4343/0.4343
                 logGa = linear_fit(2); 
                 Gb = (current(j+5)+current(j+6)+current(j+7)+current(j+8)+current(j+9))/5;
                 
                 var=(logG - fit).^2;
                 var_=mean(var);
                 
                 if var_>0.008
                     flag=0;
                 end
             else                
                 flag=0;
             end
         end
         
         if flag

            i
            var_
            
            figure(1)
                    
            plot(x,current,'b--o','LineWidth',0.0001)
            xlabel('d (AA)')
            ylabel('G/G0')
            xlim([-5 5])
            ylim([0 4])        

            fig.UserData = strcat(FileNames,'-(',num2str(i),').txt'); 
            figure(2)
                    
            plot(x,y,'b--o','LineWidth',0.0001)
            hold on
            plot(d,logG,'go')
            plot(d,fit,'-r')
            xlabel('d (AA)')
            ylabel('log(G/G0)')
            xlim([-5 5])
            ylim([-3 1.5])
            hold off
            
            figure(1)

            uiwait(fig)
         end
         
    end        


end
    
    



 
  

