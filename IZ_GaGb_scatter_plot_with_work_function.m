% Main program to treat the quantum conductance data.

% We use this to make log(Ga/G0)-Gb 2D density plot directly out of IZ data
% files.
% It visits all data IZ files in a folder, finds the all log(Ga/G0)-Gb 
% data pairs and records them in 2 two-column vectors, and plot the 2D histogram for the data in the whole folder.
% List of data pairs for D1 and D2 contacts can be collected in Data1 after execution.
% List of data pairs for D3 (multi-atom) contacts can be collected in Data2 after execution.

clear all; clc; close all;

% Put direction of folder, empty if already inside the folder
dir_ = [];
label = '0T_'; % name of the folder

% Read all the .blq files inside the folder
Files=dir(fullfile(dir_,'*.blq'));

% Change the Bias voltge here if it's different.
Bias=0.05; 
% Bias=0.05;% En V, for Au.
factor=1/0.91; % Correction of the current scale, necessary because of the resistance in series. We put it to bring the peak in the usual histogram to 1 G0.
ff=1.3318;% Correction in the length scale, for the uncertainty in piezo calibration. We put it so that the obtained work function for Au is around 5 eV.
size=1; % Size of Z in DAC 6 used during the experiment. Usually with size 1.
offset=0.0097; %Ag 
offset=0.0016; %Au
% Offset in the IV converter. Important for the low conductance part of the data (for the correct value of Ga). 
% For 10^6 data: 0.001418; for 10^5 data: 0.00623
% You can check the offset with the 'IZ_data_checker' program

% Number of points inside the 2D histogram
% Ny=100;
% Nx=100; % Nx=118 for size 1; Nx=? for size 2
% M=zeros(Ny,Nx);

% Border of the log(Ga/G0)-Gb plot.
XXmin=-3;
XXmax=0;
YYmax=2.5;
YYmin=0;

cc=0.55;   %Parameter that will be used to determine the jump position in the IZ curves. See line 80.

%'aaa' is a parameter depending on the noise level in log(G/G0). Usually at 0. 
aaa=0;

%Initialize vector 'Result', for writing the log(Ga/G0)-Gb data pairs
% Data1=[]; % List of data pairs for D1 and D2 contacts can be collected after execution.
% Data2=[]; % List of data pairs for D1 and D2 contacts can be collected after execution.
Result=[]; %Column 1 for log(Ga); column 2 for Gb; Column 3 for work function in eV.
n=0;


for k=1:length(Files)% Visiting each of the files
    
    file_number=k
    FileNames=Files(k).name;
    data = ReducedblqreaderV6([dir_ FileNames]); 

%     if k>5
%         factor=1/0.91*1000;
%     else
%         factor=1/0.91; 
%     end

    % The following lines are reserved for changing parameters. Data sets 1-3 at 20T are with
% %     Bias=0.1 and a larger offset.
%     if k<4
%         Bias=0.1;
%         factor=1/0.91;
%         offset=0.0075;
%     else
%         Bias=0.05;
%         factor=1/0.91;
%         offset=0.00155;
%     end


    for i = 1:length(data) %Reading the IZ curves one by one
        
        current(1:2048)= data(i).data(:,2)/(7.748091E-5*Bias)*factor;
        % Initialize x vector (reserved for z position in A, initially a vector of indices.)
        x = 1:2048;
        % y vector for log(G/G0)
        y = real(log10(current+offset)); 
        % offset for 10^6 data: 0.001418; for 10^5 data: 0.00623
        % check the offset with the 'IZ_data_checker' program

         % Add a flag, to not to take into account the curves where we never reach 1G0.
        if mod(i,2)&(current(1)<0.05) % We are interested in the tunneling tail of the approaching curves, so only approaching curves are visited.
            j = 1; % j is a vector indice, to be determined at the jump position.
            % Jump position: the point where G increases by more than 0.3 within 5 data
            % points.Or where G reaches cc (parameter defined in line 38) .
            % Two ways of defining jump position leads to similar results.
            flag=1;
            while ((current(j+2)-current(j))<0.4) %ALternative: current(j)<cc
                j = j+1;
                if (j>1800)|(current(j)>0.8) % Add a flag, to not to take into account the curves where we never reach 1G0.
                    flag=0;
                    break
                end
            end
            
            %Plot the first 20 curves, to check if the data is correct
            if i<21
                figure(1)
                hold on
                plot(x,y,'LineWidth',0.0001)
                xlabel('index')
                ylabel('log(G/G0)')
            end
            
            m=0;
            logG=[];
            d=[];
            fit=[];
            
            jj=j+4;
            while (current(jj)>0.65)&(jj<(j+30))
                jj=jj+1;
            end
            if jj<j+25
                flag=0;
            end
            
            % We are going to find the linear part of logG vs d, for the
            % tunneling conductance.
            if flag & (j>2)
                jj = j-2;
%                 if current(j-2)>0.45
%                     break
%                 end
                while y(jj)>-2.8 % If the low conductance data is too noisy, take only the points beyond the noise level. Line 41.
                    m=m+1;
                    logG(m)=y(jj);
                    d(m)=(jj-j)/1024*0.1*size*14*1.2361*10*ff; % d(m) in Angstrom
                    jj=jj-1;
                    if jj==0
                        break
                    end
                end

              
              if (m>10)&(y(jj+1)<-2.5)&(y(jj+2)<-2.5)  % We are only to calculate the work function if the linear part is long enough
                  linear_fit = polyfit(d,logG,1);
                  fit = polyval(linear_fit,d); % The fitting line

                  kk = linear_fit(1); % slope
                  work_function = kk*kk/0.4343/0.4343;
                  logGa = linear_fit(2); 
                  Gb = (current(j+5)+current(j+6)+current(j+7)+current(j+8)+current(j+9))/5;

                  var=(logG - fit).^2;
                  var_=mean(var);
                  if var_>0.008
                      flag=0;
                  end
                  
                  % Figure to check the fit is correct. 
%                   if i< 10
%                     figure(2)
% 
%                     hold on
% 
%                     plot(d,logG)
%                     plot(d,fit)
%                     xlabel('d (AA)')
%                     ylabel('log(G/G0)')
%                   end  
              else
                  flag=0;
              end
            end
            
            
            if flag
                  
              if (m>10)&(work_function<6)&(work_function>1)&(logGa>-3)&(logGa<0)
                    n=n+1;
                    Result(n,1) = logGa;
                    Result(n,2) = Gb;
                    Result(n,3) = work_function;
              end
                    
            end        
              
        end
           
    end
    
end


% Draw the (logGa,Gb) scatter plot.

figure(2)
h1 = gscatter(Data1(:,1),Data1(:,2));
hold on
h2 = gscatter(Data2(:,1),Data2(:,2));
haxis = gca;
xlim = haxis.XLim;
ylim = haxis.YLim;
 



