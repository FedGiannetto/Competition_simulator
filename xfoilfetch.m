clear
close all

% Modified from JoshTheEngineer: https://www.youtube.com/watch?v=bWjo3N9COz4


%% INPUT ARGUMENTS

airfoil_dat = 'DAE-21.dat';
airfoilName = airfoil_dat(1:end-4);
AOAstart    = '-1'; % deg
AOAend      = '8';
incr        = '0.25';
ReStart     = 100000;
incrRe      = 100000;
ReEnd       = 900000;

Re = ReStart:incrRe:ReEnd;
Re = num2str(Re(:));

% Data output
saveFlnm = 'Save_Airfoil.txt';

if ~exist('polars', 'dir')
    mkdir('polars')
end

% values of alpha (deg)
alfas = (str2double(AOAend) - str2double(AOAstart))/ str2double(incr);

[num_ansys,n] = size(Re); % Number of analysis (Reynolds)

airfoil_data = cell(num_ansys,6);
%alfa CL CD CM CL_alfa Re

%% ANALISI
% The input of maual commands in XFOIL is simulated

for i=1:1:num_ansys
    % Removes previous analysis
    if (exist(saveFlnm,'file'))
        delete(saveFlnm);
    end
    % Create the airfoil
    fid = fopen('xfoil_input.txt','w');
    fprintf(fid,['LOAD ' airfoil_dat '\n']);
    % Panel check
    fprintf(fid,'PANE\n');
    fprintf(fid,'\n\n');
    % Enters viscous analysis and sets Re
    fprintf(fid,'OPER\n');
    fprintf(fid,'VISC\n');
    fprintf(fid,[Re(i,:) '\n']);
    % Activates file save [saveFlnm]
    fprintf(fid,'Pacc\n');
    fprintf(fid,[saveFlnm '\n']);
    fprintf(fid,'\n');
    % Angles sequence
    fprintf(fid,['Aseq ' AOAstart ' ' AOAend '\n']);
    fprintf(fid,[incr '\n']);
    
    % Close file
    fclose(fid);
    % Run XFoil using input file
    cmd = 'xfoil.exe < xfoil_input.txt';
    [status,result] = system(cmd);
    
%% POLAR EXTRACTION
    fid = fopen(saveFlnm);
    dataBuffer = textscan(fid,'%f %f %f %f %f %f %f','HeaderLines',12,...
                                'CollectOutput',1,...
                                'Delimiter','');
    fclose(fid);
    
    % Struct creation with airfoil data

    alfa = dataBuffer{1,1}(:,1);
    CL = dataBuffer{1,1}(:,2);
    CD = dataBuffer{1,1}(:,3);
    CM = dataBuffer{1,1}(:,5);
    CL_alfa = polyfit(alfa,CL,1);
    airfoil_data(i,:) = {alfa CL CD CM CL_alfa(1) str2double(Re(i,:))};
    
end

figure(1)
hold on
grid on
xlabel('alfa');
ylabel('CL');
for i=1:1:num_ansys
    plot(airfoil_data{i,1},airfoil_data{i,2},'o-')
end
figure(2)
hold on
xlabel('CD');
ylabel('CL');
grid on
for i=1:1:num_ansys
    plot(airfoil_data{i,3},airfoil_data{i,2},'o-')
end

% Simulation data to file
save([pwd '\polars\' airfoilName],'airfoil_data')
