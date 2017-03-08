function [E, pos, vdir, inum] = LoadSRIMFile(Filename, type)
% Loads positions from SRIM input/output files
% Positions are read in as Angstroms and converted to um
%   type = [output or input]
%    
%   E (MeV)
%   pos: 
%           1st column = Depth (A)
%           2nd column = y position (A)
%           3rd column = z position (A)
%   inum - ion number, used to compare to input ions
%
%
%   Shaun Kerr, 2017
%   smkerr@gmail.com
%%
switch type
    case 'input'
        headerlines = 10;
        formatSpec = '%2s%7s%15s%14s%18s%18s%14s%14s%s%[^\n\r]';
    case 'output'
        headerlines = 12;
        formatSpec = '%1s%5f%3f%13f%14f%11f%11f%11f%10f%f%[^\n\r]';
end



file = importdata(Filename, ' ', headerlines);

LongFile = 0;
if length(file.data)>=10000
    LongFile = 1; 
end   

%[x,y,z] = importSRIMfile(Filename);

%anum = file{:,3};  % Atomic Number

scale1=10000; % convert  A to um
scale2= 10000; % output
switch type
    case 'input'
        inum = file.data(:,1);  % Ion Number
        E = file.data(:,3)/10^6; % Energy (eV--> MeV)
        pos(:,1) = file.data(:,4)/scale1; % Depth (Angstroms)
        pos(:,2) = file.data(:,5)/scale1; % Lateral Pos. y (Angstroms) 
        pos(:,3) = file.data(:,6)/scale1; % Lateral Pos. z (Angstroms) 
        vdir(:,1) = file.data(:,7);  % Direction, Cos(x)
        vdir(:,2) = file.data(:,8);  % Direction, Cos(y)
        vdir(:,3) = file.data(:,9);  % Direction, Cos(z)
    case 'output'
        scan2 = 0; 
        % 1st and 2nd columns are combined in text file, so 8 columns
        if size(file.data,2) == 8 
            inum = file.textdata(13:end,1);
            inum = strrep(inum,'T', ''); % splice out 'T'
            inum = str2double(inum); 
            offset = 1;
        elseif isnan(file.data(end,9)) % starts with 9 columns, goes to 8
            inum = file.textdata(13:end,1);
            inum = strrep(inum,'T', ''); % splice out 'T'
            inum = str2double(inum); 
            offset = 1;
            scan2 = 1; 
        else
          inum = file.data(:,1);   
          offset = 0;
        end

        E = file.data(:,3-offset)/10^6; % Energy (eV--> MeV)
        pos(:,1) = file.data(:,4-offset)/scale2; % Depth (Angstroms --> um)
        pos(:,2) = file.data(:,5-offset)/scale2; % Lateral Pos. y (Angstroms--> um) 
        pos(:,3) = file.data(:,6-offset)/scale2; % Lateral Pos. z (Angstroms--> um) 
        vdir(:,1) = file.data(:,7-offset);  % Direction, Cos(x)
        vdir(:,2) = file.data(:,8-offset);  % Direction, Cos(y)
        vdir(:,3) = file.data(:,9-offset);  % Direction, Cos(z)

        % # of columns changes after 9999 particles
        if scan2
            Etemp = file.data(:,3)/10^6; % Energy (eV--> MeV)
            xtemp = file.data(:,4)/scale2; % Depth (Angstroms --> um)
            ytemp = file.data(:,5)/scale2; % Lateral Pos. y (Angstroms--> um) 
            ztemp = file.data(:,6)/scale2; % Lateral Pos. z (Angstroms--> um) 
            xDirtemp = file.data(:,7);  % Direction, Cos(x)
            yDirtemp = file.data(:,8);  % Direction, Cos(y)
            zDirtemp = file.data(:,9);  % Direction, Cos(z)

            [a b] = max(file.data(:,1));
            s = b; 
            E(1:s) = Etemp(1:s);
            pos(1:s,1) = xtemp(1:s);
            pos(1:s,2) = ytemp(1:s);
            pos(1:s,3) = ztemp(1:s);
            vdir(1:s,1) = xDirtemp(1:s);
            vdir(1:s,2) = yDirtemp(1:s);
            vdir(1:s,3) = zDirtemp(1:s);

            inum(1:s) = file.data(1:s,1);   
        end
end
    %dirx = [dirx,file.data(:,7)]; 
    %diry = [diry,file.data(:,8)];  % Direction, Cos(y)
    %dirz = [dirz,file.data(:,9)];  % Direction, Cos(z)      
    %catch
    %   disp('Error in reading output.TXT')
    %    return
    %end
        
    %
    
    display(['Loaded - ' Filename])
end