function SaveTRIMFile(Filename, Data, Description)
% Loads positions from SRIM input/output files
% Positions are read in as Angstroms and converted to um
%
%
%   Shaun Kerr, 2017
%   smkerr@gmail.com
%
%%
HeaderLines{1} = '=========== TRIM with various Incident Ion Energies/Angles and Depths =========';
HeaderLines{2} = '= This file tabulates the kinetics of incident ions or atoms.                 =';
HeaderLines{3} = '= Col.#1: Ion Number, Col.#2: Z of atom leaving, Col.#3: Atom energy (eV).    =';
HeaderLines{4} = '= Col.#4-6: Last location:  Col.#4: X= Depth into target.                     =';
HeaderLines{5} = '= Col.#5,6: Y,Z= Transverse axes                                              =';
HeaderLines{6} = '= Col.#7-9: Cosines of final trajectory.                                      =';
HeaderLines{7} = '================ Typical Data File is shown below  ============================';


TableHeader{1} = 'Event  Atom  Energy  Depth   Lateral-Position   ----- Atom Direction ----';
TableHeader{2} = 'Name   Numb   (eV)    _X_(A)   _Y_(A)  _Z_(A)   Cos(X)   Cos(Y)   Cos(Z)';
    
fileID = fopen(Filename,'w');
 for i = 1:length(HeaderLines)
     fprintf(fileID, '%s\r\n',HeaderLines{i});
 end
fprintf(fileID, '%s\r\n',Description)
 for i = 1:length(TableHeader)
     fprintf(fileID, '%s\r\n',TableHeader{i});
 end
 
 for i = 1:length(Data)
     % 6 spaces
     SpaceString = repmat(' ', 6-numel(num2str(i)),1);
     fprintf(fileID, '%s%.0f  %.0f %.6e %.6e   %.3e  %.3e  %.7f %.7f %.7f\r\n',SpaceString,...
         i,Data(i,1),Data(i,2),Data(i,3),Data(i,4),Data(i,5),Data(i,6),Data(i,7),Data(i,8));
 end
%1  1 1.001000E+06 0.000000E+00   0.000E+00  0.000E+00  1.0000000 0.0000000 0.0000000


 fclose(fileID);
