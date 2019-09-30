function [ ] = CombustionDataProcess( )
%CombustionDataProcess Generates surface fit variables given combustion
%data from RPA. 
% WARNING: Known to work for RPA Lite v1.2. Other versions may require
% updated interface

addpath(fullfile('..', 'Supporting Functions'))

psi_to_Pa = 6894.76;

[FileName,~,~] = uigetfile('*','Select combustion data source.');
[~,name,~] = fileparts(FileName);
savefilename = [name '.mat'];

file_ID = fopen(FileName);
N_header = 8;
for ii = 1:N_header
    fgetl(file_ID);
end
data = textscan(file_ID,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');

OF = data{1};
Pc = data{2}*psi_to_Pa;
CombData.OF_range = [min(OF), max(OF)];
CombData.Pc_range = [min(Pc), max(Pc)];
CombData.n_OF = length(unique(OF));
CombData.n_Pc = length(unique(Pc));
Tc = data{6};
M = data{7}/1000;
gamma = data{8};
c_star = data{10};
CombData.OF = reshape(OF,CombData.n_Pc,CombData.n_OF);
CombData.Pc = reshape(Pc,CombData.n_Pc,CombData.n_OF);
CombData.Tc = reshape(Tc,CombData.n_Pc,CombData.n_OF);
CombData.M = reshape(M,CombData.n_Pc,CombData.n_OF);
CombData.gamma = reshape(gamma,CombData.n_Pc,CombData.n_OF);
CombData.c_star = reshape(c_star,CombData.n_Pc,CombData.n_OF);
save(savefilename,'CombData')

end

