clc
clear
addpath('Results', 'origin_script') ;

% ----------------------------------------------------------------------------------------------------
MT = 2;   MR = 2; mod = '64QAM'; 
% ---------------------------------------------------------------------------------------------------- 
plot11 = ['Plot1_MMSE_UNCODED', num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot12 = ['Plot1_SD_UNCODED',   num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot13 = ['Plot1_FSD_UNCODED',  num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot14 = ['Plot1_IFSD_UNCODED', num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot15 = ['Plot1_LSFSD_UNCODED',num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];

% -- collect 10 files
number = 3;

% -- load simulation data and consolidate
   [BER_MMSE_64QAM_2x2, BER_SD_64QAM_2x2, ...
    BER_FSD_64QAM_2x2, BER_IFSD_64QAM_2x2, ...
    BER_LSFSD_64QAM_2x2]=plot1_BER_Compile (plot11,plot12, plot13, plot14, plot15, number);

   [COUNT_FSD_64QAM_2x2, COUNT_IFSD_64QAM_2x2, ...
    COUNT_LSFSD_64QAM_2x2]=plot3_COUNT_Compile (plot11,plot12, plot13, plot14, plot15, number);
% save data into mat file and export

FileName = 'Compiled_Data_Origin/plot1_64QAM_2x2.mat';
save(FileName,'BER_MMSE_64QAM_2x2','BER_SD_64QAM_2x2', 'BER_FSD_64QAM_2x2',...
              'BER_IFSD_64QAM_2x2', 'BER_LSFSD_64QAM_2x2');

FileName = 'Compiled_Data_Origin/plot3_64QAM_2x2.mat';
save(FileName,'COUNT_FSD_64QAM_2x2','COUNT_IFSD_64QAM_2x2', 'COUNT_LSFSD_64QAM_2x2');


MT = 4;   MR = 4; mod = '64QAM'; 
% ---------------------------------------------------------------------------------------------------- 
plot11 = ['Plot1_MMSE_UNCODED', num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot12 = ['Plot1_SD_UNCODED',   num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot13 = ['Plot1_FSD_UNCODED',  num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot14 = ['Plot1_IFSD_UNCODED', num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot15 = ['Plot1_LSFSD_UNCODED',num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];

% -- collect 10 files
number = 3;

% -- load simulation data and consolidate
   [BER_MMSE_64QAM_4x4, BER_SD_64QAM_4x4, ...
    BER_FSD_64QAM_4x4, BER_IFSD_64QAM_4x4, ...
    BER_LSFSD_64QAM_4x4]=plot1_BER_Compile (plot11,plot12, plot13, plot14, plot15, number);

   [COUNT_FSD_64QAM_4x4, COUNT_IFSD_64QAM_4x4, ...
    COUNT_LSFSD_64QAM_4x4]=plot3_COUNT_Compile (plot11,plot12, plot13, plot14, plot15, number);
% save data into mat file and export

FileName = 'Compiled_Data_Origin/plot1_64QAM_4x4.mat';
save(FileName,'BER_MMSE_64QAM_4x4','BER_SD_64QAM_4x4', 'BER_FSD_64QAM_4x4',...
              'BER_IFSD_64QAM_4x4', 'BER_LSFSD_64QAM_4x4');

FileName = 'Compiled_Data_Origin/plot3_64QAM_4x4.mat';
save(FileName,'COUNT_FSD_64QAM_4x4','COUNT_IFSD_64QAM_4x4', 'COUNT_LSFSD_64QAM_4x4');



MT = 6;   MR = 6; mod = '64QAM'; 
% ---------------------------------------------------------------------------------------------------- 
plot11 = ['Plot1_MMSE_UNCODED', num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot12 = ['Plot1_SD_UNCODED',   num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot13 = ['Plot1_FSD_UNCODED',  num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot14 = ['Plot1_IFSD_UNCODED', num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];
plot15 = ['Plot1_LSFSD_UNCODED',num2str(MT),'x', num2str(MR), '_',num2str(mod), '_RunID_'];

% -- collect 10 files
number = 3;

% -- load simulation data and consolidate
   [BER_MMSE_64QAM_6x6, BER_SD_64QAM_6x6, ...
    BER_FSD_64QAM_6x6, BER_IFSD_64QAM_6x6, ...
    BER_LSFSD_64QAM_6x6]=plot1_BER_Compile (plot11,plot12, plot13, plot14, plot15, number);

   [COUNT_FSD_64QAM_6x6, COUNT_IFSD_64QAM_6x6, ...
    COUNT_LSFSD_64QAM_6x6]=plot3_COUNT_Compile (plot11,plot12, plot13, plot14, plot15, number);
% save data into mat file and export

FileName = 'Compiled_Data_Origin/plot1_64QAM_6x6.mat';
save(FileName,'BER_MMSE_64QAM_6x6','BER_SD_64QAM_6x6', 'BER_FSD_64QAM_6x6',...
              'BER_IFSD_64QAM_6x6', 'BER_LSFSD_64QAM_6x6');

FileName = 'Compiled_Data_Origin/plot3_64QAM_6x6.mat';
save(FileName,'COUNT_FSD_64QAM_6x6','COUNT_IFSD_64QAM_6x6', 'COUNT_LSFSD_64QAM_6x6');
