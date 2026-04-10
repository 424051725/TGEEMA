%% 2D simulation
clear
clc
%% setting parameters
n=400;       % sample size
m=4;         % sample size
max_rank=5;  % max CP rank, candidate models are rank 1 to rank max_rank
fold = 5;    % cross-validation fold
workcorr='equicorr'; %  real corelation matrix, options: {'equicorr', 'AR1', 'unstructured'}
link='normal';

%% generate data and run main code
for pic_name = 1:1:1  % pic_name can be choosen from 1--7, see figures in pic_64 folder
    for repeat = 1:1:1
        run('data_2d.m');
        run('simulation_func.m');
    end
end


