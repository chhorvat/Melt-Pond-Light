%drive_fig_2
clear
close all

nx = 100; % Number of lateral grid points
dx = 5; % Grid spacing in meters

%% old way
% load(['../Fig-3/pseudo-pond-stats-' num2str(nx) '-old.mat'])
% 
% 
% 
% [I_ic_save,I_mp_save,I_net_save] = deal(cell(num_it,1)); 
% 
% pool = gcp; 
% 
% if ~pool.Connected
% 
% parpool; 
% 
% end

%%
% parfor i = 1:num_it
%  
%   i
%   
%   [I_ic_save{i},I_mp_save{i},I_net_save{i}] = get_light_field_noSSL(x_save{i},y_save{i},z_save{i},pond_save{i});   
%     
% end
% 
% 
% %
% save(['pseudo-pond-stats-' num2str(nx) 'new.mat'],'*_save','num_it')
% 
% %
% clearvars -except nx

load(['pseudo-pond-stats-' num2str(nx) '-old.mat'])

plot_fig_S3;



