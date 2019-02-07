% filename='/Users/msachs2/Documents/Code/outputs/fastMD_output/test/testfile.h5';
filename='/home/xshang/Codes/2018_GLE_DPD/esamc/matlab/testfile.h5';
potential_traj = h5read(filename,'/traj/potential');
momentum_traj = h5read(filename,'/traj/momentum');
position_traj = h5read(filename,'/traj/position');
force_traj = h5read(filename,'/traj/force');
laplace_traj = h5read(filename,'/traj/laplace');
%force_traj = h5read(filename,'/traj/force');
 
kinE = .5*sum(momentum_traj.^2,1);
potE =  sum(potential_traj,1);
 
sdim = 3;% space dimension
N = size(momentum_traj,1)/sdim;
N_t = size(momentum_traj,2);
 
%%
figure;
plot(sum(momentum_traj,1))
hold on
plot(sum(momentum_traj(1:3:end,:),1))
plot(sum(momentum_traj(2:3:end,:),1))
title('Total Momentum')
 
% %%plot(sum(momentum_traj))
% figure;
% plot(sum(position_traj(1:2:end,:),1))
% hold on
% plot(sum(position_traj(2:2:end,:),1))
% title('Center of mass')
% %%
% figure;
% for i = 1:size(position_traj,1)
%     plot(position_traj(i,:))
%     hold on
% end
% title('Position 1,2')
% 
% figure
% plot(potE)
% hold on
% plot(kinE)
% plot(kinE+potE)
% title('Energy')
% legend('Potential','Kinetic', 'Total Energy')
%%
%%
% figure
% cum_avTemp = cumsum(sum(momentum_traj.^2,1)/((N-1)*sdim))./(1:N_t);
% plot(cum_avTemp)
% title('Kinetic temperature')
%%
%%
figure
KT = sum(momentum_traj.^2,1)/((N-1)*sdim);
plot(KT)
title('Kinetic temperature')
%%
%%
% figure
% cum_avConfigTemp = abs(cumsum(sum(force_traj.^2,1))./cumsum(sum(laplace_traj,1)));
% plot(cum_avConfigTemp(100:end))
% %set(gca, 'YScale', 'log')
% title('Configurational temperature')
% %%
% %%
% figure
% cum_avConfigTemp2 = cumsum(sum(-force_traj.*position_traj,1)/(N*sdim))./(1:N_t);
% plot(cum_avConfigTemp2(10:end))
% %set(gca, 'YScale', 'log')
% title('Configurational temperature')
%%
%%
figure
CT = sum(force_traj.^2,1)./sum(laplace_traj,1);
plot(CT)
%set(gca, 'YScale', 'log')
title('Configurational temperature')
%%
%%
