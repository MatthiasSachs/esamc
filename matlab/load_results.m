%  filename='/Users/msachs2/Documents/Code/outputs/fastMD_output/test/testfile.h5';
filename='/home/xshang/Codes/2019_GLE_DPD/esamc/matlab/testfile.h5';
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
Nd = (N-1)*sdim; % Number of degrees of freedom in DPD
 
%%
% figure;
% plot(sum(momentum_traj,1))
% hold on
% plot(sum(momentum_traj(1:3:end,:),1))
% plot(sum(momentum_traj(2:3:end,:),1))
% title('Total Momentum')
  
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
% figure
% KT = sum(momentum_traj.^2,1)/(Nd);
% plot(KT)
% title('Average Kinetic Temperature')
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
% figure
% CT = sum(force_traj.^2,1)./sum(laplace_traj,1);
% plot(CT)
% %set(gca, 'YScale', 'log')
% title('Average Configurational Temperature')
%%
%%
% figure
% PE = sum(potential_traj.^2,1)/(Nd);
% plot(PE)
% %set(gca, 'YScale', 'log')
% title('Average Potential Energy')
%%
%%
pow = 0;
dt = 0.05*1.15^(pow);
KT = sum(momentum_traj.^2,1)/(Nd);
CT = sum(force_traj.^2,1)./sum(laplace_traj,1);
len = length(KT);
dt*len
xx(1) = dt;
xx(2) = mean(KT(round(0.2*len):end));
xx(3) = mean(CT(round(0.2*len):end));
x10(pow+1,:) = xx;
%%
%%

% x=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10;

