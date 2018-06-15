filename='/Users/msachs2/Documents/Code/outputs/fastMD_output/test/testfile.h5';
potential_traj = h5read(filename,'/traj/potential');
momentum_traj = h5read(filename,'/traj/momentum');
position_traj = h5read(filename,'/traj/position');

kinE = .5*sum(momentum_traj.^2,1);
potE =  sum(potential_traj,1);


%%
figure;
plot(sum(momentum_traj,1))
hold on
plot(sum(momentum_traj(1:3:end,:),1))
plot(sum(momentum_traj(2:3:end,:),1))
title('Total Momentum')

%%plot(sum(momentum_traj))
figure;
plot(sum(position_traj(1:2:end,:),1))
hold on
plot(sum(position_traj(2:2:end,:),1))
title('Center of mass')
%%
figure;
for i = 1:size(position_traj,1)
    plot(position_traj(i,:))
    hold on
end
title('Position 1,2')

figure
plot(potE)
hold on
plot(kinE)
plot(kinE+potE)
title('Energy')
legend('Potential','Kinetic', 'Total Energy')
