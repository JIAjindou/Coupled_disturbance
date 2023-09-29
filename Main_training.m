clear all
clc
close all

%% Initial parameters 
size_case = 10000;
n = 1;  % The dimensionality of system state
m = 1;  % The dimensionality of uncertainty

state = zeros(size_case, n);
uncertainty = zeros(size_case, n);
time_stamp = zeros(size_case, 1);


%% load data
for i = 1:1:100
    for j = 1:1:100
        time_now = (i-1)*0.01;
        time_stamp((i-1)*100+j, 1) = time_now;
        
        state_now = 10 - 0.2 * j;
        state((i-1)*100+j, :) = state_now;

        uncertainty((i-1)*100+j, :) = [-state_now^2 + 50, -10, 0.5] * [1, time_now, time_now^2]';
    end
end

%% 0ptimization problem to solve Theta
% Params need to adjust
delta = 0.01;
p = 2;

s1 = (p+1)^(m+n);
s2 = (p+1)^(m);
Theta = ones(n,s1); 

% Calcualte the closed form solution of the LS problem
  sum1 = zeros(s1, s1);
  sum2 = zeros(n, s1);
for k = 1:1:size_case
    xi = xi_fun(time_stamp(k, 1), p);
    sum1 = sum1 + B_X_fun(state(k, :), p )* xi * xi' * B_X_fun(state(k, :), p )';
    sum2 = sum2 + uncertainty(k, :)' * xi' * B_X_fun(state(k, :), p )';
end
  Theta = sum2 * pinv(sum1  + delta .* eye(s1))

    %% Test the learning results;
    Learning_dis = zeros(size_case, n);
    for i =1:1:size_case
      Learning_dis(i, :) = (Theta * B_X_fun(state(i,:), p) * xi_fun(time_stamp(i, 1), p))';
    end

    figure(1)
    plot(1:1:size_case, uncertainty(:,1));
    hold on 
    plot(1:1:size_case, Learning_dis(:,1));
    grid on
    legend('True','Learned');
    ylabel('Uncertainty');
    
    Learn_Error = sum(abs(uncertainty(:,1)-Learning_dis(:,1)))/size_case


