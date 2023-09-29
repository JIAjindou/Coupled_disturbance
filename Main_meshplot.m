clear all
clc
close all

%% Initial variable 
n = 1;  % The dimensionality of system state
m = 1;  % The dimensionality of uncertainty
data_size = 10000;
time_stamp = zeros(data_size,1);
state = zeros(data_size,1);

delta = 0.01;
var_noise = 0.01;

sqr_data_size = sqrt(data_size);
time_end = 4;
state_begin = -2;
state_end = 2;

MAE = zeros(3, 1);
RMSE = zeros(3, 1);

for fun_k = 1:1:3
    if fun_k == 2
        p = 3;
    else
        p = 5;
    end
for j = 1:1:sqr_data_size
   state((j-1)*sqr_data_size+1:j*sqr_data_size) = (state_begin + (j-1)*(state_end-state_begin)/sqr_data_size).*ones(sqr_data_size,1);
    if mod(j,2) == 1
          time_stamp((j-1)*sqr_data_size+1:j*sqr_data_size) = (0:time_end/sqr_data_size:time_end-time_end/sqr_data_size)';
    else
          time_stamp((j-1)*sqr_data_size+1:j*sqr_data_size) = (time_end-time_end/sqr_data_size:-time_end/sqr_data_size:0)';
    end
end
uncertainty = fun_nonlinear(state, time_stamp, fun_k);

%% The state is disrupted by noise

noise = wgn(data_size, 1, 10*log10(var_noise));
state_noise = state + noise;

uncertainty_noise = fun_nonlinear(state_noise, time_stamp, fun_k);

%% Training dataset and testing dataset
randIndex = randperm(size(state,1));

%Scramble the data order
state_new = state(randIndex, :);
uncertainty_new = uncertainty(randIndex,:);
time_stamp_new = time_stamp(randIndex,:);

size_train = data_size/2;                  % The used dataset size of each case for training.

state_train = state_new(1:size_train, :);
time_stamp_train = time_stamp_new(1:size_train, :);
uncertainty_train = uncertainty_new(1:size_train, :);

rand_test = sort(randIndex(size_train+1:2*size_train));

state_test = state(rand_test, :);
time_stamp_test = time_stamp(rand_test,:);
uncertainty_test = uncertainty(rand_test,:);

%Scramble the noised data order
state_noise_new = state_noise(randIndex, :);
uncertainty_noise_new = uncertainty_noise(randIndex,:);

state_noise_train = state_noise_new(1:size_train, :);
uncertainty_noise_train = uncertainty_noise_new(1:size_train, :);

state_noise_test = state_noise(rand_test, :);
uncertainty_noise_test = uncertainty_noise(rand_test,:);

%% Initial variables of optimization
s1 = (p+1)^(m+n);
s2 = (p+1)^(m);
Theta = ones(n,s1); 

%% Calcualte the closed form solution
  sum1 = zeros(s1, s1);
  sum2 = zeros(n, s1);
for k = 1:1:size_train
    xi = xi_fun(time_stamp_train(k, 1), p);
    sum1 = sum1 + B_X_fun(state_noise_train(k, :), p )* xi * xi' * B_X_fun(state_noise_train(k, :), p )';
    sum2 = sum2 + uncertainty_noise_train(k, :)' * xi' * B_X_fun(state_noise_train(k, :), p )';
end
  Theta =sum2 * pinv(sum1  + delta .* eye(s1));

    %% Test the learning results
    len_tes = size(state_noise_test, 1);
    Learning_uncer = zeros(len_tes, 1);
    for i =1:1:len_tes
      Learning_uncer(i, 1) = (Theta * B_X_fun(state_noise_test(i,:), p) * xi_fun(time_stamp_test(i, 1), p))';
    end
     
    error = uncertainty_noise_test-Learning_uncer;
    MAE(fun_k) = sum(abs(error))/len_tes
    RMSE(fun_k) = sqrt(error'*error/len_tes)
    
    %% plot3
    figure(2)
    subplot(2,3,fun_k)
    plot3(time_stamp_test, state_test, uncertainty_test);
    grid on
%     legend('True');
    zlabel('Uncertainty');
    xlabel('time');
    ylabel('state');
    
    subplot(2,3,fun_k+3)
    plot3(time_stamp_test, state_noise_test, Learning_uncer);
    grid on
%     legend('Learned');
    zlabel('Uncertainty');
    xlabel('time');
    ylabel('state');
    hold on
    
    %% mesh plot
    step_mesh = 0.05;
    t = 0:step_mesh:time_end;
    x = state_begin:step_mesh:state_end;
    [time_mesh, state_mesh] = meshgrid(t,x);
%     uncertainty_mesh = fun_nonlinear(state_mesh, time_mesh, fun_k);
    uncertainty_mesh = zeros(size(state_mesh,1),size(time_mesh,1));
    learned_uncertainty_mesh = zeros(size(state_mesh,1),size(time_mesh,1));
    
    for ii = 1:1:size(state_mesh,1)
%         noise_current = wgn(1, 1, 10*log10(var_noise));
        noise_current = 0;
        state_current = state_begin + (ii-1)*step_mesh + noise_current;
        for jj =  1:1:size(time_mesh,1)
             time_current = (jj-1)*step_mesh;
             uncertainty_mesh(ii, jj) = fun_nonlinear(state_current, time_current, fun_k);
             learned_uncertainty_mesh(ii, jj) = (Theta * B_X_fun(state_current, p) * xi_fun(time_current, p))';
        end
    end
     
    figure(3)
    subplot(2,3,fun_k)
    mesh(time_mesh, state_mesh, uncertainty_mesh);
    zlabel('$\Delta$');
    xlabel('T');
    ylabel('x');
    
    subplot(2,3,fun_k+3)
    mesh(time_mesh, state_mesh, learned_uncertainty_mesh);
    zlabel('$\Delta$');
    xlabel('T');
    ylabel('x');
    hold on
    
end
    
%% The tesing function
function uncertainty = fun_nonlinear(x, time, k)

if k == 1
%% function 1  p=5 
d = time.^3;  
uncertainty = -1/9 .* sin(x) .* d;
elseif k == 2
%% function 2  p=3 
d = time.^2;
uncertainty = x - 1/12 .* x.^3 - 1/4 .* d;
elseif k == 3
%% function 3  p=5 
d = time;
uncertainty = sin(x) .* sin(d);
end

end


