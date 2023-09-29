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

sqr_data_size = sqrt(data_size);
time_end = 4;
state_begin = -2;
state_end = 2;

xname = {'1', '2', '3', '4', '5', '6', '7', '8'};
yname = {'100','10','1','0.1','0.01','0.001','0.0001','0.00001' };

Learning_error_fun1 = zeros(8, 8); 
Learning_error_fun2 = zeros(8, 8); 
Learning_error_fun3 = zeros(8, 8); 

MAE = zeros(3, 1);
RMSE = zeros(3, 1);

for iii = 1:1:8
    var_noise = 10^(-iii + 3);
for jjj = 1:1:8
     p = jjj;
for fun_k = 1:1:3

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
    MAE(fun_k) = sum(abs(error))/len_tes;
    RMSE(fun_k) = sqrt(error'*error/len_tes);
end
    Learning_error_fun1(iii, jjj) = MAE(1);
    Learning_error_fun2(iii, jjj) = MAE(2);
    Learning_error_fun3(iii, jjj) = MAE(3);
end
end

figure(1)
subplot(1, 3, 1)
h1 = heatmap(xname, yname, Learning_error_fun1);
subplot(1, 3, 2)
h2 = heatmap(xname, yname, Learning_error_fun2);
subplot(1, 3, 3)
h3 = heatmap(xname, yname, Learning_error_fun3);
colormap(gca, 'parula')

%% The tesing function
function uncertainty = fun_nonlinear(x, time, k)

if k == 1
%% function 1  p=5 ×î¼Ñ
d = time;
uncertainty = sin(x) .* sin(d);
elseif k == 2
%% function 2  p=3 ×î¼Ñ
d = time.^2;
uncertainty = x - 1/12 .* x.^3 - 1/4 .* d;
elseif k == 3
%% function 3  p=5 ×î¼Ñ
d = time.^3;  
uncertainty = -1/9 .* sin(x) .* d;
end

end


