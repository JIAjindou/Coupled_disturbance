%% The function of $\bm{\mathcal{B}}\left(\bm{x}\right)$
%input p-decide the approximate accuracy
function B_X = B_X_fun(x, p)
 n = length(x);
 m = 1;   % The matched disturbance
 s1 = (p+1)^(m+n);
 s2 = (p+1)^(m);
 PI = zeros((p+1)^(n), 1);
 B_X = zeros(s1, s2);
 
 for i = 1:1:(p+1)^(n)
     k_str     = dec2base(i-1,(p+1));
     k_length  = length(k_str);
     PI(i,1) = 1;
     for j = 1:1:k_length
      PI(i,1) =  PI(i,1)*Cheby_poly(str2double(k_str(k_length+1-j)), x(j));
     end
 end
 
 for j = 1:1:s2
    B_X(:, j) = [zeros(1, (j-1)*(p+1)^n) PI' zeros(1, s1-j*(p+1)^n)];
 end
end