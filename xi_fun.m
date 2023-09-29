%% The function of $\bm{\mathcal{B}}\left(\bm{x}\right)$
%input p-decide the approximate accuracy
function xi = xi_fun(time_stamp, p)
 m = 1;
 s2 = (p+1)^(m);
 xi = zeros(s2, 1);
 
 for i = 1:1:s2
      xi(i,1) =  Cheby_poly(i-1, time_stamp);
 end
 
end