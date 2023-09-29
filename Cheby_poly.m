%% The Chebyshev polynominal
function poly = Cheby_poly(order, var)
  if order == 0
    poly = 1;
  elseif order == 1
    poly = var;
%    elseif order == 10
%     poly = 512*var^10 - 1280*var^8 + 1120*var^6 -400 *var*4 + 50 *var^2 -1;
  else
    poly = 2*var*Cheby_poly(order-1, var)- Cheby_poly(order-2, var);
  end
end



