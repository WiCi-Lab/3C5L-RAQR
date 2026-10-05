function N = make_odd_grid_size(value,min_value)
N=max(min_value,round(value));
if mod(N,2)==0, N=N+1; end
end
