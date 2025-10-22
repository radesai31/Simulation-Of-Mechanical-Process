clear all;
tstart = 0;
tend = 1;
N = 100
dt = (tend-tstart)/N

t = tstart: dt : tend; % of array of times

f_ty = @()
