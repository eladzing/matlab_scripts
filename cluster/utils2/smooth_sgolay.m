function [smoothed_data n_range] = smooth_sgolay(data, ORDER, R_FRAC)
% Savitzky-Golay smoothening with adaptive window size (window size
% proportional to the radius. Useful on noisy profiles.

if (R_FRAC == 0)
    smoothed_data = data;
    n_range = 1:length(data);
    return;
end

SG2 = zeros([1 floor(length(data)/2)]);
for n = 1:length(data)
    WindowLength = floor(n*R_FRAC);
    WindowLength = WindowLength - (1-mod(WindowLength, 2));
    
    if (WindowLength < 3)
        SG2(n) = data(n);
        continue;
    end
    
    HalfWin  = ((WindowLength+1)/2) -1;
    if ((n + HalfWin) > length(data))
        smoothed_data = SG2(1:n-1);
        n_range = 1:n-1;
        return;
    end
    
    [b g] = sgolay(ORDER, WindowLength);
    SG2(n) = dot(g(:,1), data(n - HalfWin: n + HalfWin));
end    

error('Unreachable Code Reached!')