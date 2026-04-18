function C_d = sigmoid(t, tmax)
% Returns sigmoidal signal climbing from 0 to 1 between 0 and tmax.
% inputs:  t: time,   tmax: duration of increase

if t < 0
    C_d = 0;
elseif t > tmax
    C_d = 1;
else
    C_d = (1 - cos(t / tmax * pi)) / 2;
end
end