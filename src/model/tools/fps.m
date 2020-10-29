function A = fps(time,FPS)
dt = mean(diff(time));
A = max(round((1/FPS)/dt),1);
end

