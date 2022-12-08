function h = timetitle(t)
%TIMETITLE set the current figure title to time T (s): h = timetitle(T)
set(gcf,'Name',['T = ',num2str(t,3),' (s)']);
end

