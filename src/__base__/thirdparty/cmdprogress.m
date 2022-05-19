classdef cmdprogress < handle
% class for command-line progress-bar notification.
% Use example:
%   pb = cmdprogress('Doing stuff...');
%   for k = 1 : 10
%       pb.print(k,10)
%       % do stuff
%   end
%
% Author: Itamar Katz, itakatz@gmail.com

    properties
        last_msg_len = 0;
        dot_length = 10;
    end
    
    methods
        %--- ctor
        function obj = cmdprogress(msg)
            fprintf('>> %s ', msg)
        end
        %--- print method
        function print(obj, n, tot)
            fprintf('%s', char(8*ones(1, obj.last_msg_len))) % delete last info_str
            info = sprintf('%d/%d', n, tot);
            dots = dotter(obj,n,tot);
            info_str = strcat(dots," (",info,")");
            fprintf('%s', info_str);
            %--- assume user counts monotonically
            if n == tot
                %cout('\n')
            end
            obj.last_msg_len = strlength(info_str);
        end
        %--- dtor
        function delete(~)
            fprintf('\n')
        end
        %--- dot conversion
        function str = dotter(obj,n,tot)
            maxdots = obj.dot_length;
            dots = ceil((n/tot)*maxdots);
            booldot = zeros(maxdots,1);
            booldot(1:dots) = 1;
            str = [];
            for ii = 1:maxdots
                if booldot(ii) == 1
                    str = strcat(str,".");
                else, str = strcat(str," ");
                end
            end
            str = strcat(str," ");
        end
    end
end