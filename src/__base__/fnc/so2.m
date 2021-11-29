%so(3) Lie algebra of the group SO(3)
function S = so2(t, deg)
    
    assert((isreal(t)) | isa(t, 'sym'), ...
        'SMTB:rotz:badarg', 'theta must be a real scalar or symbolic');

    if nargin > 1 && strcmp(deg, 'deg')
        t = t *pi/180;
    end
    
    tx = t;
    
%     % make almost zero elements exactly zero
%     if ~isa(t, 'sym')
%         if abs(st) < eps
%             st = 0;
%         end
%         if abs(ct) < eps
%             ct = 0;
%         end
%     end
    
    % create the rotation matrix
    S = [0  -tx;
         tx   0];
end