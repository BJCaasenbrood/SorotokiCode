function y = destack(x,varargin)

if isempty(varargin)
    Ndim = 2;
else
    Ndim = varargin{1};
    if Ndim < 0
        error(['Destacking number should be positive, N = ',num2str(Ndim)]);
    end
end

Nc = numel(x)/Ndim;
if ~isflint(Nc)
    error(['Length of x should be a multiple of N = ',num2str(Ndim)]);
end

y = zeros(Nc,Ndim);

for ii = 1:Ndim
    y(:,ii) = x(ii:Ndim:end);
end

end

function    isf = isflint( m )
%   floating double only
%   
%   http://www.mathworks.com/company/newsletters/news_notes/pdf/Fall96Cleve.pdf
%   JanSimon: http://www.mathworks.se/matlabcentral/answers/67247-about-isinteger- ...
%       command-confusion#answer_78914   
    assert( isa( m, 'double' ) && isvector( m )     ...
        ,   'isflint:IllegalInput'                  ...
        ,   'The input should be a double vector'   )
    isf =   all( abs(   m ) <= flintmax ) ...
        &&  all( floor( m ) == m        ) ; 
end