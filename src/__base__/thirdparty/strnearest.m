function [index,dist]=strnearest(key,list,varargin)
% function [index,dist]=strnearest(key,list)
% based on fileexchange strdist.m
%
% compares one input string to a list of target strings using a converging
% Levenshtein algorithm, contracting the bounds of its search to the best 
% match found so far.  Returns the index into the list of the best match, 
% along with the Levenshtein distance.  If multiple matches are found at the
% same distance, index will be a vector.
%
% function [...]=strnearest(key,list,args)
% allows the user to switch the function's behavior.
% Currently recognized arguments include:
%
%  any integer: sets a threshold for the maximum distance to match.  If no
% matches are found, the function will return an empty index.
%
%  'first': stop after first exact match, or first match under threshold
%  
%  'case':  ignore case (default is to match case)
%
%  'editor': compute editor distance: substitutions have weight 2
%
% Inherently vectorized for one input; recursively vectorizes the other 
% input, returning a cell array of index vectors and a vector of distances.
%d=strdist(r,b,krk,cas) computes Levenshtein and editor distance 
%between strings r and b with use of Vagner-Fisher algorithm.
%   Levenshtein distance is the minimal quantity of character
%substitutions, deletions and insertions for transformation
%of string r into string b. An editor distance is computed as 
%Levenshtein distance with substitutions weight of 2.
%d=strdist(r) computes numel(r);
%d=strdist(r,b) computes Levenshtein distance between r and b.
%If b is empty string then d=numel(r);
%d=strdist(r,b,krk)computes both Levenshtein and an editor distance
%when krk=2. d=strdist(r,b,krk,cas) computes a distance accordingly 
 %with krk and cas. If cas>0 then case is ignored.
%
%Example.
% disp(strdist('matlab'))
%    6
% disp(strdist('matlab','Mathworks'))
%    7
% disp(strdist('matlab','Mathworks',2))
%    7    11
% disp(strdist('matlab','Mathworks',2,1))
%    6     9
%
% Modified by BK to vectorize both r and b (trivially: recurse)
thresh=NaN;
firstmatch=false;
cas=0;
krk=1;
if isa(key,'cell')
  for i=1:length(key)
    [index{i},dist(i)]=strnearest(key{i},list,varargin{:});
  end
else
  while ~isempty(varargin)
    arg=varargin{1}; varargin=varargin(2:end);
    if isnumeric(arg)
      if arg==fix(arg)
        thresh=arg;
      else
        warning('Ignoring non-integer numeric optarg')
      end
    else
      switch arg(1:3)
        case 'fir'
          firstmatch=true;
        case 'cas'
          cas=1;
        case 'edi'
          krk=2;
        otherwise
          warning(['Ignoring unrecognized optarg ' arg])
      end
    end
  end
  
  if cas>0
    key=upper(key);
  end
  index=[];
  
  if isnan(thresh)
    dist=max([length(key),length(list{1})]);
  else
    dist=thresh;
  end
  for T=1:length(list)
    if cas>0
      bb=upper(list{T});
    else
      bb=list{T};
    end
    eql=char({key;bb}); % equal length
    keyy=eql(1,:);      % space-padded key
    bb=eql(2,:);        % space-padded candidate
    luma=numel(keyy);
    dl=dist*ones([luma+1,luma+1]);
    dl(1,:)=0:luma;   dl(:,1)=0:luma;
    %Distance
    for i=1:luma
      for j=max([1,i-dist]):min([luma,i+dist])
        kr=krk*(~(keyy(min([j luma]))==bb(i)));
        dl(i+1,j+1)=min([dl(i,j)+kr,dl(i,j+1)+1,dl(i+1,j)+1]);
      end
    end
    if dl(end,end)==dist
      index=[index T];
      if firstmatch && dist<=thresh
        break
      end
    elseif dl(end,end)<dist
      index=T;
      dist=dl(end,end);
      if firstmatch && dist==0
        break
      end
    end
  end
end
