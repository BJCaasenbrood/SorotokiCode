function y = dTranslate(x,move)
y = x - repmat(move(:)',length(x),1);
end

