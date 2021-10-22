function [varargout] = getSvH(H,n)

assert(n == nargout);

H = sum(H, 3);
sv = svd(H);

validNOut = min(nargout,length(sv));
for i = 1:validNOut
    varargout{i} = sv(i);
end
for i = validNOut+1:nargout
    varargout{i} = NaN;
end

end