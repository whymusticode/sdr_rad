function [P,E] = pdfCust(DX,S)

%   if size(M,2) == 1
%     DX = X-repmat(M,1,size(X,2));  
    E = 0.5*sum(DX.*(S\DX),1);
    d = size(S,1);
    E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
    P = exp(-E);
%   elseif size(X,2) == 1
%     DX = repmat(X,1,size(M,2))-M;  
%     E = 0.5*sum(DX.*(S\DX),1);
%     d = size(M,1);
%     E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
%     P = exp(-E);
%   else
%     DX = X-M;  
%     E = 0.5*DX'*(S\DX);
%     d = size(M,1);
%     E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
%     P = exp(-E);
%   end
end