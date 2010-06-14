function C = sq_dist(a, b, Q)
%% Computes the squared distance between a and b
% Based on the algorithm from D.H. Hald, S.E. Grützmeier, S. Harder, C.I.
% Blücher: An ICA algorithm based on MCMC

if nargin < 1 || nargin > 3 || nargout > 1
    error('Wrong number of arguments.');
end

if nargin == 1 || isempty(b)         % input arguments are taken to be
    b = a;                           % identical if b is missing or empty
end

[D, n] = size(a);
[d, m] = size(b);
if d ~= D
    error('Error: column lengths must agree.');
end

if nargin < 3
    C = zeros(n,m);
    for d = 1:D
        C = C + (repmat(b(d,:), n, 1) - repmat(a(d,:)', 1, m)).^2;
    end
else
    if [n m] == size(Q) %#ok<BDSCA>
        C = zeros(D,1);
        for d = 1:D
            C(d) = sum(sum((repmat(b(d,:), n, 1) - repmat(a(d,:)', 1, m)).^2.*Q));
        end
    else
        error('Third argument has wrong size.');
    end
end

end
