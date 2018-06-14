% %
function varargout = locpoly(y, x, x_0,degree, derivative, width, kernel)
% Local polynomial regression
% Author: Liu Wei(117020208005@2017.swufe.edu.cn)
% Version: 2017-12-03
% univariate local polynomial regression, fit regression function and
% any order-derivatives of regression function.
% Local Polinomial Smoothing (LPS) for 1-dim curve fitting.
%
%Syntax:
%[m0, dm0, d2m0,...] = locpoly(y, x, x_0,degree, derivative, width, kernel)
%
% Input:
% y -- dependent variable.
% x -- independent variable.
% x_0 -- points of interests (optional), default as x_0 = x.
% degree -- degree of polynomial(optional), default as degree = 1,that is local linear.
% derivative -- derivative of interest(optional), default as 0, that is regression function itself.
% width -- bandwidth of kernel function(optional), default as optimal
% bandwidth of normal kernel or 1.
% kernel -- kernel function (optional), default as kernel = 'normal', with options: 
% normal
% epanechinikov
% box
% triangle.
%
% Output:
% m0 -- fitted values at x_0, if derivative =0 , only return m0.
% dm0 -- fitted first derivative values at x_0, if derivative=1, return m0 and dm0.
% dKm0 -- fitted K-order derivative values at x_0, if derivative = K,
% return m0,dm0,...,dKm0.
%
% Ref:
% Fan (1996) Local Polynomial Modelling And Its Application.

if(size(x,2) ~= 1)
    x = x';
end
if(size(y,2) ~= 1)
    y =y';
end
if(~exist('x_0', 'var') || isempty(x_0))
   x_0 = x; 
end
if(~exist('degree', 'var') ||  isempty(degree))
    degree = 1;
end
if(~exist('derivative', 'var') || isempty(derivative))
   derivative = 0; 
end
% Default window parameter is optimal for normal distribution
n = length(y);
if (~exist('width', 'var') ||  isempty(width))
    med = median(x);
    sig = median(abs(x-repmat(med,n,1))) / 0.6745;
    if sig>0
        width = sig * (4/(3*n))^(1/6) / 2;
    else
        width= 1;
    end
end

if(~exist('kernel', 'var')  || isempty(kernel))
   kernel = 'normal'; 
end
if(nargout > (derivative+1))
   error('output parameters can not exceed derivative+1 !') 
end
if(degree < derivative)
   error('degree must be greater than derivative!') 
end
% start local regression
n0 = length(x_0);
Alpha = zeros(n0, degree+1);
for i = 1:n0
    Alpha(i,:) = pointwise_reg(y, x, x_0(i), kernel, width, degree);
end
if(degree >= derivative)
    for j = 1:(derivative+1)
        varargout{j} = Alpha(:,j);
    end
end

% -----------------------------
% The following are functions that define smoothing kernels k(z).
% Each function takes a single input Z and returns the value of
% the smoothing kernel.  These sample kernels are designed to
% produce outputs that are somewhat comparable (differences due
% to shape rather than scale), so they are all probability
% density functions with unit variance.
%
% The density estimate has the form
%    f(x;k,h) = mean over i=1:n of k((x-y(i))/h) / h

function f = normal(z)
%NORMAL Normal density kernel.
%f = normpdf(z);
f = exp(-0.5 * z .^2) ./ sqrt(2*pi);

function f = epanechinikov(z)
%EPANECHINIKOV Epanechinikov's asymptotically optimal kernel.%'
a = sqrt(5);
z = max(-a, min(z,a));
f = .75 * (1 - .2*z.^2) / a;

function f = box(z)
%BOX    Box-shaped kernel
a = sqrt(3);
f = (abs(z)<=a) ./ (2 * a);

function f = triangle(z)
%TRIANGLE Triangular kernel.
a = sqrt(6);
z = abs(z);
f = (z<=a) .* (1 - z/a) / a;

% pointwise estimates for E(Y|X=x0) and its derivatives, including higher
% order derivatives
function halpha = pointwise_reg(y, x, x0, kernel, width, degree)
n = length(y);
X = zeros(n, degree+1);
X(:,1) = ones(n,1);
for j = 1:degree
    X(:,j+1) = (x-x0).^j;
end
W = diag(feval(kernel, (x-x0)/width)/width);
halpha = (X'*W*X)\(X'*W*y);
