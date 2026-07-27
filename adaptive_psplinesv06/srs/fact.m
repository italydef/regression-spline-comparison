function y = fact(x)

% Returns x!
% Syntax: y = fact(x)

n=max(size(x));
y=ones(size(x));
for i=1:n
if x(i) ~= 0
y(i)=prod(1:x(i));
end
end
