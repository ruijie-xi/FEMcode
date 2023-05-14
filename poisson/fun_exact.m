function result = fun_exact(coords,dx,dy)
x = coords(1,:); y = coords(2,:);
if dx==0 && dy==0
    result = exp(x+y);
elseif dx==1 && dy==0
    result = exp(x+y);
elseif dx==0 && dy==1
    result = exp(x+y);
else
    warning('cannot handle this type of derivative!');
end