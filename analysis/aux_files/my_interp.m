function y = my_interp(vector, factor)

y = [];
for m = 1 : length(vector)
    
    x = vector (m);
    y = cat(2, y, repmat(x, 1, factor));
    
end

y = y';
