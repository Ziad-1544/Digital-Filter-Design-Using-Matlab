function mse = calc_mse(y, x)


    if length(y) ~= length(x) % lengths must be equal
        error('Both vectors must be same lengths');
    end
    
    mse = (1/(length(y))) * sum((y - x).^2);

end
