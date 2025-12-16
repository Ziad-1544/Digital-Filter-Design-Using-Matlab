function loss_percent = calc_energylost(x, y)

    L = max(length(x), length(y)); % which of the two is bigger

    % fill shorter vector with zeroes 
    x(end+1:L) = 0; % extends X to length L the rest being zeroes
    y(end+1:L) = 0; % same if Y is the shorter one

    % Compute energies
    Ex = sum(abs(x).^2); % x is original signal before filtering
    Ey = sum(abs(y).^2); % y is signal AFTER filtering

    % Energy loss percentage
    loss_percent = (1 - Ey/Ex) * 100;

end
