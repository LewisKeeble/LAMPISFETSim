function y = BuMeters(x, whichway, dim, N)
    % if whichway = 0, bu to meters
    % if whichway = 1, meters to bu
    if(whichway==0)
        y = x*dim/N;
    else
        y = x*N/dim;
    end
end