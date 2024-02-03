function ind = resample(w)
    % Perform stochastic resampling Ripley (1988)
    % Implementation by G. Hendeby
    
    N = numel(w);
    qs = cumsum(w);
    u = fliplr(cumprod(rand(1,N).^(1 ./ (N:-1:1))));

    i = 1;
    ind = zeros(1, N);
    for p = 1:N
        while qs(i) < u(p)
            i = i + 1;
        end
            ind(p) = i;
    end
end