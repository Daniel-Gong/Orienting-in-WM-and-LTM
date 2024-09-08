function rng=myRange(x)
    %replacement for stats toolbox range function
    rng = max(x) - min(x);
end