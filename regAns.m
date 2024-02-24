function [W_E1] = regAns(W_TO_guess)
    A = 0.0933;
    B = 1.0329;
    W_E1 = 10 .^ ((log10(W_TO_guess) - A) / B);
end
