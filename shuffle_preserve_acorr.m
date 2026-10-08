function x_surrogate = shuffle_preserve_acorr(x)
    % Ensure x is a column vector
    x = x(:);
    N = length(x);
    
    % 1. Transform to frequency domain
    X = fft(x);
    
    % 2. Generate random phases (preserving symmetry for real-valued signals)
    if mod(N, 2) == 0
        % Even length
        rand_phases = rand(N/2 - 1, 1) * 2 * pi;
        phases = [0; rand_phases; 0; -flipud(rand_phases)];
    else
        % Odd length
        rand_phases = rand((N - 1)/2, 1) * 2 * pi;
        phases = [0; rand_phases; -flipud(rand_phases)];
    end
    
    % 3. Apply random phases to the Fourier amplitudes
    X_shuffled = abs(X) .* exp(1i * (angle(X) + phases));
    
    % 4. Transform back to the time domain (real part discards numerical imaginary noise)
    x_surrogate = real(ifft(X_shuffled));
end