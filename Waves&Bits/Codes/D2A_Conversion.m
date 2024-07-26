function audio = D2A_Conversion(bitstream ,n)

scalefactor = 2^(n-1);

binary_vector_r = char(bitstream + '0');

binary_matrix = reshape(binary_vector_r, [], n); % Reshape back to matrix
audio_integers = bin2dec(binary_matrix); % Convert binary to decimal
audio_reconstructed = typecast(uint16(audio_integers), 'int16'); % Typecast to int16
audio_reconstructed_normalized = double(audio_reconstructed) / (scalefactor-1); % Normalize to [-1, 1]
audio = audio_reconstructed_normalized;

end