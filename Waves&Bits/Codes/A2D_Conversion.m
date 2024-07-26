function bitstream = A2D_Conversion(audio,n)

scalefactor = 2^(n-1);

audio_normalized = int16(audio(:,1)/max(abs(audio(:,1))) * (scalefactor-1)); % Normalize
audio_binary = dec2bin(typecast(audio_normalized(:), 'uint16'), n); % Convert to binary
binary_vector = audio_binary(:)'; % Reshape to vector for transmission

bitstream = binary_vector - '0';

end