function bitstream = A2D_Conversion_Quant(audiofile,n)

scalefactor = 2^(n-1);

audio = audiofile(:,1);
mxa = max(abs(audio));
audio = audio + mxa;
audio = audio/max(abs(audio));
audio = audio*scalefactor;

bitstream  = dec2bin(audio,n) - '0';
bitstream = bitstream';
bitstream = bitstream(:)';

end