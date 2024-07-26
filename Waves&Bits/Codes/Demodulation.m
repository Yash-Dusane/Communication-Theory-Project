function outp = Demodulation(noisy_rec_ip, carrier, fs, n)
    out_r = noisy_rec_ip.*carrier;
    outp = 2*lowpass(out_r,fs,fs*n/2);
end