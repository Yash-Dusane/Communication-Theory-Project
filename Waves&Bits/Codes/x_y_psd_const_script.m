clc, clearvars
[audio, fs] = audioread('project.wav');
sgtit = {"Memoryless Channel";"Memory Channel";"Rectangular Pulse"; "Raised Cosine Pulse";"PSD (on log scale) - "; "Constellation Plot : "};

n = 16;
bitstream = A2D_Conversion(audio, n);

encoded_stream = Encoder(bitstream);

symbols = encoded_stream;%(1:2000);
nsymbols = length(symbols);

m = 4;
rof = 0.25;
length = 10;


% -------------------------- Inputs

transmit_filter_choice = 0;   % 0 : Raised Cosine, 1 : Rectangular
dispersive = 0;               % 0 : Memoryless,    1 : Memory

% --------------------------

if(~transmit_filter_choice), [transmit_filter,dummy] = raised_cosine(rof,m,length); end
if(transmit_filter_choice), transmit_filter = ones(m,1); end

transmit_output = LineCoding(symbols, transmit_filter, m);

figure(1)
subplot(4,1,1), stem(bitstream); title('Input Bitstream'); legend('x_1')
subplot(4,1,2), stem(symbols); title('Symbols generated from Encoder'); legend('x_2')
sgtitle([ sgtit{(numel(transmit_filter)~=m)+3}, sgtit{dispersive+1}])

t2 = (cumsum(ones(size(transmit_output)))-1)/m;
subplot(4,1,3), plot(t2,transmit_output); xlabel('t/T'); title("Line coded sequence"); legend('x_3')
signal{1} = transmit_output;
% autocorr_result = xcorr(signal1;)
% 
% fft_autocorr = fft(autocorr_result);
% nm = numel(fft_autocorr)-1;
% freq = (-nm/2:nm/2)*fs*n/(2*(nm+1));
% figure(4);
% subplot(411),plot(freq, fftshift(10*log10(abs(fft_autocorr))));
% xlabel('Frequency');
% ylabel('Magnitude (dB)');
% title('Sx_3(f)')
% sgtitle('PSD (on log scale) - Rectangular Pulse, Memory Channel');

fc = 1000000;

modulating_sig=cos(2*pi*fc*t2/(fs*n/2));

transmit_mod_op = Amp_Modulation(transmit_output, modulating_sig);

subplot(4,1,4),plot(t2,transmit_mod_op); title("modulated sequence"); legend('x_4')
signal{2} = transmit_mod_op;
% autocorr_result = xcorr(signal2);
% 
% fft_autocorr = fft(autocorr_result);
% nm = numel(fft_autocorr)-1;
% freq = (-nm/2:nm/2)*fs*n/(2*(nm+1));
% figure(4);
% subplot(412),plot(freq, fftshift(10*log10(abs(fft_autocorr))));
% xlabel('Frequency');
% ylabel('Magnitude (dB)');
% title('Sx_4(f)')


a = 0.5;
b = 4;
Tb = m;
sigma = 0.01;

[noisy_rec_ip, channel] = Channel_Tx(transmit_mod_op, dispersive, sigma, a, b, Tb);

figure(2);
sgtitle([ sgtit{(numel(transmit_filter)~=m)+3}, sgtit{dispersive+1}])
t2 = (cumsum(ones(size(noisy_rec_ip)))-1)/m;
subplot(4,1,1),plot(t2,noisy_rec_ip);title("noisy modulated signal");legend('y_4')

modulating_sig = [modulating_sig; zeros(numel(noisy_rec_ip) - numel(modulating_sig),1)];
signal{3} = noisy_rec_ip;
% autocorr_result = xcorr(signal3);
% 
% fft_autocorr = fft(autocorr_result);
% nm = numel(fft_autocorr)-1;
% freq = (-nm/2:nm/2)*fs*n/(2*(nm+1));
% figure(4);
% subplot(413),plot(freq, fftshift(10*log10(abs(fft_autocorr))));
% xlabel('Frequency');
% ylabel('Magnitude (dB)');
% title('Sy_4(f)')

outp = Demodulation(noisy_rec_ip, modulating_sig, fs, n);


subplot(4,1,2),plot(t2,outp); title("demodulated signal");legend('y_3')
signal{4} = outp;
% autocorr_result = xcorr(signal4);
% 
% fft_autocorr = fft(autocorr_result);
% nm = numel(fft_autocorr)-1;
% freq = (-nm/2:nm/2)*fs*n/(2*(nm+1));
% figure(4);
% subplot(414),plot(freq, fftshift(10*log10(abs(fft_autocorr))));
% xlabel('Frequency');
% ylabel('Magnitude (dB)');
% title('Sy_3(f)')


[receive_output, rx_sampling_times, rx_sampled_outputs, rx_symbols] = LineDecoding(outp, transmit_filter, m, channel, nsymbols);


t3 = (cumsum(ones(size(receive_output)))-1)/m;
subplot(4,1,3),plot(t3,receive_output,'b');xlabel('t/T'); title("line decoding")
hold on; 
subplot(4,1,3),stem((rx_sampling_times-1)/m,rx_sampled_outputs,'r');legend('y_2')
hold off;


err_pr=(sum(symbols~=rx_symbols))/numel(rx_symbols);
disp(['error probability is ', num2str(err_pr)]);
subplot(4,1,3), stem(rx_symbols,'r'), title("Symbol decoding"); legend('y_2')

rx_output_bits = Decoder(rx_symbols);

figure(2);
subplot(414),stem(rx_output_bits);title("Decoded Bitstream");legend('y_1');

audio_rec = D2A_Conversion(rx_output_bits, n);

figure(3)
tit = {'Sx_3f';'Sx_4f';'Sy_4f';'Sy_3f'};
for k = 1:4
autocorr_result = xcorr(signal{k});
fft_autocorr = fft(autocorr_result);
nm = numel(fft_autocorr)-1;
freq = (-nm/2:nm/2)*fs*n/(2*(nm+1));
subplot(4,1,k),plot(freq, fftshift(10*log10(abs(fft_autocorr))));
xlabel('Frequency');
ylabel('Magnitude (dB)');
title(tit{k})
sgtitle([sgtit{5}, sgtit{(numel(transmit_filter)~=m)+3}, sgtit{dispersive+1}]);
end

figure(4), scatter(real(rx_sampled_outputs),imag(rx_sampled_outputs));
title([sgtit{6}, sgtit{(numel(transmit_filter)~=m)+3}, sgtit{dispersive+1}])
xlabel('Real'), ylabel('Imag')


%sound(audio_rec,fs);
