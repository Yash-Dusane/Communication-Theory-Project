clc, clearvars
[audio, fs] = audioread('project.wav');

n = 16;
scalefactor = 2^(n-1);

audio_normalized = int16(audio(:,1)/max(abs(audio(:,1))) * (scalefactor-1)); % Normalize
audio_binary = dec2bin(typecast(audio_normalized(:), 'uint16'), n); % Convert to binary
binary_vector = audio_binary(:)'; % Reshape to vector for transmission

bitstream = binary_vector - '0';

gain = 1;
signal = {1,4};

encoded_stream = Encoder(bitstream);
symbols = gain*encoded_stream;
nsymbols = length(symbols);

oversampling_factor = 4;
m = oversampling_factor;

rof = 0.25;
length = 10;

[transmit_filter,dummy] = raised_cosine(rof,m,length);
%transmit_filter = ones(m,1);
dispersive = 0;



nsymbols_upsampled = 1+(nsymbols-1)*m;
symbols_upsampled = zeros(nsymbols_upsampled,1);
symbols_upsampled(1:m:nsymbols_upsampled)=symbols;

transmit_output = (conv(symbols_upsampled,transmit_filter));
figure(1)
subplot(4,1,1), stem(bitstream(1:2000)); title('Input Bitstream');
legend('x_1')


figure(1)
subplot(4,1,2), stem(symbols); title('Symbols generated from Encoder');
legend('x_2')


t2 = (cumsum(ones(size(transmit_output)))-1)/m;
figure(1);
subplot(4,1,3)
plot(t2,transmit_output);
xlabel('t/T'); title("Line coded sequence");
legend('x_3')
 % signal = transmit_output;
 % [pxx, f] = periodogram(signal, hamming(numel(signal)), numel(signal), fs);
 % 
 %    % Return PSD
 %    psd = pxx;
 % 
 %    % Plot PSD (optional)
 %    figure(3);
 %    subplot(4,1,1)
 %    plot(f, 10*log10(psd));
 %    title('Power Spectral Density of x3 (line coded s/g)');
 %    xlabel('Frequency (Hz)');
 %    ylabel('Power/Frequency (dB/Hz)');
signal{1} = transmit_output;
autocorr_result = xcorr(signal{1});

fft_autocorr = fft(autocorr_result);
nm = numel(fft_autocorr)-1;
freq = (-nm/2:nm/2)*fs*n/(2*(nm+1));
figure(4);
subplot(411),plot(freq, fftshift(10*log10(abs(fft_autocorr))));
xlabel('Frequency');
ylabel('Magnitude (dB)');
title('Sx_3(f)')
sgtitle('PSD (on log scale) - Rectangular Pulse, Memory Channel');

fc = 1000000;
modulating_sig = zeros(numel(t2),1);
for klo = 1:numel(t2)
    modulating_sig(klo)=cos(2*pi*fc*t2(klo)/(fs*n/2));
end

transmit_mod_op = modulating_sig.*transmit_output;

figure(1);
subplot(4,1,4)
plot(t2,transmit_mod_op);
title("modulated sequence");
legend('x_4')
    % signal2 = transmit_mod_op;
    % [pxx2, f2] = periodogram(signal2, hamming(numel(signal2)), numel(signal2), fs);
    % 
    % Return PSD
    % psd2 = pxx2; 
    % Plot PSD (optional)
    % figure(3);
    % subplot(4,1,2)
    % plot(f2, 10*log10(psd2));
    % title('Power Spectral Density of x4 (modulated s/g)');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

signal2 = transmit_mod_op;
autocorr_result = xcorr(signal2);

fft_autocorr = fft(autocorr_result);
nm = numel(fft_autocorr)-1;
freq = (-nm/2:nm/2)*fs*n/(2*(nm+1));
figure(4);
subplot(412),plot(freq, fftshift(10*log10(abs(fft_autocorr))));
xlabel('Frequency');
ylabel('Magnitude (dB)');
title('Sx_4(f)')


a = 0.5;
b = 4;
Tb = fs/1000;

if dispersive == 0, channel = 1;
else
channel = [a; zeros(b*Tb-1,1); 1-a];
end

% receive_input = transmit_output;
receive_input = conv(transmit_mod_op,channel);


% EbN0_dB_range = 0:0.25:15; 
% EbN0_range = 10.^(EbN0_dB_range/10);
% 
% % prob=zeros(numel(EbN0_range),1);
% Eb=1;
% 
% err_pr = zeros(numel(EbN0_range),1);
% for hvc = 1:numel(EbN0_range)
% 
% N0 = Eb/EbN0_range(hvc);
% sigma = sqrt(N0/2);
% variance = sigma^2;



sigma = 0.1;
noisy_rec_ip  = receive_input + sigma*randn(numel(receive_input), 1);

t2 = (cumsum(ones(size(noisy_rec_ip)))-1)/m;
figure(2);
subplot(4,1,1)
plot(t2,noisy_rec_ip);
title("noisy modulated signal");
legend('y_4')

    % signal3 = noisy_rec_ip;
    %  [pxx3, f3] = periodogram(signal2, hamming(numel(signal2)), numel(signal2), fs);
    % % Return PSD
    % psd3 = pxx3;
    % 
    % % Plot PSD (optional)
    % figure(3);
    % subplot(4,1,3)
    % plot(f3, 10*log10(psd3));
    % title('Power Spectral Density of y4 (noisy modulated s/g)');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    signal3 = noisy_rec_ip;
autocorr_result = xcorr(signal3);

fft_autocorr = fft(autocorr_result);
nm = numel(fft_autocorr)-1;
freq = (-nm/2:nm/2)*fs*n/(2*(nm+1));
figure(4);
subplot(413),plot(freq, fftshift(10*log10(abs(fft_autocorr))));
xlabel('Frequency');
ylabel('Magnitude (dB)');
title('Sy_4(f)')

modulating_sig = [modulating_sig; zeros(numel(noisy_rec_ip) - numel(modulating_sig),1)];

out_r = noisy_rec_ip.*modulating_sig;
outp = 2*lowpass(out_r,48000,48000*n/2);
figure(2);
subplot(4,1,2)
plot(t2,outp);
title("demodulated signal");
legend('y_3')
% signal4 = outp;
%  [pxx4, f4] = periodogram(signal4, hamming(numel(signal4)), numel(signal4), fs);
%     % Return PSD
%     psd4 = pxx4;
%     % Plot PSD (optional)
%     figure(3);
%     subplot(4,1,4)
%     plot(f4, 10*log10(psd4));
%     title('Power Spectral Density of y3 (demodulated s/g)');
%     xlabel('Frequency (Hz)');
%     ylabel('Power/Frequency (dB/Hz)');

signal4 = outp;
autocorr_result = xcorr(signal4);

fft_autocorr = fft(autocorr_result);
nm = numel(fft_autocorr)-1;
freq = (-nm/2:nm/2)*fs*n/(2*(nm+1));
figure(4);
subplot(414),plot(freq, fftshift(10*log10(abs(fft_autocorr))));
xlabel('Frequency');
ylabel('Magnitude (dB)');
title('Sy_3(f)')

%outp = receive_input;

receive_filter = flipud(transmit_filter);
% figure(4), plot(receive_filter)
receive_output = (1/m)*conv(outp,receive_filter,'same');
t3 = (cumsum(ones(size(receive_output)))-1)/m;
figure(2);
subplot(4,1,3)
plot(t3,receive_output,'b');
xlabel('t/T'); title("line decoding")
hold on; 
% legend('y_2')
pulse = conv(transmit_filter,channel,'same');
rx_pulse = conv(pulse,receive_filter,'same')/m;
[maxval, maxloc] = max(rx_pulse);
rx_sampling_times = maxloc:m:(nsymbols-1)*m+maxloc;
rx_sampled_outputs = receive_output(rx_sampling_times);
subplot(4,1,3),stem((rx_sampling_times-1)/m,rx_sampled_outputs,'r');
legend('y_2')
hold off;

rx_sampled_outputs = rx_sampled_outputs/gain;
figure(5), scatter(real(rx_sampled_outputs),imag(rx_sampled_outputs));
title('Constellation Plot : Rectangular Pulse, Memory Channel')
xlabel('Real'), ylabel('Imag')
rx_output_bits=zeros(numel(rx_sampled_outputs),1);

for klo = 1:numel(rx_sampled_outputs)
    if(rx_sampled_outputs(klo) <0.5)
        rx_output_bits(klo)=0;
    end
    if(rx_sampled_outputs(klo) <1.5 & rx_sampled_outputs(klo)>=0.5)
        rx_output_bits(klo)=1;
    end
    if(rx_sampled_outputs(klo) <2.5 & rx_sampled_outputs(klo)>=1.5)
        rx_output_bits(klo)=2;
    end
    if(rx_sampled_outputs(klo) >=2.5)
        rx_output_bits(klo)=3;
    end
end



figure(2);
subplot(4,1,3)
stem(rx_output_bits,'r'), title("Symbol decoding");
legend('y_2')
N = numel(rx_output_bits);
count=0;var=0;
for k = 1:N
    var = symbols(k)/gain-rx_output_bits(k);
    if(var~=0)
        count=count+1;
       % disp(k);
    end
end

err_pr=count/N;

disp(['error probability is ', num2str(err_pr)]);

final = zeros(2*numel(rx_output_bits),1);
j=1;
g=1;

while j <= numel(rx_output_bits)
if(rx_output_bits(j)>=3)
    final(g)=1;
    final(g+1)=1;
end
if(rx_output_bits(j)>=2 && rx_output_bits(j)<3)
    final(g)=1;
    final(g+1)=0;
end
if(rx_output_bits(j)>=1 && rx_output_bits(j)<2)
    final(g)=0;
    final(g+1)=1;
end
if(rx_output_bits(j)<1)
    final(g)=0;
    final(g+1)=0;
end
j=j+1;
g=g+2;
end

% while j <= numel(final) - 8
% fgo=bin2dec(binstr(final(j:j+8)));
% j = j+8;
% end
%filename = 'j.wav';
%audiowrite(filename,fgo,fs);

binary_vector_r = char(final + '0');

binary_matrix = reshape(binary_vector_r, [], n); % Reshape back to matrix
audio_integers = bin2dec(binary_matrix); % Convert binary to decimal
audio_reconstructed = typecast(uint16(audio_integers), 'int16'); % Typecast to int16
audio_reconstructed_normalized = double(audio_reconstructed) / (scalefactor-1); % Normalize to [-1, 1]

figure(2);
subplot(414),stem(final); title("Decoded Bitstream");
legend('y_1');



sound(audio_reconstructed_normalized,fs);
