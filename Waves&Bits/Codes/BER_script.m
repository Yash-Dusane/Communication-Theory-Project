clc, clearvars
[audio, fs] = audioread('project.wav');
sgtit = {"Memoryless Channel";"Memory Channel";"Rectangular Pulse"; "Raised Cosine Pulse";"'BER vs Eb/N0 :"};

n = 16;
bitstream = A2D_Conversion(audio, n);

encoded_stream = Encoder(bitstream);

symbols = encoded_stream;%(1:2000);
nsymbols = length(symbols);

m = 4;
rof = 0.25;
len = 10;


% -------------------------- Inputs

transmit_filter_choice = 1;   % 0 : Raised Cosine, 1 : Rectangular
dispersive = 0;               % 0 : Memoryless,    1 : Memory

% --------------------------

if(~transmit_filter_choice), [transmit_filter,dummy] = raised_cosine(rof,m,len); end
if(transmit_filter_choice), transmit_filter = ones(m,1); end

transmit_output = LineCoding(symbols, transmit_filter, m);

fc = 1000000;
t2 = (cumsum(ones(size(transmit_output)))-1)/m;
modulating_sig=cos(2*pi*fc*t2/(fs*n/2));

transmit_mod_op = Amp_Modulation(transmit_output, modulating_sig);

a = 0.5;
b = 4;
Tb = m;

EbN0_dB_range = 0:1:30; 
EbN0_range = 10.^(EbN0_dB_range/10);
Eb=1;

err_pr = zeros(numel(EbN0_range),1);

for hvc = 1:numel(EbN0_range)

N0 = Eb/EbN0_range(hvc);
sigma = sqrt(N0/2);
variance = sigma^2;

[noisy_rec_ip, channel] = Channel_Tx(transmit_mod_op, dispersive, sigma, a, b, Tb);
modulating_sig = [modulating_sig; zeros(numel(noisy_rec_ip) - numel(modulating_sig),1)];
t2 = (cumsum(ones(size(noisy_rec_ip)))-1)/m;
outp = Demodulation(noisy_rec_ip, modulating_sig, fs, n);
[receive_output, rx_sampling_times, rx_sampled_outputs, rx_symbols] = LineDecoding(outp, transmit_filter, m, channel, nsymbols);

err_pr(hvc) =(sum(symbols~=rx_symbols))/numel(rx_symbols);

end

% rx_output_bits = Decoder(rx_symbols);
% audio_rec = D2A_Conversion(rx_output_bits, n);

figure(1)
semilogy(EbN0_dB_range,err_pr/2,'LineWidth',2);
ylabel("P(error)");
xlabel("Eb/N0");
title([sgtit{5}, sgtit{(numel(transmit_filter)~=m)+3}, sgtit{dispersive+1}])

