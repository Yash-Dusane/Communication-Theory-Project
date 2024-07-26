function [receive_output, rx_sampling_times, rx_sampled_outputs, rx_symbols] = LineDecoding(outp, transmit_filter, m, channel, nsymbols)

receive_filter = flipud(transmit_filter);
receive_output = (1/m)*conv(outp,receive_filter,'same');

pulse = conv(transmit_filter,channel,'same');
rx_pulse = conv(pulse,receive_filter,'same')/m;

[~, maxloc] = max(rx_pulse);
rx_sampling_times = maxloc:m:(nsymbols-1)*m+maxloc;
rx_sampled_outputs = receive_output(rx_sampling_times);

rx_symbols=zeros(numel(rx_sampled_outputs),1);

for klo = 1:numel(rx_sampled_outputs)
    if(rx_sampled_outputs(klo) <0.5)
        rx_symbols(klo)=0;
    end
    if(rx_sampled_outputs(klo) <1.5 & rx_sampled_outputs(klo)>=0.5)
        rx_symbols(klo)=1;
    end
    if(rx_sampled_outputs(klo) <2.5 & rx_sampled_outputs(klo)>=1.5)
        rx_symbols(klo)=2;
    end
    if(rx_sampled_outputs(klo) >=2.5)
        rx_symbols(klo)=3;
    end
end

end