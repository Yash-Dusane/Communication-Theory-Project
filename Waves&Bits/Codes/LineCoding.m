function transmit_output = LineCoding(symbols, transmit_filter, m)

nsymbols = length(symbols);

nsymbols_upsampled = 1+(nsymbols-1)*m;
symbols_upsampled = zeros(nsymbols_upsampled,1);
symbols_upsampled(1:m:nsymbols_upsampled)=symbols;

transmit_output = (conv(symbols_upsampled,transmit_filter));

end