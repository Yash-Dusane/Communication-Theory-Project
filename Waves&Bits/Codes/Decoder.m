function rx_output_bits = Decoder(rx_symbols)

rx_output_bits = zeros(2*numel(rx_symbols),1);
j=1;
g=1;

while j <= numel(rx_symbols)
if(rx_symbols(j)>=3)
    rx_output_bits(g)=1;
    rx_output_bits(g+1)=1;
end
if(rx_symbols(j)>=2 & rx_symbols(j)<3)
    rx_output_bits(g)=1;
    rx_output_bits(g+1)=0;
end
if(rx_symbols(j)>=1 & rx_symbols(j)<2)
    rx_output_bits(g)=0;
    rx_output_bits(g+1)=1;
end
if(rx_symbols(j)<1)
    rx_output_bits(g)=0;
    rx_output_bits(g+1)=0;
end

j=j+1;
g=g+2;

end

end