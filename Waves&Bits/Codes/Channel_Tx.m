function [noisy_rec_ip, channel] = Channel_Tx(transmit_mod_op, dispersive, sigma, a, b, Tb)


if dispersive == 0, channel = 1;
else
channel = [a; zeros(b*Tb-1,1); 1-a];
end

receive_input = conv(transmit_mod_op,channel);

noisy_rec_ip  = receive_input + sigma*randn(numel(receive_input), 1);

end