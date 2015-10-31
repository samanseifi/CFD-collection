function [ u_init ] = WavePacket( x, nx, k )
for i = 1:nx
	if (x(i) >= 0) && (x(i) <= 1.0)
        u_init(i) = sin(k*pi*x(i));
    else
        u_init(i) = 0.0;
    end
end

end

