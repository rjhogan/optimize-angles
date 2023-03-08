quads = {'gauss-legendre','gauss-laguerre','gauss-jacobi-5','optimized-fw0.02','optimized-fw0.02-ratio4'};
only_two = [0 0 0 0 1];
for iq = 1:length(quads)
  quadrature{iq} = loadnc(['../run/quadrature_' quads{iq} '.nc']);
end
cols = ['kbrgc'];
clf
set(gcf,'defaultlinelinewidth',1.5);
for io = 1:8
  subplot(4,2,io)
  for iq = 1:length(quads)
    q=quadrature{iq};
    if only_two(iq)
      if io==2
	for in = 1:io
	  plot(q.mu(in,1).*[1 1],q.weight(in,1).*[0 1],cols(iq));
	  hold on
	end
      end
    else
      for in = 1:io
	plot(q.mu(in,io).*[1 1],q.weight(in,io).*[0 1],cols(iq));
	hold on
      end
    end
  end
  xlim([0 1]);
end
