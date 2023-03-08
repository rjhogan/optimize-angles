% Print out numbers for paper
qstr = {'gauss-legendre','gauss-laguerre','gauss-jacobi-5','optimized-fw0.02'};
for iq = 1:length(qstr)
  qdat{iq} = loadnc(['../run/quadrature_' qstr{iq} '.nc']);
end

for io = 1:4
  prefix = [num2str(io) ' & $\mu$ &'];
  for in = 1:io
    disp(sprintf('%s %10.8f & %10.8f & %10.8f & %10.8f \\\\', ...
		prefix, qdat{1}.mu(in,io), qdat{2}.mu(in,io), qdat{3}.mu(in,io), qdat{4}.mu(in,io)));
    prefix = '  &       &';
  end
  prefix = '  & $w$   &';
  for in = 1:io
    disp(sprintf('%s %10.8f & %10.8f & %10.8f & %10.8f \\\\',...
		prefix, qdat{1}.weight(in,io), qdat{2}.weight(in,io), qdat{3}.weight(in,io), qdat{4}.weight(in,io)));
    prefix = '  &       &';
  end

end
