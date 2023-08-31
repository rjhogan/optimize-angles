% Print out numbers for paper
qstr = {'gauss-legendre','gauss-laguerre','gauss-jacobi-5','optimized-fw0.02'};
qstr = {'gauss-legendre','gauss-laguerre','gauss-jacobi-5','optimized-default','optimized-default-ratio2'};
qstr = {'gauss-laguerre','gauss-jacobi-5','optimized-default','optimized-default-ratio2','optimized-default-jacobi-5-ratio-prior3'};

transform_weights = 1;

for iq = 1:length(qstr)
  qdat{iq} = loadnc(['../run/quadrature/quadrature_' qstr{iq} '.nc']);

  if transform_weights
    for iorder = 1:length(qdat{iq}.order)
      wt = qdat{iq}.weight(1:qdat{iq}.order(iorder),iorder).*qdat{iq}.mu(1:qdat{iq}.order(iorder),iorder);
      wt = wt./sum(wt);
      qdat{iq}.weight(1:qdat{iq}.order(iorder),iorder) = wt;      
    end
  end
end
qdat{4}.mu=[[0;0;0;0] qdat{4}.mu];
qdat{4}.weight=[[0;0;0;0] qdat{4}.weight];
qdat{5}.mu=[[0;0;0;0] qdat{5}.mu];
qdat{5}.weight=[[0;0;0;0] qdat{5}.weight];

% To compare Gauss-Jacobi-5 to Table 25.8 of Abramowitz and Stegun
%qdat{2}.weight = qdat{2}.weight./6;

do_row_colour = 1;
do_ratio = 1;

precision = '%10.8f';
precision = '%12.10f';

if do_row_colour
  disp('\hline');
end

for io = 1:4
  prefix = [num2str(io) ' & $\mu$ &'];
  for in = 1:io
    if io == 1
      mustr4='';
      mustr5='';
    elseif do_ratio & in > 1
      mustr4=sprintf([precision ' ($\\times$%d)'],qdat{4}.mu(in,io),round(qdat{4}.mu(in,io)./qdat{4}.mu(1,io)));
      mustr5=sprintf([precision ' ($\\times$%d)'],qdat{5}.mu(in,io),round(qdat{5}.mu(in,io)./qdat{5}.mu(1,io)));
    else
      mustr4=sprintf(precision,qdat{4}.mu(in,io));
      mustr5=sprintf(precision,qdat{5}.mu(in,io));
    end
    if do_row_colour
      disp('\rowcolor{white}');
    end
    disp(sprintf(['%s ' precision ' & ' precision ' & ' precision ' & %s & %s \\\\'], ...
		 prefix, qdat{1}.mu(in,io), qdat{2}.mu(in,io), qdat{3}.mu(in,io), mustr4, mustr5));
    prefix = '  &       &';
  end
  prefix = '  & $w$   &';
  for in = 1:io
    if do_row_colour
      disp('\rowcolor{black!10}');
    end
    disp(sprintf(['%s ' precision ' & ' precision ' & ' precision ' & ' precision ' & ' precision ' \\\\'],...
		prefix, qdat{1}.weight(in,io), qdat{2}.weight(in,io), qdat{3}.weight(in,io), qdat{4}.weight(in,io), qdat{5}.weight(in,io)));
    prefix = '  &       &';
  end
  if do_row_colour
    disp('\hline');
  end
end
