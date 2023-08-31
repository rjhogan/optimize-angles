% This Matlab script plots the simulations previously performed in the
% "run" directory, producing some of the figures 3 and 4 of the Hogan
% paper

directory='../run';
dataset='evaluation2';
scenario='present';
supspace = ' '; % Fix bad PDF rendering in some Matlab installations
input_file = [directory '/optical-depth/ecckd-1.0_' dataset '_lw_climate_fsck-32b_optical-depth_' scenario '.nc'];
ref_file = [directory '/fluxes/' dataset '_gauss-jacobi-5_64_fluxes.nc'];
quadratures = {'gauss-legendre','gauss-laguerre','gauss-jacobi-5',...
	       'optimized-default','optimized-default-ratio','elsasser','lacis'};
% Legends
leg_full = {'Gauss-Legendre','Gauss-Laguerre','Gauss-Jacobi-5','Optimized',...
       'Optimized-IR','Elsasser','Lacis & Oinas'}
leg = {'Legendre','Laguerre','Jacobi-5','Optimized',...
       'Optimized-IR','Elsasser','Lacis'}
stys = {'-','-','-','-','--','-','-'};
colours = [0 0 0;
	   0 0 1;
	   1 0 0;
	   0 0.7 0;
	   0.3 1 0;
	   1 0 1;
	   0 0.9 0.9];
reorder = 1:length(quadratures);
show_hr_error_bar = [1 0 1 0 0 0 0];
streams_offset = [0 0 0 0 1 0 2];

in = loadnc(input_file);
ref= loadnc(ref_file);

nlay = size(ref.heating_rate,1);
hr_weight = (sqrt(in.pressure_hl(2:end,:))-sqrt(in.pressure_hl(1:end-1,:))) ...
	    ./ (ones(nlay,1)*sqrt(in.pressure_hl(end,:)));
pmid = 0.5.*(in.pressure_hl(2:end,:)+in.pressure_hl(1:end-1,:));
% Tropospheric weight
index = find(pmid<100e2);
hr_weight_trop = hr_weight;
hr_weight_trop(index) = 0.0;
hr_weight_trop = hr_weight_trop./(ones(nlay,1)*sum(hr_weight_trop));
% Stratospheric weight
index = find(pmid>=100e2);
hr_weight_strat = hr_weight;
hr_weight_strat(index) = 0.0;
hr_weight_strat = hr_weight_strat./(ones(nlay,1)*sum(hr_weight_strat));

for iquad = 1:length(quadratures)
  qdata = loadnc([directory '/quadrature/quadrature_' quadratures{iquad} '.nc']);
  orders=double(qdata.order(find(qdata.order<=16)));
  if length(orders)>1
    is_one(iquad) = 0;
    if orders(1) > 1
      syms{iquad} = '*--';
    else
      syms{iquad} = 'o-';
    end
  else
    is_one(iquad) = orders(1);
    syms{iquad} = '*';
  end
  streams{iquad} = orders.*2;
  for iord = 1:length(orders)
    file = [directory '/fluxes/' dataset '_' quadratures{iquad} '_' ...
		      num2str(orders(iord)) '_fluxes.nc'];
    data = loadnc(file);
    % Error variances
    toa_errvar{iquad}(iord) = mean((data.flux_up_toa-ref.flux_up_toa).^2);
    surf_errvar{iquad}(iord)= mean((data.flux_dn_surf-ref.flux_dn_surf).^2);
    hr_errvar{iquad}(iord)= mean(sum(hr_weight.*(data.heating_rate-ref.heating_rate).^2,1));
    hr_errvar_trop{iquad}(iord)= mean(sum(hr_weight_trop.*(data.heating_rate-ref.heating_rate).^2,1));
    hr_errvar_strat{iquad}(iord)= mean(sum(hr_weight_strat.*(data.heating_rate-ref.heating_rate).^2,1));
    hr_errvar_lay{iquad}(:,iord)= mean(hr_weight.*(data.heating_rate-ref.heating_rate).^2,2);
    bias_hr{iquad}(:,iord) = mean(data.heating_rate-ref.heating_rate,2);
    std_hr{iquad}(:,iord)  = std(data.heating_rate-ref.heating_rate,0,2);
    ci_hr{iquad}(:,iord)   = std_hr{iquad}(:,iord).*1.96;
  end
end

flux_weight = 0.02;

for iquad = 1:length(quadratures)
  flux_errstd{iquad} = sqrt(toa_errvar{iquad}+surf_errvar{iquad});
  toa_errstd{iquad} = sqrt(toa_errvar{iquad});
  surf_errstd{iquad} = sqrt(surf_errvar{iquad});
  hr_errstd{iquad}   = sqrt(hr_errvar{iquad});
  hr_errstd_trop{iquad}  = sqrt(hr_errvar_trop{iquad});
  hr_errstd_strat{iquad} = sqrt(hr_errvar_strat{iquad});
  tot_errstd{iquad}  = sqrt(flux_weight.*(toa_errvar{iquad}+surf_errvar{iquad}) ...
			    +hr_errvar{iquad});
end

% Figure 3 of the Hogan paper
figure(1)
clf
set(gcf,'paperposition',[0.5 0.5 27 20],'defaultlinelinewidth',1)

subplot(2,3,1)
for iquad = 1:length(quadratures)
  loglog(streams{iquad}, flux_errstd{iquad}, syms{iquad},'color',colours(iquad,:));
  hold on
end
xlabel('Number of streams')
ylabel(['RMSE in surface/TOA irradiances (W m^{' supspace '-2})'])
xlim([1.8 36])
ylim([1e-7 1e1]);set(gca,'ytick',10.^[-8:1])
set(gca,'xtick',[2 4 8 16 32],'xminortick','off','xminorgrid','off');
set(gca,'yminortick','on','yminorgrid','off')
grid on
legend(leg,'location','southwest')
text(0,1.01,'\bf(a) Irradiances','units','normalized','verticalalignment','bottom')

subplot(2,3,2)
for iquad = 1:length(quadratures)
  loglog(streams{iquad}, hr_errstd_trop{iquad}, syms{iquad},'color',colours(iquad,:));
  hold on
end
xlabel('Number of streams')
ylabel(['RMSE in heating rate,  {\itp} >100 hPa (K d^{' supspace '-1})'])
xlim([1.8 36])
ylim([1e-8 1e-1]);set(gca,'ytick',10.^[-8:1])
set(gca,'xtick',[2 4 8 16 32],'xminortick','off','xminorgrid','off');
set(gca,'yminortick','on','yminorgrid','off')
grid on
text(0,1.01,'\bf(b) Tropospheric heating rates','units','normalized','verticalalignment','bottom')

subplot(2,3,3)
for iquad = 1:length(quadratures)
  loglog(streams{iquad}, hr_errstd_strat{iquad}, syms{iquad},'color',colours(iquad,:));
  hold on
end
xlabel('Number of streams')
ylabel(['RMSE in heating rate,  {\itp} <100 hPa (K d^{' supspace '-1})'])
xlim([1.8 36])
set(gca,'xtick',[2 4 8 16 32],'xminortick','off','xminorgrid','off');
set(gca,'yminortick','on','yminorgrid','off')
grid on
text(0,1.01,'\bf(c) Stratospheric heating rates','units','normalized','verticalalignment','bottom')
ylim([1e-7 1e0]);set(gca,'ytick',10.^[-8:1])

% For graphical abstract
figure(10)
clf
set(gcf,'paperposition',[0.5 0.5 27 20],'defaultlinelinewidth',1)

subplot(2,3,1)
for iquad = 1:length(quadratures)
  loglog(streams{iquad}, flux_errstd{iquad}, syms{iquad},'color',colours(iquad,:));
  hold on
end
xlabel('Number of streams')
ylabel(['RMS error in clear-sky fluxes (W m^{' supspace '-2})'])
xlim([1.8 36])
ylim([1e-7 1e1]);set(gca,'ytick',10.^[-8:1])
set(gca,'xtick',[2 4 8 16 32],'xminortick','off','xminorgrid','off');
set(gca,'yminortick','on','yminorgrid','off')
grid on
pos=get(gca,'position')
legend(leg_full,'location','eastoutside')
drawnow
set(gca,'position',pos);


pmidmed = median(pmid,2)./100;

gg=[0.25 0.25 0.25;
    0 0 1;
    1 0 0;
    0 1 0;
    0 1 1];

% Figure 4 of the Hogan paper
figure(2)
clf
set(gcf,'paperposition',[0.5 0.5 32 10]);
xspan = [1 0.2 0.05];
set(gcf,'defaultlinelinewidth',1.5);
for io = 1:3
  h=subplot(1,3,io);
  pos=get(h,'position');
  set(h,'position',pos+[0 0.05 0 -0.05]);
  if io == 1
    for iq=1:length(bias_hr)
      plot(-1,-1,stys{iq},'color',colours(iq,:))
      hold on
    end
  end
  for iq = 1:length(bias_hr)
    iqq = reorder(iq);
    if show_hr_error_bar(iqq)
      h=fill([bias_hr{iqq}(:,io)-ci_hr{iqq}(:,io);...
	      flip(bias_hr{iqq}(:,io)+ci_hr{iqq}(:,io))],...
	     [pmidmed;flip(pmidmed)],'r','facecolor',gg(iqq,:),'edgecolor','none');
      set(h,'facealpha',0.2);
      hold on
    end
  end
   for iq = 1:length(bias_hr)
    iqq = reorder(iq);
    if ~is_one(iqq) & show_hr_error_bar(iqq)
      plot(bias_hr{iqq}(:,io)-ci_hr{iqq}(:,io),pmidmed,'k:','color',gg(iqq,:),'linewidth',0.5)
      plot(bias_hr{iqq}(:,io)+ci_hr{iqq}(:,io),pmidmed,'k:','color',gg(iqq,:),'linewidth',0.5)
    end
   end
  for iq = 1:length(bias_hr)
    iqq = reorder(iq);
    if ~is_one(iqq) | io==is_one(iqq)
      if is_one(iqq)
	semilogy(bias_hr{iqq}(:,1),pmidmed,stys{iqq},'color',colours(iqq,:));
      else
	ind=io-streams_offset(iqq);
	if ind > 0
	  semilogy(bias_hr{iqq}(:,io-streams_offset(iqq)),pmidmed,stys{iqq},'color',colours(iqq,:));
	end
      end
    end
  end
  set(gca,'ydir','reverse');set(gca,'yscale','log')
  ylim([0.01 1e3]);
  xlim(xspan(io).*[-1 1]);
  set(gca,'xtick',xspan(io).*[-1 -0.5 0 0.5 1])
  grid on
  plot([0 0],[1e-2 1e3], 'k','linewidth',0.5)
  xlabel('Heating rate bias (K d^{-1})');
  ylabel('Pressure (hPa)')
  if io == 1
    legend(leg,'location','southeast')
  end
  set(gca,'yminorgrid','off','layer','top')
  set(gca,'Xticklabelrotationmode','manual')
  title(['(' 'a'+io-1 ') ' num2str(io.*2) ' streams']);
end

