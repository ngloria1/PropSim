%% Test GCrit and FindG

clear
close all

cd('..\')
addpath(fullfile(pwd, 'Supporting Functions'))

Temperature = 275.81;
N2O_i = N2O_Properties(Temperature);

tic
%% Test GCrit
[crit_flow, mass_flux] = TwoPhaseN2OFlow(Temperature, 1);
fprintf('%.3g s for GCrit execution\n', toc)
fprintf('%.0f NaNs in oxidizer mass flux returned by GCrit\n', ...
    sum(sum(isnan(mass_flux.G))))

figure
subplot(2,1,1)
plot(crit_flow.p_up_norm, crit_flow.G_crit)
xlabel('Normalized upstream pressure (p_{upstream}/p_{vapor})')
xlim([1, 1.5])
ylabel({'Injector Oxidizer Mass Flux','(kg/m^2s)'})
subplot(2,1,2)
plot(crit_flow.p_up_norm, crit_flow.p_down_norm_crit)
xlabel('Normalized upstream pressure (p_{upstream}/p_{vapor})')
xlim([1, 1.5])
ylabel({'Critical normalized downstream pressure','(p_{upstream}/p_{vapor})'})
suptitle(sprintf('Critical Flow Conditions for T_{upstream}=%.1f K', ...
    Temperature))


%% Test FindG
n_plot = 3;
for p_up = N2O_i.Pvap*linspace(1.0,1.5, n_plot)
    p_down = (0:0.01:1)*p_up;
    for ii = 1:length(p_down)
        timer1 = tic;
        [G_test(ii),G_crit_test(ii),~] = FindG(Temperature, p_up, p_down(ii));
        t1 = toc(timer1);
    end
    
    figure
    plot(p_down./p_up,G_test,'x',p_down./p_up,G_crit_test,'o')
    xlabel('p_{down}')
    ylabel('Oxidizer Mass Flux [kg/m^2*s]')
    legend({'Actual','Critical'})
    title(sprintf('Test: Oxidizer Mass Flux Curve at p_{1,norm} = %.3g',...
        p_up/N2O_i.Pvap))
end
fprintf('%.3g s for FindG (after initial loading) execution\n', t1)

cd('.\Test Cases')