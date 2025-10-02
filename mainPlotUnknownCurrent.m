clear all

load CookAssignemnt1UnknownCurrent

figure(5);
clf;

for j = 1 : size(vStep, 2)
    
    subplot(2,1,1);
    hold on;
    grid on;
    plot(t,vStep(:,j), 'LineWidth', 2);
    axis tight;
    ylabel('V step (mV)');
    
    subplot(2,1,2);
    hold on;
    grid on;
    plot(t,iUnknownCurrent(:,j), 'LineWidth', 2);
    axis tight;
    xlabel('Time (msec)');
    ylabel('Current (mA)');
end