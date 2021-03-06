function foo = popultationsPlot(t, X)
%{
    Input: X the soltion generated by ode45
    Return all the populations plots vs time in matrix plot format
    
%}
    Ms = X(:, 1); M1 = X(:, 2); M2 = X(:, 3); 
    Is = X(:, 4); I1 = X(:, 5); I2 = X(:, 6); 
    Sm1 = X(:, 7); Sm2 = X(:, 8);
    Ym1 = X(:, 9); Ym2 = X(:, 10); 
    Rec = X(:, 11);
    z = X(:,12); w = X(:,13);

    figure(1)
    hold off
    subplot(4, 3, 1)
    plot(t, Ms,'k')
    title('Ms')

    subplot(4, 3, 2)
    plot(t, M1,'k')
    title('M1')

    subplot(4, 3, 3)
    plot(t, M2,'k')
    title('M2')

    subplot(4, 3, 4)
    plot(t, Is,'k')
    title('Is')

    subplot(4, 3, 5)
    plot(t, I1,'k')
    title('I1')

    subplot(4, 3, 6)
    plot(t, I2,'k')
    title('I2')

    subplot(4, 3, 7)
    plot(t, Sm1,'k')
    title('Sm1')

    subplot(4, 3, 8)
    plot(t, Sm2,'k')
    title('Sm2')

    subplot(4, 3, 9)
    plot(t, Ym1,'k')
    title('Ym1')

    subplot(4, 3, 10)
    plot(t, Ym2,'k')
    title('Ym2')

    subplot(4, 3, 11)
    plot(t, Rec,'k')
    title('R')
    %% Dengue Clasic Fever and Dengue Hemorragic Fever likening
    figure(2)
    hold off
    plot(t, z,'b',t, w,'r')
    legend('CD','HDF')
    xlabel('t')

    %% Dengue Clasic Fever and Dengue Hemorragic Fever likening
    figure(3)
    hold off
    
    % Load data sets of DF and HDF for 2010 in Hermosillo, Sonora

    datFD=load('dengue_c_her2010.dat');
    diaFD=datFD(:,1)-40386;
    casosFD=datFD(:,2);
    plot(diaFD, casosFD,'b*');
    grid on
    hold on

    datFHD=load('dengue_h_her2010.dat');
    diaFHD=datFHD(:,1)-40384;
    casosFHD=datFHD(:,2);

    plot(diaFHD,casosFHD,'ro')
    hold on
    plot(t, z,'b',t, w, 'r')
    legend('CD','HDF')
    xlabel('t')
    axis([0 160 0 70])
    
    title('Confirmed DF and DHF cases, Hermosillo 2010')
    xlabel('date')
    ylabel('Incidence')
    legend('FD','FHD')




foo = true;
