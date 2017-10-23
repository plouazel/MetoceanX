function plotSomeMetData(met,WaveType)
colorstr={[0 0 1], [1 0 1],[1 0 0],[1 1 0],[0 1 0],[0 1 1],[.6 .6 0],[0 .6 .6]};

if ~exist('WaveType','var')
    WaveType='combo';
end
if strcmp(WaveType,'combo')
    metfield = 'all';
    WaveStr='';
elseif strcmp(WaveType,'chop')
    metfield = WaveType;
    WaveStr = 'Wind-Sea';
elseif strcmp(WaveType,'swell')
    metfield = WaveType;
    WaveStr = 'Swell';
end
nMet = length(met);
%%%%%%%%%%%%%%%%%%%%% ---------- Hs v Tp -----------%%%%%%%%%%%%%%%%%%%%
    figure('name',[WaveStr ' H_s vs T_p'])
    str1='$H_s$ vs $T_p$' ;
    hold on
    for ii=1:nMet
        plot(met(ii).(metfield).Tp,met(ii).(metfield).Hs,'linestyle','none','marker','.','markeredgecolor',colorstr{ii},'markersize',8-ii)%'markerfacecolor',colorstr{ii}
        legstr1{ii} = met(ii).name;
    end
    tp1s=2:.01:10;
    hsbreak=(tp1s/1.4/7.07).^2*9.81;
    plot(tp1s,hsbreak,'g-')
    legend(legstr1{:})
    titlestr=sprintf('%s   ',str1);
    grid on
    xlabel([WaveStr 'Wave Period (s)'])
    ylabel([WaveStr  'Wave Height (m)'])
    title(titlestr,'FontSize',12,'interpreter','latex','FontWeight','demi')
    hold off
    %print( [figdir sprintf('HsTp_N%d%s',nbins,ixstr) '.png'],'-dpng','-r300')
%%%%%%%%%%%%%%%%%%%%% ---------- Vspd v Vdir -----------%%%%%%%%%%%%%%%%%%%%
    
    %polar(pi,max(met(ii).all.vspd))
    str1='$V_{wind}$ (m/s)   vs  $\theta_{wind}$';
    figure('name', 'WindSpd vs WindDir')
    for ii=1:nMet
        h1=polar(met(ii).all.vdir*pi/180,met(ii).all.vspd,'k.');
        set(h1,'linestyle','none','marker','.','markeredgecolor',colorstr{ii},'markersize',8-2*ii);
                hold on
    end
    grid on
    hold off
    view([90 -90])
    legend(legstr1{:})
    titlestr=sprintf('%s       , %s which the wind is blowing',str1,'FROM');
    %set(gca,'thetaticks',[0:45:315],'thetaticklabels',{'N','NE','E','SE','S','SW','W','NW'})
    ht=title(titlestr,'FontSize',12,'interpreter','latex','FontWeight','demi');
    hAxes2=gca;
    axlimits=axis;
    currentAxesPosition=get(hAxes2,'Position'); % Get the axes' (polar plot) current position
    set(hAxes2,'Position',currentAxesPosition+[0 -0.04 0 0]);   % Shift that position downward
    title2=get(hAxes2,'Title'); % Get handle to the title object
    currentTitlePosition=get(title2,'Position');    % Get the title's current position
    set(title2,'Position',currentTitlePosition+[abs(axlimits(1))*.1 0 0]); % Shift the position upward
%%%%%%%%%%%%%%%%%%%%% ---------- Hs v Wdir -----------%%%%%%%%%%%%%%%%%%%%
    figure('name', [WaveStr ' WaveHs vs WaveDir'])
    str1=['$H_{s}$   vs    $\theta_{' WaveStr '}$'];
    for ii=1:nMet
        h1=polar(met(ii).(metfield).wdir*pi/180,met(ii).(metfield).Hs,'k.');
        set(h1,'linestyle','none','marker','.','markeredgecolor',colorstr{ii},'markersize',8-2*ii)
        hold on
    end 
    grid on
    hold off
    legend(legstr1{:})
    
    view([90 -90])
    titlestr=sprintf('%s      , %s which the waves travel',str1,'TOWARDS');
    title(titlestr,'FontSize',12,'interpreter','latex','FontWeight','demi')
        hAxes2=gca;
    axlimits=axis;
    currentAxesPosition=get(hAxes2,'Position'); % Get the axes' (polar plot) current position
    set(hAxes2,'Position',currentAxesPosition+[0 -0.04 0 0]);   % Shift that position downward
    title2=get(hAxes2,'Title'); % Get handle to the title object
    currentTitlePosition=get(title2,'Position');    % Get the title's current position
    set(title2,'Position',currentTitlePosition+[abs(axlimits(1))*.1 0 0]); % Shift the position upward
    %print( [figdir sprintf('FatigueBins_HsWth_N%d%s',nbins,ixstr) '.png'],'-dpng','-r300')
%%%%%%%%%%%%%%%%%%%%% ---------- Vspd v Wdir -----------%%%%%%%%%%%%%%%%%%%%
    
%     figure('name', [ ' WindSpd vs ' WaveStr ' WaveDir'])
%     polar(pi,max(met(ii).all.vspd));
%     str1=['$V_{wind}$   vs   $H_{s ' WaveStr '}$'];    
%     for ii=1:nMet
%         h2=polar(met(ii).(metfield).wdir*pi/180,met(ii).all.vspd,'k.');
%         set(h2,'linestyle','none','marker','.','markeredgecolor',colorstr{ii},'markersize',8-2*ii);
%         hold on
%     end 
%     grid on
%     hold off
%     ylabel('Wind Speed (m/s)')
%     xlabel([WaveStr  'Wave Height (m)'])
%     legend(legstr1{:})
%     titlestr=sprintf('%s      , %s which the waves travel',str1,'TOWARDS');
%     title(titlestr,'FontSize',12,'interpreter','latex','FontWeight','demi');
%     view([90 -90])
%     hAxes2=gca;
%     axlimits=axis;
%     currentAxesPosition=get(hAxes2,'Position'); % Get the axes' (polar plot) current position
%     set(hAxes2,'Position',currentAxesPosition+[0 -0.04 0 0]);   % Shift that position downward
%     title2=get(hAxes2,'Title'); % Get handle to the title object
%     currentTitlePosition=get(title2,'Position');    % Get the title's current position
%     set(title2,'Position',currentTitlePosition+[abs(axlimits(1))*.1 0 0]);
    
%%%%%%%%%%%%%%%%%%%%% ---------- Hs v Vspd -----------%%%%%%%%%%%%%%%%%%%%
    figure('name',[WaveStr ' H_s vs Vspd'])
    str1='$H_s$ vs $V_{wind}$' ;
    hold on
    for ii=1:nMet
        plot(met(ii).all.vspd,met(ii).(metfield).Hs,'linestyle','none','marker','.','markeredgecolor',colorstr{ii},'markersize',8-ii)%'markerfacecolor',colorstr{ii}
        legstr1{ii} = met(ii).name;
    end
    legend(legstr1{:})
    titlestr=sprintf('%s   ',str1);
    grid on
    xlabel(['Wind Speed (m/s)'])
    ylabel([WaveStr  'Wave Height (m)'])
    title(titlestr,'FontSize',12,'interpreter','latex','FontWeight','demi')
    hold off
    
end
