%% Visualize slopes Well to Tunnel

load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\concatenated_source_table.mat")

% pair=199;
pair=find(AFR_table.fi==3 & AFR_table.well_names=="K11" & AFR_table.tunnel_names=="G10-u1");
row=AFR_table(pair,:);

for i=1:3000
    if i==2900
        break
    end
    if ~isnan(row.slope{1}(i)) && row.slope{1}(i)>0.1 && length(nonzeros(row.well_AFR{1}(i:i+99)))>4 && i<3000-100
        well=row.well_AFR{1}(i:i+99);
        tunnel=row.tunnel_AFR{1}(i:i+99);
        scatter(well,tunnel)
        hold on
        plot([0:100],row.slope{1}(i)*[0:100]+row.intercept{1}(i))
        hold off
        title("T="+string(i/10)+" to "+string((i+99)/10)+ " s")
        subtitle("slope="+string(row.slope{1}(i))+" rsq="+string(row.rsq{1}(i)))
        xlabel("Well AFR")
        ylabel("Tunnel AFR")
        xlim([min(well),max(well)])
        ylim([min(tunnel),max(tunnel)])
        saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim" + ...
            "\AFR Output Full Index slide 10sLM\linear regressions\FID 3 K11 G10\t "+string(i*100)+" ms.png")
        close
    end
end

%% Visualize slopes Tunnel to Well

load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\concatenated_target_table.mat")

% pair=199;
pair=find(AFR_table.fi==3 & AFR_table.well_names=="K11" & AFR_table.tunnel_names=="M7-u1");
row=AFR_table(pair,:);

for i=1:3000
    if i==2900
        break
    end
    if ~isnan(row.slope{1}(i)) && row.slope{1}(i)>0.1 && length(nonzeros(row.well_AFR{1}(i:i+99)))>4 && i<3000-100
        well=row.well_AFR{1}(i:i+99);
        tunnel=row.tunnel_AFR{1}(i:i+99);
        scatter(tunnel,well,"filled",SizeData=50)
        hold on
        plot([-10:1000],row.slope{1}(i)*[-10:1000]+row.intercept{1}(i),LineWidth=2)
        hold off
        title("T="+string(i/10)+" to "+string((i+99)/10)+ " s",FontSize=20)
        subtitle("slope="+string(row.slope{1}(i))+" rsq="+string(row.rsq{1}(i)),FontSize=16)
        xlabel("Tunnel AFR","FontSize",16)
        ylabel("Well AFR","FontSize",16)
        ylim([-10,max(well)])
        xlim([-10,max(tunnel)])
        saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim" + ...
            "\AFR Output Full Index slide 10sLM\linear regressions\FID 3 G8 H12 FB\t "+string(i*100)+" ms.png")
        close
    end
end

%% Plot all well to tunnel

load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\concatenated_source_table.mat")

% pair=199;
% pair=find(AFR_table.fi==3 & AFR_table.well_names=="K11" & AFR_table.tunnel_names=="G10-u1");
% row=AFR_table(pair,:);


AFR_table=AFR_table(AFR_table.fi==3,:);
AFR_table=AFR_table([44,91,58,33,102,58],:);

for r=2:size(AFR_table,1)
    row=AFR_table(r,:);
    for i=1:5:3000
        if i>=2900
            break
        end
        if ~isnan(row.slope{1}(i)) && row.slope{1}(i)>0.1 && length(nonzeros(row.well_AFR{1}(i:i+99)))>4 && i<3000-100
            well=row.well_AFR{1}(i:i+99);
            tunnel=row.tunnel_AFR{1}(i:i+99);
            to_keep=find(well~=0 | tunnel~=0); %keep nonzero elements in both
            hold on
            plot([-10:1000],row.slope{1}(i)*[-10:1000]+row.intercept{1}(i),'LineWidth',4,'Color','r')
            plot([-10:1000],[-10:1000],"--",'LineWidth',4,'Color',[192/255,192/255,192/255])
            scatter(well(to_keep),tunnel(to_keep),"b","filled")
            hold off
            title("T="+string(i/10)+" to "+string((i+99)/10)+ " s",'FontSize',24)
            subtitle("slope="+string(round(row.slope{1}(i),2))+" rsq="+string(round(row.rsq{1}(i),2)),'FontSize',18)
            xlabel("Well AFR (Hz)")
            ylabel("Tunnel AFR (Hz)")
            max_AFR=max([well,tunnel]);
            xlim([-10,max_AFR])
            ylim([-10,max_AFR])
            axis square
            ax=gca;
            ax.FontSize=16;
            if max_AFR>200
                xticks([0,50,100,150,200,250,300,350,400,450,500,550])
                yticks([0,50,100,150,200,250,300,350,400,450,500,550])
            elseif max_AFR>100 && max_AFR<=200
                xticks([0,20,40,60,80,100,120,140,160,180,200])
                yticks([0,20,40,60,80,100,120,140,160,180,200])
            else
                xticks([0,10,20,30,40,50,60,70,80,90,100])
                yticks([0,10,20,30,40,50,60,70,80,90,100])
            end
            mydir="D:\Brewer lab data\HFS\Temporal Analysis\No Stim" + ...
                "\AFR Output Full Index slide 10sLM\linear regressions\FID " + row.fi + "\" + row.well_names+" "+row.tunnel_names;
            if ~exist(mydir,'dir')
                mkdir(mydir)
            end
            saveas(gcf,mydir+"\t "+string(i*100)+" ms.png")
            close
        end
    end
    close all
end
%% Plot all tunnel to well

load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\concatenated_target_table.mat")

% pair=199;
% pair=find(AFR_table.fi==3 & AFR_table.well_names=="K11" & AFR_table.tunnel_names=="G10-u1");
% row=AFR_table(pair,:);


AFR_table=AFR_table(AFR_table.fi==3,:);
AFR_table=AFR_table([39],:);

for r=1%:size(AFR_table,1)
    row=AFR_table(r,:);
    for i=1:5:3000
        if i>=2900
            break
        end
        if ~isnan(row.slope{1}(i)) && row.slope{1}(i)>0.1 && length(nonzeros(row.well_AFR{1}(i:i+99)))>4 && i<3000-100
            well=row.well_AFR{1}(i:i+99);
            tunnel=row.tunnel_AFR{1}(i:i+99);
            to_keep=find(well~=0 | tunnel~=0); %keep nonzero elements in both
            
            hold on
            plot([-10:1000],row.slope{1}(i)*[-10:1000]+row.intercept{1}(i),'LineWidth',4,'Color','r')
            plot([-10:1000],[-10:1000],"--",'LineWidth',4,'Color',[192/255,192/255,192/255])
            scatter(tunnel(to_keep),well(to_keep),"b","filled")
            hold off
            title("T="+string(i/10)+" to "+string((i+99)/10)+ " s")
            subtitle("slope="+string(round(row.slope{1}(i),2))+" rsq="+string(round(row.rsq{1}(i),2)))
            xlabel("Tunnel AFR (Hz)")
            ylabel("Well AFR (Hz)")
            max_AFR=max([well,tunnel]);
            xlim([-10,max_AFR])
            ylim([-10,max_AFR])
            axis square
            ax=gca;
            ax.FontSize=16;
            if max_AFR>200
                xticks([0,50,100,150,200,250,300,350,400,450,500,550])
                yticks([0,50,100,150,200,250,300,350,400,450,500,550])
            elseif max_AFR>100 && max_AFR<=200
                xticks([0,20,40,60,80,100,120,140,160,180,200])
                yticks([0,20,40,60,80,100,120,140,160,180,200])
            else
                xticks([0,10,20,30,40,50,60,70,80,90,100])
                yticks([0,10,20,30,40,50,60,70,80,90,100])
            end
            mydir="D:\Brewer lab data\HFS\Temporal Analysis\No Stim" + ...
                "\AFR Output Full Index slide 10sLM\linear regressions\FID " + row.fi + "\" + row.tunnel_names+" "+row.well_names;
            if ~exist(mydir,'dir')
                mkdir(mydir)
            end
            saveas(gcf,mydir+"\t "+string(i*100)+" ms.png")
            close
        end
    end
    close all
end