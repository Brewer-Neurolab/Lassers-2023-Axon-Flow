function [forward_t_table,backward_t_table]=pre_post_stim_paired_ttest_210225(table1ff,table1fb,table2ff,table2fb)

% Assumes tables 1 and 2 are already organized

if all(~ismember({'Electrode Pairs','Spike Area','Conduction Time','Unit Pairs','Direction'}, table1ff.Properties.VariableNames))...
        && all(~ismember({'Electrode Pairs','Spike Area','Conduction Time','Unit Pairs','Direction'}, table2ff.Properties.VariableNames))
    error('This table does not have the appropriate variables: Electrode Pairs,Spike Area,Conduction Time,Unit Pairs,Direction.')
end

sub_conx=unique(table1ff.Subregion);

t=tiledlayout(2,length(sub_conx));

for i=1:length(sub_conx)*2
    ax(i)=subplot(2,5,i);
end

% fft2t=table('Subregion','h','p','ci','stats');
% fbt2t=table('Subregion','h','p','ci','stats');
% fftrt=table('Subregion','h','p','ci','stats');
% fbtrt=table('Subregion','h','p','ci','stats');
% fftlt=table('Subregion','h','p','ci','stats');
% fbtlt=table('Subregion','h','p','ci','stats');

for i=1:length(sub_conx)
    if length([table1ff.("Spike Area"){table1ff.Subregion==sub_conx(i)}])~=length([table2ff.("Spike Area"){table2ff.Subregion==sub_conx(i)}])
        error('The tables conain a different number of spike area elements. Try zero padding the data.')
    end
    if length([table1fb.("Spike Area"){table1fb.Subregion==sub_conx(i)}])~=length([table2fb.("Spike Area"){table2fb.Subregion==sub_conx(i)}])
        error('The tables conain a different number of spike area elements. Try zero padding the data.')
    end
    
    range=table1ff.Subregion==sub_conx(i);
    
    t1_ff_spikes=[table1ff.("Spike Area"){range}];
    t2_ff_spikes=[table2ff.("Spike Area"){range}];
    
    if (~isempty(t1_ff_spikes)&&~isempty(t2_ff_spikes))&&length(t1_ff_spikes)>2
        [hff(i),pff(i),ciff(i,:),statsff{i}]=ttest(t1_ff_spikes,t2_ff_spikes);
        [hffr(i),pffr(i),ciffr(i,:),statsffr{i}]=ttest(t1_ff_spikes,t2_ff_spikes,'Tail','right');
        [hffl(i),pffl(i),ciffl(i,:),statsffl{i}]=ttest(t1_ff_spikes,t2_ff_spikes,'Tail','left');
        
        if isnan(hff(i))
            hff(i)=0;
        end
        
        if isnan(hffr(i))
            hffr(i)=0;
        end
        
        if isnan(hffl(i))
            hffl(i)=0;
        end
        
        if isnan(pff(i))
            pff(i)=1;
        end
        
        if isnan(pffr(i))
            pffr(i)=1;
        end
        
        if isnan(pffl(i))
            pffl(i)=1;
        end
        
        xffax=[];
        for j=find(range,1,'first'):find(range,1,'last')
            xffax=[xffax;table1ff.("Unit Pairs"){j}];
        end        
        
        xffax=categorical(xffax);
        
        subplot(ax(i))
        hold on
        scatter(xffax,t1_ff_spikes)
        scatter(xffax,t2_ff_spikes)
        title({strcat('Feedforward: ', sub_conx(i)),strcat('Two tailed h: ',string(hff(i)),' p: ', string(pff(i))),strcat('Downreg h: ',string(hffr(i)),' p: ', string(pffr(i))),strcat('Upreg h: ',string(hffl(i)),' p: ', string(pffl(i)))})
        ylabel("Spikes")
        legend({'Pre-Stim','Post-Stim'})
        %dim=[i*.15 .7 .3 .15];
        %str={strcat('h:',string(hff)),strcat('p: ',string(pff)),strcat('ci: [',string(ciff(1)),',',string(ciff(2)),']'),strcat('tstat: ',string(statsff.tstat))};
        %annotation('textbox',dim,'String',str,'FitBoxToText','on')
        hold off
    
    else
        %placeholders
        hff(i)=-1;
        pff(i)=-1;
        ciff(i,:)=[0 0];
        statsff{i}=struct;
        
        hffr(i)=-1;
        pffr(i)=-1;
        ciffr(i,:)=[0 0];
        statsffr{i}=struct;
        
        hffl(i)=-1;
        pffl(i)=-1;
        ciffl(i,:)=[0 0];
        statsffl{i}=struct;
        subplot(ax(i))
        title({strcat('Feedforward: ', sub_conx(i)),strcat('Two tailed h: NaN',' p: NaN'),strcat('Downreg h: NaN',' p: NaN'),strcat('Upreg h: NaN',' p: NaN')})
    end
    
    t1_fb_spikes=[table1fb.("Spike Area"){range}];
    t2_fb_spikes=[table2fb.("Spike Area"){range}];
    
    if (~isempty(t1_fb_spikes)&&~isempty(t2_fb_spikes))&&length(t1_fb_spikes)>2

        [hfb(i),pfb(i),cifb(i,:),statsfb{i}]=ttest(t1_fb_spikes,t2_fb_spikes);
        [hfbr(i),pfbr(i),cifbr(i,:),statsfbr{i}]=ttest(t1_fb_spikes,t2_fb_spikes,'Tail','right');
        [hfbl(i),pfbl(i),cifbl(i,:),statsfbl{i}]=ttest(t1_fb_spikes,t2_fb_spikes,'Tail','left');

        % Convert NaN to []
        % Sorry for hard code


        if isnan(hfb(i))
            hfb(i)=0;
        end

        if isnan(hfbr(i))
            hfbr(i)=0;
        end

        if isnan(hfbl(i))
            hfbl(i)=0;
        end

        % p vals
        if isnan(pfb(i))
            pfb(i)=1;
        end

        if isnan(pfbr(i))
            pfbr(i)=1;
        end

        if isnan(pfbl(i))
            pfbl(i)=1;
        end

        xfbax=[];
        for j=find(range,1,'first'):find(range,1,'last')
            xfbax=[xfbax;table1fb.("Unit Pairs"){j}];
        end

        xfbax=categorical(xfbax);

        subplot(ax(i+5))
        hold on
        scatter(xfbax,t1_fb_spikes)
        scatter(xfbax,t2_fb_spikes)
        title({strcat('Feedback: ', sub_conx(i)),strcat('Two tailed h: ',string(hfb(i)),' p: ', string(pfb(i))),strcat('Downreg h: ',string(hfbr(i)),' p: ', string(pfbr(i))),strcat('Upreg h: ',string(hfbl(i)),' p: ', string(pfbl(i)))})
        ylabel("Spikes")
        legend({'Pre-Stim','Post-Stim'})
        hold off
    else
        %placeholders
        hfb(i)=-1;
        pfb(i)=-1;
        cifb(i,:)=[0 0];
        statsfb{i}=struct;
        
        hfbr(i)=-1;
        pfbr(i)=-1;
        cifbr(i,:)=[0 0];
        statsfbr{i}=struct;
        
        hfbl(i)=-1;
        pfbl(i)=-1;
        cifbl(i,:)=[0 0];
        statsfbl{i}=struct;
        
        subplot(ax(i+5))
        title({strcat('Feedforward: ', sub_conx(i)),strcat('Two tailed h: NaN',' p: NaN'),strcat('Downreg h: NaN',' p: NaN'),strcat('Upreg h: NaN',' p: NaN')})
    end
    
    
    sgtitle(strcat(table1ff.("Recording Name"){1}(1:34),' Paired T-Test'))
    %figure('units','normalized','outerposition',[0 0 1 1])
end

%saveas(gcf,fullfile('D:\Brewer lab data\Yash',strcat(table1ff.("Recording Name"){1}(1:34),' Paired T-Test')),'png')
saveas(gcf,fullfile('D:\Brewer lab data\ttest files',strcat(table1ff.("Recording Name"){1}(1:34),' Paired T-Test')),'png')

directionality_of_t=["Two Tail";"Right/Downreg";"Left/Upreg"];
h=[hff; hffr; hffl];
p=[pff; pffr; pffl];
ci=[{ciff};{ciffr};{ciffl}];
stats=[statsff;statsffr;statsffl];
fid=[table1ff.fid(1)';table1ff.fid(1)';table1ff.fid(1)'];
direction=["feedforward";"feedforward";"feedforward"];

forward_t_table=table(directionality_of_t,h,p,ci,stats,fid,direction);

h=[hfb; hfbr; hfbl];
p=[pfb; pfbr; pfbl];
ci=[{cifb};{cifbr};{cifbl}];
stats=[statsfb;statsfbr;statsfbl];
fid=[table1fb.fid(1)';table1fb.fid(1)';table1fb.fid(1)'];
direction=["feeback";"feeback";"feeback"];

backward_t_table=table(directionality_of_t,h,p,ci,stats,fid,direction);

hff=hff';
pff=pff';
hffr=hffr';
pffr=pffr';
hffl=hffl';
pffl=pffl';

hfb=hfb';
pfb=pfb';
hfbr=hfbr';
pfbr=pfbr';
hfbl=hfbl';
pfbl=pfbl';

fft2t=table(sub_conx,hff,pff,ciff,statsff');
fbt2t=table(sub_conx,hfb,pfb,cifb,statsfb');
fftrt=table(sub_conx,hffr,pffr,ciffr,statsffr');
fbtrt=table(sub_conx,hfbr,pfbr,cifbr,statsfbr');
fftlt=table(sub_conx,hffl,pffl,ciffl,statsffl');
fbtlt=table(sub_conx,hfbl,pfbl,cifbl,statsfbl');


% Save the plot and results
% Hard code directory here
%save(fullfile('D:\Brewer lab data\Yash',strcat('ttests ',table1ff.("Recording Name"){1}(1:34))),'fft2t','fbt2t','fftrt','fbtrt','fftlt','fbtlt')
save(fullfile('D:\Brewer lab data\ttest files',strcat('ttests ',table1ff.("Recording Name"){1}(1:34))),'fft2t','fbt2t','fftrt','fbtrt','fftlt','fbtlt')


end