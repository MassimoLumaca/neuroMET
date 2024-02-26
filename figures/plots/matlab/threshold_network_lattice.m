% Script adapted by the author from CONNSPM CONN toolbox (Whitfield-Gabrieli (http://www.nitrc.org/projects/conn)
% Clear workspace, command window, and close any open figures
clear,clc, close all

connectome={'structural'}; % 'structural' or 'functional'
%%
cd('/Users/au540169/Documents/MATLAB/CPM')


% Use 'matrix' to create a matrix ZZ with NaN values in the diagonal

if strcmp(connectome,'structural')
    load('s_conn_vviq.mat')
    Z = matrix;
elseif strcmp(connectome,'functional')
    load('f_conn_vviq.mat')
else
    error('Invalid connectome. Use either ''structural'' or ''functional''.');

end


ZZ=zeros(size(Z));

% Rank-based normalization
for n=1:size(Z,3)
    z=Z(:,:,n);
    z(1:size(z,1)+1:end)=nan; 
    idx=find(~isnan(z));
    [st,sidx]=sort(z(idx)); z(idx(sidx))=((1:numel(idx))-.5)'/numel(idx);
    ZZ(:,:,n)=z;
end
N=size(ZZ,1);

%explores variable thresholding
THR=0:.01:.50;

% Initialize matrices to store network measures
net_d=zeros(size(ZZ,3),numel(THR));net_c=net_d;net_k=net_d;
net_d_rand=zeros(1,numel(THR));net_c_rand=net_d_rand;net_k_rand=net_d_rand;
net_d_latt=zeros(1,numel(THR));net_c_latt=net_d_latt;net_k_latt=net_d_latt;

h = waitbar(0,'Compute network measures based on thresholds...');
total_iterations = numel(THR);

% Compute network measures based on thresholds
for nthr=1:numel(THR)
    for n=1:size(ZZ,3)
        thr=THR(nthr);
        C=(ZZ(:,:,n)>=1-thr);
        [~,~,~,d,c,k]=conn_network_efficiency(C);
        net_c(n,nthr)=c;
        net_d(n,nthr)=d;
        net_k(n,nthr)=k;
    end
    thr=THR(nthr);
    [~,~,~,d,c,k]=conn_network_efficiency([N,mean(net_k(:,nthr))],2e1);
    net_c_rand(nthr)=c;
    net_d_rand(nthr)=d;
    net_k_rand(nthr)=k;
    [~,~,~,d,c,k]=conn_network_efficiency([N,mean(net_k(:,nthr))],-2e1);
    net_c_latt(nthr)=c;
    net_d_latt(nthr)=d;
    net_k_latt(nthr)=k;
    % Update the waitbar
    waitbar(nthr/total_iterations, h, sprintf('Processing Threshold %d/%d', nthr, total_iterations));
end

% Close the waitbar
close(h);

%plots network measures
x={net_d,net_c};%y={net_k};
x0={{net_d_rand,net_d_latt},{net_c_rand,net_c_latt}};%y0={net_k_rand};
y={net_k}; y0={net_k_rand};ylabels={{'\fontsize{16}Cost (K)'}};

xlabels={{'\fontsize{16}\color{black}Global efficiency'},{'\fontsize{16}\color{black}Local efficiency'}};
% Create a figure for plotting data
figure('name','Network theory: explore variable threshold range','numbertitle','off','units','norm','position',[.3 .3 .6 .6]);clf;

for n1=1:numel(x)
    thr=mean(y{1},1);
    mnet_d=mean(x{n1},1);
    snet_d=std(x{n1},0,1);
    n=size(x{n1},1);
    m=size(x{n1},2);
    alpha=.05;
    subplot(numel(x),2,2*n1-1);
    hold on;
    h=patch([thr,fliplr(thr)],[mnet_d+2*snet_d,fliplr(mnet_d-2*snet_d)],'k');
    %h=patch([thr,fliplr(thr)],[mnet_d+spm_invTcdf(1-alpha/m/2,n-1)*snet_d/sqrt(n),fliplr(mnet_d-spm_invTcdf(1-alpha/m/2,n-1)*snet_d/sqrt(n))],'k');
    set(h,'facecolor',.75*[1,1,1],'edgecolor','none');
    h1=plot(thr,mnet_d,'k-','linewidth',2);
    h2=[];linetypes={'--',':','-.','.-'}; linecolors={.5*[1,1,1],.5*[1,1,1],.5*[1,1,1],.5*[1,1,1]}; 
    for n2=1:numel(x0{n1})
        thr0=mean(y0{1},1);
        mnet_d0=mean(x0{n1}{n2},1);
        h2(n2)=plot(thr0,mnet_d0,[linetypes{n2}],'linewidth',2,'color',linecolors{n2});
    end
    hold off;
    set(gcf,'color','w');set(gca,'fontsize',16,'xcolor',.75*[1,1,1],'ycolor',.75*[1,1,1]);
    ylabel(xlabels{n1});xlabel(ylabels{1});
    if n1==1, legend([h1,h2],{'Data','Random graph','Lattice'}); end
end

% Create a subplot for the difference with lattice and random
h=subplot(1,2,2); set(h,'position',get(h,'position')*[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 .25 0 .5]);
mnet_d_1=mnet_d-mean(x0{1}{2},1); %title('Difference with lattice')
mnet_d_2=mnet_d-mean(x0{n1}{1},1); %title('Difference with random');
mnet_d=mnet_d_1+mnet_d_2; title('GE_{data}-GE_{lattice} + LE_{data}-LE_{random}');
snet_d_1=std(x{1},0,1);
snet_d_2=std(x{2},0,1);
snet_d=sqrt(snet_d_1.^2+snet_d_2.^2);
hold on;

h=patch([thr,fliplr(thr)],[mnet_d+2*snet_d,fliplr(mnet_d-2*snet_d)],'k');
set(h,'facecolor',.75*[1,1,1],'edgecolor','none');
h1=plot(thr,mnet_d,'k-','linewidth',2);
hold off;
set(gcf,'color','w');set(gca,'fontsize',16,'xcolor',.75*[1,1,1],'ycolor',.75*[1,1,1]);

xlabel(ylabels{1});

