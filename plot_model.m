n=5;
cntrsx=(1:100:100*n)';
cntrsy=14*ones(n,1);
strns=[1 12 25 37 50];

%%

i=2;
factor=4.5;
fig_prop(10,8);
for j=1:5
    hold on
    k=strns(1,j);
    pos=[cntrsx(j,1) cntrsy(j,1)-(X(k,i))/factor 2*(X(k,i))/factor 2*(X(k,i))/factor];
    rectangle('Position',pos,'Curvature',[1 1],'FaceColor',CT(j,:));
    axis equal;
    xlim([-10 450]);
    set(gca,'xtick',[]);
    text(cntrsx(j,1)+(X(k,i))/factor-5,cntrsy(j,1)+2*(X(k,i))/factor,num2str(j),'FontSize',20);
end

for j=1:4
    k=strns(1,j);
    l=k+1;
    fill([cntrsx(j,1)+2*(X(k,i))/factor,cntrsx(j,1)+2*(X(k,i))/factor,cntrsx(j+1,1),cntrsx(j+1,1)],[cntrsy(j,1)-A(k,l)/80,cntrsy(j,1)+A(k,l)/80,cntrsy(j,1)+A(l,k)/80,cntrsy(j,1)-A(l,k)/80],CT(6,:));
    axis equal;
end

pause(0.1)

%%

n=5;
cntrsx=(1:2000:2000*n)';
cntrsy=14*ones(n,1);
cntrs=[cntrsx,cntrsy];
strns=[1 12 25 37 50];

%%

factor=4.5;
factor1=4;
fig_prop(10,8);
for i=3:10:100
    clf;
    for j=1:5
        hold on
        k=strns(1,j);
        pos=[cntrsx(j,:) cntrsy(j,:)-(X(k,i))/factor 2*(X(k,i))/factor 2*(X(k,i))/factor];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor',CT(j,:));
        axis equal;
        xlim([-1000 14000]);
        ylim([-4000 4000]);
        set(gca,'xtick',[]);
    end
    
    for j=1:4
        hold on
        k=strns(1,j);
        l=k+1;
        fill([cntrsx(j,:)+2*(X(k,i))/factor,cntrsx(j,:)+2*(X(k,i))/factor,cntrsx(j+1,:),cntrsx(j+1,:)],[cntrsy(j,:)-A(k,l)/factor1,cntrsy(j,:)+A(k,l)/factor1,cntrsy(j,:)+A(l,k)/factor1,cntrsy(j,:)-A(l,k)/factor1],CT(6,:));
        axis equal;
        xlim([-1000 14000]);
        ylim([-4000 4000]);
    end
    
    pause(0.5)
end