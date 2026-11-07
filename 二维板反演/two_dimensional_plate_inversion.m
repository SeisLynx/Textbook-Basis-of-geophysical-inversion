function result = two_dimensional_plate_inversion(data_file, opts)
% 二维板磁异常乘法型阻尼最小二乘反演，无需 Optimization Toolbox。
% 直接运行 two_dimensional_plate_inversion，或传入实测CSV文件。
% p = [x0,z0,2b,2L,alpha_deg,i_deg,M]，z向下为正，距离m、角度deg。
% 详见同目录 README_二维板反演.md。
if nargin<1, data_file=''; end
if nargin<2, opts=struct(); end
defaults=struct('initial',[900 1100 500 700 50 60 2000], ...
    'active',1:5,'component','ppt','unit','ppt','noise',0.01, ...
    'background',0,'max_iter',200,'plot',true, ...
    'out',fullfile(fileparts(mfilename('fullpath')),'inversion_results_matlab'));
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(opts,names{k}), opts.(names{k})=defaults.(names{k}); end
end
assert(numel(opts.initial)==7 && all(isfinite(opts.initial)), 'initial必须含7个有限参数');
assert(~isempty(opts.active) && all(ismember(opts.active,1:7)) && numel(unique(opts.active))==numel(opts.active),'active应为不重复的1至7');
assert(ismember(opts.component,{'ppt','dx','dz','dt'}),'未知component');
assert(ismember(opts.unit,{'ppt','si'}),'未知unit');
assert(opts.noise>=0,'noise应非负');
truth=[];
if isempty(data_file)
    truth=[1000 1000 400 800 60 60 2000];
    x=(0:20:2000)'; z=zeros(size(x));
    clean=XASU(truth,x,z,opts);
    sd=opts.noise*max(abs(clean));
    sigma=ones(size(x))*max(sd,1e-10);
    rng(2026); obs=clean+opts.background+sd*randn(size(x));
else
    data=readmatrix(data_file,'NumHeaderLines',1);
    assert(ismember(size(data,2),[3 4]),'CSV必须有3或4列');
    x=data(:,1); z=data(:,2); obs=data(:,3);
    sigma=ones(size(x));
    if size(data,2)==4, sigma=data(:,4); end
end
assert(all(isfinite([x;z;obs;sigma])) && all(sigma>0),'输入必须有限，std必须为正');
assert(numel(x)>numel(opts.active),'测点数必须大于自由参数个数');
[p,pred,history,status,sv]=MAQT(x,z,obs,sigma,opts);
if ~exist(opts.out,'dir'), mkdir(opts.out); end
param_names={'x0_m','z0_m','width_2b_m','length_2L_m','alpha_deg','i_deg','M'};
writetable(array2table([x z obs sigma],'VariableNames',{'x_m','z_m','signal','std'}),fullfile(opts.out,'observations.csv'));
writetable(array2table([x z obs pred obs-pred],'VariableNames',{'x_m','z_m','observed','predicted','residual'}),fullfile(opts.out,'fit.csv'));
writetable(array2table(history,'VariableNames',[{'iteration','weighted_SSE','lambda'},param_names]),fullfile(opts.out,'history.csv'));
result=struct('parameters',p,'parameter_names',{param_names},'initial',opts.initial, ...
    'truth',truth,'active',opts.active,'component',opts.component,'unit',opts.unit, ...
    'status',status,'background',opts.background,'max_iter',opts.max_iter, ...
    'rmse',sqrt(mean((obs-pred).^2)), ...
    'scaled_jacobian_singular_values',sv);
fid=fopen(fullfile(opts.out,'result.json'),'w','n','UTF-8');
assert(fid>=0,'无法写入result.json');
fprintf(fid,'%s',jsonencode(result)); fclose(fid);
save(fullfile(opts.out,'result.mat'),'result','x','z','obs','pred','sigma','history','opts');
disp(result); disp(array2table(p,'VariableNames',param_names));
if opts.plot
    fig=figure('Color','w','Position',[80 80 1200 800]);
    subplot(2,2,1); plot(x,obs,'.',x,pred,'-',x,XASU(opts.initial,x,z,opts)+opts.background,'--');
    legend('Observed','Inverted','Initial'); xlabel('x (m)');
    if strcmp(opts.unit,'si'), ylabel('nT'); else, ylabel('Textbook amplitude'); end
    grid on;
    subplot(2,2,2); plot(x,obs-pred); xlabel('x (m)'); ylabel('Residual'); grid on;
    subplot(2,2,3); hold on;
    models=[opts.initial;p]; labels={'Initial','Inverted'};
    if ~isempty(truth), models=[models;truth]; labels{end+1}='True'; end
    model_bottom=max(z);
    for k=1:size(models,1)
        c=plate_corners(models(k,:)); c=c([1 2 4 3 1],:);
        model_bottom=max(model_bottom,max(c(:,2)));
        plot(c(:,1),c(:,2),'LineWidth',1.5,'DisplayName',labels{k});
    end
    plot(x,z,'k.','HandleVisibility','off'); axis equal; set(gca,'YDir','reverse');
    ylim([min(z)-100 model_bottom+100]);
    xlabel('x (m)'); ylabel('Depth z (m)'); legend('Location','best'); grid on;
    subplot(2,2,4); semilogy(history(:,1),max(history(:,2),1e-30),'o-');
    xlabel('Iteration'); ylabel('Weighted SSE'); grid on;
    axes_handles=findall(fig,'Type','axes');
    for ax_index=1:numel(axes_handles)
        axes_handles(ax_index).Toolbar.Visible='off';
    end
    exportgraphics(fig,fullfile(opts.out,'inversion.png'),'Resolution',180);
end
end

function c=plate_corners(p)
u=[cosd(p(5)) sind(p(5))];
t=p(1:2)-p(4)/2*u; b=p(1:2)+p(4)/2*u;
c=[t+[-p(3)/2 0]; t+[p(3)/2 0]; b+[-p(3)/2 0]; b+[p(3)/2 0]];
end

function ok=valid_model(p,z)
lo=[-2000 10 10 10 5 -89 1]; hi=[4000 5000 2500 4000 175 89 10000];
c=plate_corners(p);
ok=all(p>=lo & p<=hi) && min(c(:,2))>max(z)+1;
end

function out=XASU(p,x,z,opts)
% PPT第15页公式，atan2处理角度象限；log形式避免乘除溢出。
a0=deg2rad(p(5)); inc=deg2rad(p(6));
if strcmp(opts.component,'ppt')
    a=p(2)-z+p(4)*sin(a0)/2; b=p(2)-z-p(4)*sin(a0)/2;
    c=x-p(1)+p(4)*cos(a0)/2; d=x-p(1)-p(4)*cos(a0)/2;
    e=c+p(3)/2; f=c-p(3)/2; g=d+p(3)/2; h=d-p(3)/2;
    logu=log(b.^2+f.^2)+log(a.^2+g.^2)-log(b.^2+e.^2)-log(a.^2+h.^2);
    v=atan2(p(3)*b,b.^2+c.^2-p(3)^2/4);
    w=atan2(p(3)*a,a.^2+d.^2-p(3)^2/4);
    out=p(7)/(2*pi)*sin(a0)*(sin(a0-inc)*logu/2+cos(a0-inc)*(v-w));
else
    pts=plate_corners(p); r=zeros(numel(x),4); phi=r;
    for k=1:4
        dx=pts(k,1)-x; dz=pts(k,2)-z;
        r(:,k)=hypot(dx,dz); phi(:,k)=atan2(dz,dx);
    end
    q=phi(:,1)-phi(:,2)-phi(:,3)+phi(:,4);
    lr=log(r(:,2))+log(r(:,3))-log(r(:,1))-log(r(:,4));
    hx=p(7)/(2*pi)*sin(a0)*(lr*cos(a0-inc)-sin(a0-inc)*q);
    hz=p(7)/(2*pi)*sin(a0)*(lr*sin(a0-inc)+cos(a0-inc)*q);
    switch opts.component
        case 'dx', out=hx;
        case 'dz', out=hz;
        case 'dt', out=hx*cos(inc)+hz*sin(inc);
    end
end
if strcmp(opts.unit,'si'), out=out*(4*pi*1e2); end
end

function J=DETU(p,x,z,opts)
% 中心差分，步长随参数尺度变化。
J=zeros(numel(x),numel(opts.active));
for k=1:numel(opts.active)
    j=opts.active(k); h=1e-4*max(abs(p(j)),1);
    pp=p; pm=p; pp(j)=pp(j)+h; pm(j)=pm(j)-h;
    J(:,k)=(XASU(pp,x,z,opts)-XASU(pm,x,z,opts))/(2*h);
end
end

function [F,r]=QF(obs,pred,sigma,bo)
r=(obs-pred-bo)./sigma; F=r'*r;
end

function [p,pred,history,status,sv]=MAQT(x,z,obs,sigma,opts)
p=reshape(opts.initial,1,[]); assert(valid_model(p,z),'初始板体应位于全部测点下方且在边界内');
lambda=0.1; pred=XASU(p,x,z,opts); [F,r]=QF(obs,pred,sigma,opts.background);
history=[0 F lambda p]; status='max_iter';
for it=1:opts.max_iter
    J=DETU(p,x,z,opts)./sigma;
    scale=max(sqrt(sum(J.^2,1)),1e-15); Js=J./scale;
    numerical_floor=100*eps*max(1,max(abs(obs-opts.background)));
    if max(abs(Js'*r))<1e-9 || F<1e-20 || max(abs(obs-pred-opts.background))<=numerical_floor
        status='gradient_or_exact_fit'; break;
    end
    accepted=false;
    for retry=1:24
        A=Js'*Js+lambda*eye(numel(opts.active));
        R=chol(A); y=R\(R'\(Js'*r));
        trial=p; trial(opts.active)=trial(opts.active)+(y'./scale);
        if valid_model(trial,z)
            tp=XASU(trial,x,z,opts); [tf,tr]=QF(obs,tp,sigma,opts.background);
            if tf<F
                relative=(F-tf)/max(F,1e-30);
                p=trial; pred=tp; F=tf; r=tr;
                history(end+1,:)=[it F lambda p]; %#ok<AGROW>
                lambda=max(lambda/10,1e-12); accepted=true; break;
            end
        end
        lambda=lambda*10;
        if lambda>1e16, break; end
    end
    if ~accepted, status='no_acceptable_step'; break; end
    if relative<1e-9, status='relative_objective'; break; end
end
J=DETU(p,x,z,opts)./sigma;
sv=svd(J./max(sqrt(sum(J.^2,1)),1e-15));
pred=pred+opts.background;
end
