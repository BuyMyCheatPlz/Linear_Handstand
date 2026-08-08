function interactive_pendulum()
% 一阶倒立摆 LQR — 交互式仪表盘 (性能优化版)
% 操作: [推摆杆]/P | [推小车]/C | [重置]/R | [空格]暂停 | [↑↓]调强度

% ── 关闭旧实例 ───────────────────────────────────────
old = findobj('Type','figure','-and','Name','一阶倒立摆 LQR 交互式仿真');
for f = old'; set(f,'CloseRequestFcn','closereq'); close(f); end
clc;

% ═══════════════════════════════════════════════════════
%  参数 & LQR
% ═══════════════════════════════════════════════════════
M=1.0; m=0.1; ell=0.5; g=9.81; d=0.05;
A=[0,1; (M+m)*g/(M*ell), -d];
B=[0; -1/(M*ell)];
Q=diag([100,1]); R=0.01;
[K,~,cle]=lqr(A,B,Q,R);
Acl=A-B*K;

fprintf('══ 倒立摆 LQR 仪表盘 ══\nK=[%.2f, %.2f]  λ=%.2f,%.2f\n\n',K(1),K(2),cle);

% ═══════════════════════════════════════════════════════
%  轻量句柄 (只存 graphics handles + 常量, 大数组另存)
% ═══════════════════════════════════════════════════════
C = struct('bg',[0.11 0.11 0.13],'panel',[0.17 0.17 0.19],...
    'txt',[0.88 0.88 0.88],'grid',[0.25 0.25 0.28],...
    'blue',[0.2 0.6 1.0],'red',[1.0 0.3 0.3],...
    'green',[0.25 0.88 0.45],'orange',[1.0 0.63 0.1]);

% ═══════════════════════════════════════════════════════
%  动态状态 (主循环直接访问, 零 guidata 开销)
% ═══════════════════════════════════════════════════════
st_x   = [10*pi/180; 0];   st_t   = 0;
st_cx  = 0;                st_cv  = 0;
st_paused = false;

pend_dist  = 0;    cart_f = 0;    cart_tleft = 0;

dist_t=[];  dist_th=[];   push_t=[];  push_th=[];
flash_timer = 0;

% 降采样相平面缓冲区
phase_buf = nan(1500,2);  phase_idx=0;  phase_cnt=0;

% ═══════════════════════════════════════════════════════
%  GUI 搭建
% ═══════════════════════════════════════════════════════
FW=1550; FH=920;

fig = figure('Name','一阶倒立摆 LQR 交互式仿真',...
    'Position',[30 30 FW FH],'NumberTitle','off',...
    'Color',C.bg,'KeyPressFcn',@keyCb,'CloseRequestFcn',@closeCb);

% --- 动画轴 ---
ax_anim = axes('Units','pixels','Position',[40 FH-440 540 400]);
hold(ax_anim,'on');
set(ax_anim,'Color',[0.06 0.06 0.08],'XColor',C.txt,'YColor',C.txt,...
    'XLim',[-0.9 0.9],'YLim',[-0.12 0.72],'DataAspectRatio',[1 1 1],...
    'XTick',[],'YTick',[],'Box','on');
title(ax_anim,'倒立摆实时动画','Color',C.txt,'FontSize',13,'FontWeight','bold');

ground = plot(ax_anim,[-2 2],[0 0],'w-','LineWidth',2);
% 地面标尺: NaN 分隔单向量, 避免矩阵维度不匹配
N_TICK=21;
gx0=linspace(-1,1,N_TICK);
gt_x = [gx0; gx0; nan(1,N_TICK)];
gt_y = [zeros(1,N_TICK); 0.015*ones(1,N_TICK); nan(1,N_TICK)];
gt = line(ax_anim, gt_x(:), gt_y(:), 'Color','w','LineWidth',0.5);

cart_w=0.16; cart_h=0.06;
cart_r = rectangle(ax_anim,'Position',[-cart_w/2 0 cart_w cart_h],...
    'FaceColor',C.blue,'EdgeColor','none','Curvature',[0.1 0.1]);
whl_l = rectangle(ax_anim,'Position',[-0.05 -0.02 0.04 0.04],...
    'Curvature',[1 1],'FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
whl_r = rectangle(ax_anim,'Position',[0.01 -0.02 0.04 0.04],...
    'Curvature',[1 1],'FaceColor',[0.3 0.3 0.3],'EdgeColor','none');

pend_h = plot(ax_anim,[0 0],[cart_h ell+cart_h],'w-','LineWidth',3.5);
pv_r=0.028;
pivot = rectangle(ax_anim,'Position',[-pv_r cart_h-pv_r 2*pv_r 2*pv_r],...
    'Curvature',[1 1],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
ms_r=0.048;
mass = rectangle(ax_anim,'Position',[-ms_r ell+cart_h-ms_r 2*ms_r 2*ms_r],...
    'Curvature',[1 1],'FaceColor',C.red,'EdgeColor','none');
plot(ax_anim,[0 0],[0 ell+0.15],':','Color',[0.3 0.3 0.3],'LineWidth',0.8);

txt_th = text(ax_anim,-0.85,0.64,'','Color',C.txt,'FontSize',11,'FontName','Consolas');
txt_td = text(ax_anim,-0.85,0.56,'','Color',C.txt,'FontSize',11,'FontName','Consolas');
txt_cx = text(ax_anim,-0.85,0.48,'','Color',C.txt,'FontSize',11,'FontName','Consolas');
txt_f  = text(ax_anim,-0.85,0.40,'','Color',C.txt,'FontSize',11,'FontName','Consolas');
txt_flash = text(ax_anim,0,0.35,'','Color',C.orange,'FontSize',28,...
    'FontWeight','bold','HorizontalAlignment','center');

% --- 相平面 ---
ax_phase = axes('Units','pixels','Position',[600 FH-440 360 400]);
hold(ax_phase,'on');
set(ax_phase,'Color',C.panel,'XColor',C.txt,'YColor',C.txt,...
    'GridColor',C.grid,'GridAlpha',0.5);
grid(ax_phase,'on');
xlabel(ax_phase,'\theta (°)','Color',C.txt);
ylabel(ax_phase,'\theta'' (°/s)','Color',C.txt);
title(ax_phase,'相平面轨迹','Color',C.txt,'FontSize',12);
set(ax_phase,'XLim',[-20 20],'YLim',[-150 150]);
plot(ax_phase,[0 0],[-300 300],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
plot(ax_phase,[-300 300],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
phase_traj = plot(ax_phase,nan,nan,'-','Color',C.blue,'LineWidth',1.5);
phase_now  = plot(ax_phase,nan,nan,'o','Color',C.red,'MarkerSize',10,'MarkerFaceColor',C.red);
phase_mark = scatter(ax_phase,nan,nan,50,C.orange,'filled','v','MarkerEdgeColor','none');
plot(ax_phase,10,0,'go','MarkerSize',8,'LineWidth',1.5);

% --- 信息面板 ---
info_pos = [980 FH-440 530 400];
ax_info = axes('Units','pixels','Position',info_pos);
set(ax_info,'Color','none','XColor','none','YColor','none','XLim',[0 1],'YLim',[0 1]);
hold(ax_info,'on');
info_lines = {sprintf('小车 M=%.1f  摆杆 m=%.2f  l=%.2f',M,m,ell),...
    sprintf('阻尼 d=%.3f  重力 g=%.2f',d,g),'',...
    sprintf('LQR: Q=diag(100,1)  R=%.2f',R),...
    sprintf('K=[%.2f, %.2f]  λ=%.2f,%.2f',K(1),K(2),cle),'',...
    '快捷键: P=推摆杆 C=推小车 R=重置',...
    '        空格=暂停 ↑↓=调强度'};
for i=1:length(info_lines)
    text(ax_info,0.05,0.95-(i-1)*0.06,info_lines{i},...
        'Color',C.txt,'FontSize',10,'FontName','Consolas','VerticalAlignment','top');
end

% --- 按钮 ---
bw=105; bh=38;
uicontrol('Style','pushbutton','String','推摆杆','FontSize',12,'FontWeight','bold',...
    'Position',[info_pos(1)+0.02*info_pos(3), info_pos(2)+0.32*info_pos(4), bw, bh],...
    'BackgroundColor',C.orange*0.7,'ForegroundColor',[0 0 0],'Callback',@pendCb);
uicontrol('Style','pushbutton','String','推小车','FontSize',12,'FontWeight','bold',...
    'Position',[info_pos(1)+0.02*info_pos(3)+bw+12, info_pos(2)+0.32*info_pos(4), bw, bh],...
    'BackgroundColor',C.green*0.7,'ForegroundColor',[0 0 0],'Callback',@cartCb);
uicontrol('Style','pushbutton','String','重置','FontSize',12,'FontWeight','bold',...
    'Position',[info_pos(1)+0.02*info_pos(3)+2*(bw+12), info_pos(2)+0.32*info_pos(4), bw, bh],...
    'BackgroundColor',[0.45 0.45 0.45],'ForegroundColor',[1 1 1],'Callback',@resetCb);
btn_pause = uicontrol('Style','pushbutton','String','暂停','FontSize',12,'FontWeight','bold',...
    'Position',[info_pos(1)+0.02*info_pos(3)+3*(bw+12), info_pos(2)+0.32*info_pos(4), bw, bh],...
    'BackgroundColor',[0.35 0.35 0.45],'ForegroundColor',[1 1 1],'Callback',@pauseCb);

slider_h = uicontrol('Style','slider','Min',0.1,'Max',5.0,'Value',1.0,...
    'Position',[info_pos(1)+72, info_pos(2)+0.32*info_pos(4)-44, 320, 22],...
    'BackgroundColor',[0.25 0.25 0.30],'Callback',@sliderCb);
uicontrol('Style','text','String','扰动强度','FontSize',10,...
    'Position',[info_pos(1)+5, info_pos(2)+0.32*info_pos(4)-48, 65, 20],...
    'BackgroundColor',C.bg,'ForegroundColor',C.txt);
slider_lbl = uicontrol('Style','text','String','1.0x','FontSize',10,...
    'Position',[info_pos(1)+398, info_pos(2)+0.32*info_pos(4)-44, 50, 20],...
    'BackgroundColor',C.bg,'ForegroundColor',C.txt);

% --- 时域图: animatedline (增量更新, 性能核心优化) ---
row_h=130; row_gap=10; row_w=FW-80;
ry1=FH-440-row_h-row_gap-10; ry2=ry1-row_h-row_gap; ry3=ry2-row_h-row_gap;
max_pts = 2000;

% θ(t)
ax_th = axes('Units','pixels','Position',[40 ry1 row_w row_h]); hold(ax_th,'on');
set(ax_th,'Color',C.panel,'XColor',C.txt,'YColor',C.txt,'GridColor',C.grid,'GridAlpha',0.4);
grid(ax_th,'on'); ylabel(ax_th,'\theta (°)','Color',C.txt);
title(ax_th,'摆杆角度  \theta(t)','Color',C.txt,'FontSize',12);
set(ax_th,'XLim',[0 15],'YLim',[-28 28]);
yline(ax_th,0,'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
al_th = animatedline(ax_th,'Color',C.blue,'LineWidth',1.8,'MaximumNumPoints',max_pts);
pend_marks = scatter(ax_th,nan,nan,55,C.orange,'filled','v','MarkerEdgeColor','none');
cart_marks = scatter(ax_th,nan,nan,65,C.green,'filled','^','MarkerEdgeColor','none');

% θ̇(t)
ax_td = axes('Units','pixels','Position',[40 ry2 row_w row_h]); hold(ax_td,'on');
set(ax_td,'Color',C.panel,'XColor',C.txt,'YColor',C.txt,'GridColor',C.grid,'GridAlpha',0.4);
grid(ax_td,'on'); ylabel(ax_td,'\theta'' (°/s)','Color',C.txt);
title(ax_td,'角速度  θ''(t)','Color',C.txt,'FontSize',12);
set(ax_td,'XLim',[0 15],'YLim',[-200 200]);
yline(ax_td,0,'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
al_td = animatedline(ax_td,'Color',C.red,'LineWidth',1.8,'MaximumNumPoints',max_pts);

% F(t)
ax_f = axes('Units','pixels','Position',[40 ry3 row_w row_h]); hold(ax_f,'on');
set(ax_f,'Color',C.panel,'XColor',C.txt,'YColor',C.txt,'GridColor',C.grid,'GridAlpha',0.4);
grid(ax_f,'on'); xlabel(ax_f,'时间 (s)','Color',C.txt); ylabel(ax_f,'力 F (N)','Color',C.txt);
title(ax_f,'控制力  F(t)','Color',C.txt,'FontSize',12);
set(ax_f,'XLim',[0 15],'YLim',[-60 60]);
yline(ax_f,0,'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
al_f = animatedline(ax_f,'Color',C.green,'LineWidth',1.8,'MaximumNumPoints',max_pts);

% ── 把回调需要的 handles 打包存入 guidata ────────────
H = struct();
H.fig=fig; H.slider=slider_h; H.slider_lbl=slider_lbl;
H.btn_pause=btn_pause; H.pulse_dur=0.10;
H.ax_th=ax_th; H.ax_td=ax_td; H.ax_f=ax_f;
H.ax_phase=ax_phase;
H.phase_traj=phase_traj; H.phase_now=phase_now; H.phase_dots=phase_mark;
H.pend_marks=pend_marks; H.cart_marks=cart_marks;
H.al_th=al_th; H.al_td=al_td; H.al_f=al_f;
H.anim=struct('ax',ax_anim,'ground',ground,'gt',gt,...
    'cart_r',cart_r,'whl_l',whl_l,'whl_r',whl_r,...
    'pend_h',pend_h,'pivot',pivot,'mass',mass,...
    'cart_w',cart_w,'cart_h',cart_h,'pv_r',pv_r,'ms_r',ms_r,...
    'txt_th',txt_th,'txt_td',txt_td,'txt_cx',txt_cx,'txt_f',txt_f,...
    'txt_flash',txt_flash);
H.C=C; H.M=M; H.m=m; H.ell=ell; H.K=K; H.Acl=Acl; H.B=B;
guidata(fig, H);

% ── UserData 用于主循环↔回调通信 ─────────────────────
ud = struct('pend_dist',0,'cart_f',0,'cart_tleft',0,...
    'do_pause',false,'do_reset',false,'do_quit',false);
set(fig,'UserData',ud);

fprintf('操作: [推摆杆]/P | [推小车]/C | [重置]/R | [空格]暂停 | [↑↓]强度\n\n');

% ═══════════════════════════════════════════════════════
%  主循环 — 直接访问工作区变量, 零 guidata 序列化开销
% ═══════════════════════════════════════════════════════
dt = 0.005;
anim_skip = 6;    % 每 6 帧 (~33 fps)
chart_skip = 12;  % 每 12 帧 (~16 fps)

fc=0; tic;
while isvalid(fig)
    % ── 处理回调信号 ────────────────────────────────
    ud = get(fig,'UserData');
    if ud.do_quit, break; end
    if ud.do_pause
        st_paused = ~st_paused;
        if st_paused, set(btn_pause,'String','继续');
        else, set(btn_pause,'String','暂停'); tic; end   % 重置计时防追帧
        ud.do_pause = false;
    end
    if ud.do_reset
        st_x=[10*pi/180;0]; st_t=0; st_cx=0; st_cv=0;
        pend_dist=0; cart_f=0; cart_tleft=0;
        dist_t=[]; dist_th=[]; push_t=[]; push_th=[];
        phase_idx=0; phase_cnt=0; phase_buf(:)=nan;
        clearpoints(al_th); clearpoints(al_td); clearpoints(al_f);
        set(pend_marks,'XData',nan,'YData',nan);
        set(cart_marks,'XData',nan,'YData',nan);
        set(phase_traj,'XData',nan,'YData',nan);
        set(phase_now,'XData',nan,'YData',nan);
        set(phase_mark,'XData',nan,'YData',nan);
        set([ax_th,ax_td,ax_f],'XLim',[0 15]);
        fprintf('  🔄 已重置\n');
        ud.do_reset = false;
    end
    % 扰动信号
    if ud.pend_dist ~= 0
        pend_dist = ud.pend_dist;
        ud.pend_dist = 0;
    end
    if ud.cart_f ~= 0
        cart_f = ud.cart_f; cart_tleft = ud.cart_tleft;
        ud.cart_f = 0; ud.cart_tleft = 0;
    end
    set(fig,'UserData',ud);

    if st_paused
        pause(0.05); tic; continue;
    end

    % ── 施加扰动 ────────────────────────────────────
    if pend_dist ~= 0
        st_x(2) = st_x(2) + pend_dist;
        dist_t(end+1)=st_t;  dist_th(end+1)=st_x(1)*180/pi;
        flash_timer=0.4;
        fprintf('  ⚡ 推摆杆! %.0f°/s\n', pend_dist*180/pi);
        pend_dist=0;
    end

    F_dist=0;
    if cart_tleft > 0
        F_dist=cart_f;
        if cart_tleft >= H.pulse_dur - dt*1.5
            push_t(end+1)=st_t;  push_th(end+1)=st_x(1)*180/pi;
            flash_timer=0.4;
        end
        cart_tleft=cart_tleft-dt;
        if cart_tleft<=0, fprintf('  🚗 推小车! F=%.1fN %.0fms\n',cart_f,H.pulse_dur*1000); end
    end

    % ── RK4 积分 ────────────────────────────────────
    x=st_x;  sdt=dt/4;
    for ss=1:4
        dx=Acl*x+B*F_dist;
        k1=dx; k2=Acl*(x+0.5*sdt*k1)+B*F_dist;
        k3=Acl*(x+0.5*sdt*k2)+B*F_dist; k4=Acl*(x+sdt*k3)+B*F_dist;
        x=x+(sdt/6)*(k1+2*k2+2*k3+k4);
    end
    theta_ddot=Acl(2,:)*st_x+B(2)*F_dist;
    st_x=x;  st_t=st_t+dt;

    % 小车位移 (显示用)
    F_ctrl=-K*x;
    cart_acc=(F_ctrl+F_dist-m*ell*theta_ddot)/(M+m);
    st_cv=st_cv+cart_acc*dt;  st_cx=st_cx+st_cv*dt;

    % ── 动画 (~33fps) ───────────────────────────────
    fc=fc+1;
    do_anim = mod(fc,anim_skip)==0;
    do_chart = mod(fc,chart_skip)==0;

    if do_anim
        a=H.anim; theta=st_x(1); cx=st_cx; cy=a.cart_h;
        tx=cx+ell*sin(theta); ty=cy+ell*cos(theta);

        gxl=cx+[-1.2 1.2];
        gx=linspace(cx-1,cx+1,N_TICK);
        set(a.ground,'XData',gxl);
        gt_xu=[gx; gx; nan(1,N_TICK)];
        set(a.gt,'XData',gt_xu(:));
        set(a.ax,'XLim',cx+[-0.9 0.9]);
        set(a.pend_h,'XData',[cx tx],'YData',[cy ty]);
        set(a.mass,'Position',[tx-H.anim.ms_r,ty-H.anim.ms_r,2*H.anim.ms_r,2*H.anim.ms_r]);
        set(a.cart_r,'Position',[cx-a.cart_w/2,0,a.cart_w,a.cart_h]);
        set(a.whl_l,'Position',[cx-0.05,-0.02,0.04,0.04]);
        set(a.whl_r,'Position',[cx+0.01,-0.02,0.04,0.04]);
        set(a.pivot,'Position',[cx-a.pv_r,cy-a.pv_r,2*a.pv_r,2*a.pv_r]);

        an=min(abs(theta)/(30*pi/180),1.0);
        set(a.pend_h,'Color',[an,1.0-an*0.7,1.0-an]);

        set(a.txt_th,'String',sprintf('θ  = %+7.2f°',theta*180/pi));
        set(a.txt_td,'String',sprintf('θ̇  = %+7.2f °/s',st_x(2)*180/pi));
        set(a.txt_cx,'String',sprintf('小车位移 = %+.3f m',cx));
        set(a.txt_f,'String',sprintf('控制力 = %+7.2f N',F_ctrl));

        af=abs(F_ctrl);
        if af>25, fcc=C.red; elseif af>8, fcc=C.orange; else, fcc=C.green; end
        set(a.txt_f,'Color',fcc);

        if flash_timer>0
            set(a.txt_flash,'String','⚡');
            flash_timer=flash_timer-dt*anim_skip;
        elseif flash_timer<=-0.05
            set(a.txt_flash,'String',''); flash_timer=0;
        else
            flash_timer=flash_timer-dt*anim_skip;
        end
    end

    % ── 图表增量更新 (~16fps) ───────────────────────
    if do_chart
        thd=st_x(1)*180/pi; tdd=st_x(2)*180/pi; ftot=F_ctrl+F_dist;

        addpoints(al_th,st_t,thd);
        addpoints(al_td,st_t,tdd);
        addpoints(al_f,st_t,ftot);

        % X 轴滚动
        if st_t>ax_th.XLim(2)
            nx=[0,max(15,st_t+2)];
            set([ax_th,ax_td,ax_f],'XLim',nx);
        end

        % Y 轴 (每 3s 检查一次)
        if mod(fc,500)==0
            [~,d]=getpoints(al_th); d=d(~isnan(d));
            if ~isempty(d), mx=max(abs(d)); if mx>ax_th.YLim(2)*0.8, ny=max(28,mx*1.3); set(ax_th,'YLim',[-ny ny]); end; end
            [~,d]=getpoints(al_td); d=d(~isnan(d));
            if ~isempty(d), mx=max(abs(d)); if mx>ax_td.YLim(2)*0.8, ny=max(200,mx*1.3); set(ax_td,'YLim',[-ny ny]); end; end
            [~,d]=getpoints(al_f); d=d(~isnan(d));
            if ~isempty(d), mx=max(abs(d)); if mx>ax_f.YLim(2)*0.8, ny=max(60,mx*1.3); set(ax_f,'YLim',[-ny ny]); end; end
        end

        % 扰动标记
        if ~isempty(dist_t)
            n=min(length(dist_t),25);
            set(pend_marks,'XData',dist_t(end-n+1:end),'YData',dist_th(end-n+1:end));
        end
        if ~isempty(push_t)
            n=min(length(push_t),25);
            set(cart_marks,'XData',push_t(end-n+1:end),'YData',push_th(end-n+1:end));
        end

        % 相平面
        phase_idx=mod(phase_idx,1500)+1;
        phase_buf(phase_idx,:)=[thd,tdd];
        phase_cnt=min(phase_cnt+1,1500);
        if phase_cnt<1500, idx=1:phase_idx;
        else, idx=[phase_idx+1:1500,1:phase_idx]; end
        set(phase_traj,'XData',phase_buf(idx,1),'YData',phase_buf(idx,2));
        set(phase_now,'XData',thd,'YData',tdd);

        if mod(fc,500)==0
            mx=max(abs(phase_buf(idx,1))); my=max(abs(phase_buf(idx,2)));
            if mx>ax_phase.XLim(2)*0.8, nx=max(20,mx*1.3); set(ax_phase,'XLim',[-nx nx]); end
            if my>ax_phase.YLim(2)*0.8, ny=max(150,my*1.3); set(ax_phase,'YLim',[-ny ny]); end
        end

        if ~isempty(dist_t)||~isempty(push_t)
            all_th=[dist_th,push_th]; n=min(length(all_th),25);
            set(phase_mark,'XData',all_th(end-n+1:end),'YData',zeros(1,n));
        end
    end

    % ── 渲染 ────────────────────────────────────────
    if do_anim||do_chart, drawnow limitrate; end

    % ── 时间同步 ────────────────────────────────────
    elapsed=toc; wait_t=st_t-elapsed;
    if wait_t>0.03, pause(wait_t*0.8);
    elseif wait_t<-1.5, tic; end
end
fprintf('仿真结束。\n');
end

%% ═══════════════════════════════════════════════════════
%  回调 (通过 UserData 发信号给主循环, 轻量)
% ═══════════════════════════════════════════════════════

function pendCb(~,~)
    H=guidata(gcbf); if isempty(H), return; end
    ud=get(H.fig,'UserData');
    mag=get(H.slider,'Value');
    ud.pend_dist = sign(randn())*mag*2.0;
    set(H.fig,'UserData',ud);
end

function cartCb(~,~)
    H=guidata(gcbf); if isempty(H), return; end
    ud=get(H.fig,'UserData');
    mag=get(H.slider,'Value');
    ud.cart_f=sign(randn())*mag*15.0;
    ud.cart_tleft=H.pulse_dur;
    set(H.fig,'UserData',ud);
end

function resetCb(~,~)
    H=guidata(gcbf); if isempty(H), return; end
    ud=get(H.fig,'UserData'); ud.do_reset=true;
    set(H.fig,'UserData',ud);
end

function pauseCb(~,~)
    H=guidata(gcbf); if isempty(H), return; end
    ud=get(H.fig,'UserData'); ud.do_pause=true;
    set(H.fig,'UserData',ud);
end

function sliderCb(~,~)
    H=guidata(gcbf); if isempty(H), return; end
    set(H.slider_lbl,'String',sprintf('%.1fx',get(H.slider,'Value')));
end

function keyCb(~,event)
    H=guidata(gcbf); if isempty(H), return; end
    ud=get(H.fig,'UserData');
    switch event.Key
        case 'space', ud.do_pause=true;
        case 'r', ud.do_reset=true;
        case 'p'
            mag=get(H.slider,'Value');
            ud.pend_dist=sign(randn())*mag*2.0;
        case 'c'
            mag=get(H.slider,'Value');
            ud.cart_f=sign(randn())*mag*15.0;
            ud.cart_tleft=H.pulse_dur;
        case 'uparrow'
            set(H.slider,'Value',min(5.0,get(H.slider,'Value')*1.5));
            set(H.slider_lbl,'String',sprintf('%.1fx',get(H.slider,'Value')));
        case 'downarrow'
            set(H.slider,'Value',max(0.1,get(H.slider,'Value')/1.5));
            set(H.slider_lbl,'String',sprintf('%.1fx',get(H.slider,'Value')));
    end
    set(H.fig,'UserData',ud);
end

function closeCb(~,~)
    H=guidata(gcbf);
    if ~isempty(H)&&isvalid(H.fig)
        ud=get(H.fig,'UserData'); ud.do_quit=true;
        set(H.fig,'UserData',ud); pause(0.05);
    end
    closereq();
end
