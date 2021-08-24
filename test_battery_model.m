clear
close all

%% Variables
syms q z s FCC0 SOH taud Rd R0 omega
vars    = {     R0,       Rd, taud, SOH};
varvals = {0.565e-3, 0.896e-3,  224, 0.91};

%% Electrochemistry model 
G=(4.1-3.6)/(0.9-0.1)/(FCC0*SOH)/s+R0+Rd/sqrt(taud*s)*tanh(sqrt(taud*s));
G=subs(G,FCC0,65.6*3600);
G0=subs(G,vars,varvals);

%% 感度関数のプロット（電気化学モデル）
h=figure(1001);
clf
set(h,'Name','Sensitivity function (electrochemistry model)')
freq=logspace(-5,-1,100);
sensitibity_plot(G,freq,vars,varvals);
ylim([-100,-20])


%% 有理関数モデルの作成
% 有限次元近似モデルの作成
n=3;
Z=@(n) Rd/(4*n-3);
Y=@(n) taud*s/(Rd*(4*n-1));
Zw= @(m) Zw_rec(Z,Y,1,m+1);
Glow=Zw(n-1);
Glow=subs(Glow,FCC0,65.6*3600);
Glow0=subs(Glow,vars,varvals);
[Glow0n,Glow0d]=numden(Glow0);
varvals_rd=sym2poly(Glow0d);
nd=length(varvals_rd);
varvals_rn=sym2poly(Glow0n);
nn=length(varvals_rn);

% 有理関数モデルに変換
varvals_r=[varvals{1},varvals{4},varvals_rd,varvals_rn];
vars_rd=fliplr([sym('a0'),sym('a',[1,nd-1])]);
vars_rn=fliplr([sym('b0'),sym('b',[1,nn-1])]);
vars_r=[R0,SOH,vars_rd,vars_rn];
varvals_r=num2cell(varvals_r);
vars_r=num2cell(vars_r);

Gr=(4.1-3.6)/(0.9-0.1)/(FCC0*SOH)/s+R0+(s.^(nn-1:-1:0)*transpose(vars_rn))/(s.^(nd-1:-1:0)*transpose(vars_rd));
Gr=subs(Gr,FCC0,65.6*3600);
Gr0=subs(Gr,vars_r,varvals_r);

%% 離散時間有限次元近似モデルの作成
%for dt=logspace(3,-2,20)
dt=1;
Glow0dt=subs(Glow0,s,(2/dt)*(z-1)/(z+1));
Glow0dt=subs(Glow0dt,z,q^(-1));
[Glow0dtn,Glow0dtd]=numden(Glow0dt);
varvals_rdtd=sym2poly(Glow0dtd);
dtnd=length(varvals_rdtd);
varvals_rdtn=sym2poly(Glow0dtn);
dtnn=length(varvals_rdtn);

% 有理関数モデルに変換(離散時間)
varvals_rdt=[varvals{1},varvals{4},varvals_rdtd,varvals_rdtn];
vars_rdtd=fliplr([sym('da0'),sym('da',[1,dtnd-1])]);
vars_rdtn=fliplr([sym('db0'),sym('db',[1,dtnn-1])]);
vars_rdt=[R0,SOH,vars_rdtd,vars_rdtn];
varvals_rdt=num2cell(varvals_rdt);
vars_rdt=num2cell(vars_rdt);

Grdt=(4.1-3.6)/(0.9-0.1)/(FCC0*SOH)*(dt/(1-z^(-1)))+R0+(q.^(dtnn-1:-1:0)*transpose(vars_rdtn))/(q.^(dtnd-1:-1:0)*transpose(vars_rdtd));
Grdt=subs(Grdt,q,z^(-1));
Grdt=subs(Grdt,FCC0,65.6*3600);
Gr0dt=subs(Grdt,vars_rdt,varvals_rdt);





%% 感度関数のプロット（有限次元・有理関数モデル）
h=figure(1002);
clf
set(h,'Name','Sensitivity function (blackbox model)')
freq=logspace(-5,-1,100);
sensitibity_plot(Gr,freq,vars_r,varvals_r);
ylim([-100,-20])
print('-depsc2','output/blackbox.eps')
print('-dpng','output/blackbox.png')


%% カウエル型ブラックボックス
n=3;
Zw= Cauer(n);
Gcau = (4.1-3.6)/(0.9-0.1)/(FCC0*SOH)/s+R0+Zw;
Gcau = subs(Gcau,FCC0,65.6*3600);
h=figure(1003);
clf
set(h,'Name','Sensitivity function (Cauer black)')
freq=logspace(-5,-1,100);
Rd_v= 0.896e-3;
taud_v = 224;
Cd_v = taud_v/Rd_v;
vars_cau    = {     R0,  SOH};
varvals_cau = {0.565e-3,  0.91};
for k=1:n
    vars_cau{end+1} = sym(sprintf('R%d',k));
    varvals_cau{end+1} = Rd_v/(4*k-3);
    vars_cau{end+1} = sym(sprintf('C%d',k));
    varvals_cau{end+1} = Cd_v/(4*k-1);
end

sensitibity_plot(Gcau,freq,vars_cau,varvals_cau);
ylim([-100,-20])
print('-depsc2','output/cau.eps')
print('-dpng','output/cau.png')
    

%%　フォスター型ブラックボックス
n=3;
Gfos = (4.1-3.6)/(0.9-0.1)/(FCC0*SOH)/s+R0;
vars_fos    = {     R0,  SOH};
varvals_fos = {0.565e-3,  0.91};
Rd_v= 0.896e-3;
taud_v = 224;
Cd_v = taud_v/Rd_v;

for k=1:n
    Rk = sym(sprintf('R%d',k));
    Ck = sym(sprintf('C%d',k));
    Gfos = Gfos + Rk/(s*Rk*Ck+1);
    vars_fos{end+1} = Rk;
    varvals_fos{end+1} = 8*Rd_v/((2*k-1)^2*pi^2);
    vars_fos{end+1} = Ck;
    varvals_fos{end+1} = Cd_v/2;
end
Gfos = subs(Gfos,FCC0,65.6*3600);
h=figure(1004);
clf
set(h,'Name','Sensitivity function (Foster black)')
freq=logspace(-5,-1,100);

sensitibity_plot(Gfos,freq,vars_fos,varvals_fos);
ylim([-100,-20])
print('-depsc2','output/fos.eps')
print('-dpng','output/fos.png')
    

%% 感度関数のプロット（有限次元・有理関数・離散時間モデル）
h=figure(1005);
clf
set(h,'Name','Sensitivity function (blackbox model)')
freq=logspace(-5,-1,100);
sensitibity_plot(Grdt,freq,vars_rdt,varvals_rdt,dt);
ylim([-100,20])
drawnow
%end

%% 有限次元近似の精度確認
h=figure(1010);
clf
set(h,'Name','Finite-dimensional approximation')
bode_sym(G0,freq);
hold on
h=bode_sym(Gr0,freq);
set(h,'Color',[1,0,0]);
set(h,'LineStyle','--');
h=bode_sym(Gr0dt,freq,dt);
set(h,'Color',[0,1,0]);
set(h,'LineStyle','--');
legend('original','approximate','approximate(dt)')


function out = Zw_rec(Z,Y,k,m)
if k==m
  out=0;
else
  out = 1/(1/Z(k)+1/(1/Y(k)+Zw_rec(Z,Y,k+1,m)));
end
end

