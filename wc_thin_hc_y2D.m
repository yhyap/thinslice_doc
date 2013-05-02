% PDE of convection-diffusion problem in 1D-2D cylindrical coordinate
% Isothermal and flat entering velocity profile
% when (j=1,r=0), (jj=1,s=0)

% Clear previous files
clear all
% Clear command window or clean screen
clc
format long

%% Global area
global pi nr nz ns nu dr ds du dz drs dss dus dzs r r0 s u u0 z P Mco Mair Mo2 Mco2 Mc3h8 Vco Vair Vco2 Vc3h8 re por tau Ya Yb Yc Tk Yas Ybs Ycs Tks Tku Yae Ybe Yce Tke fcv fcCO fcO fcOO HisRa HisRb HisRc HisRd HisRe HisRf HisTks HisTke G H Av r_gtc rho_wc rho_cord R vmean ncall History_t

%% Type of sample
% 1=JM, 2=Degussa, 3=Ford
sample=1;

%% Operating parameters and inlet conditions
% a represents CO, b represents CO2
Ya0=3000/1e6;   % ppm (initial CO mol fraction)
Yae=3000/1e6;   % ppm (entering CO mol fraction)
Yb0=0.0;        % ppm (initial CO2 mol fraction)
Ybe=0.0;        % ppm (entering CO2 mol fraction)
Yc0=500/1e6;    % ppm (initial C3H8 mol fraction)
Yce=500/1e6;    % ppm (entering C3H8 mol fraction)
Tk0=417;        % K (initial temperature)
Tke=417;        % K (entering temperature) @ 190 C
R=8.3145;       % m3.Pa/(K.mol) (gas constant)
P=101325;       % Pa (pressure) (assume isobaric throughout simulation)
vmean=2.4;      % m/s (linear mean fluid velocity) (default:0.06)
Mflow=0.64;     % mol/s, average molar flow rate of fuel+air

%% Dimension of monolith & active sites parameters
% Information of dimension and active sites for different samples
if (sample==1)
    % Dimension
    s0=100e-6;      % m (washcoat thickness)
    u0=90e-6;       % m (cordierite thickness)
    r0=600e-6-s0-u0;% m (channel radius)
    dm=0.106;       % m (diameter of monolith block)
    w_slice=27.6316;% g, weight of a monolith slice based on 630 g per mon
    cpsi=400;       % cells per square inch
    % Details of porous media (Assume parallel pore)
    re=12e-9;       % equivalent pore radius in m, 13.33e-9 from estimation
    por=0.55;       % Porosity, m3 gas/m3 total
    tau=1;          % Tortuosity
    % Active sites per catalyst/washcoat
    LPt=2.0e-6/0.3; % mol(Pt)/g(cat) Pt loading data from CO chemisorption - on wc & cord  (default: 2.2e-6)
    aBET=130;       % m2(cat)/g(cat) (BET surface area from ASAP & Autopore)
    rho_wc=1.3e6;   % g/m3 loose bulk density from Disperal data sheet, sensitivity: 500, 700 g/l, def: 0.6e3
    rho_cord=2.5e6; % g/m3 density of substrate from H&K pg.655, also measured using pycnometer
    LPtc=0.005468;  % g(Pt)/g(cat) data from weighing, JM claim: 120g/ft3 = 0.006763
    %LPtc=0.006763;
elseif (sample==2)
    % Dimension
    dm=0.116;       % m (diameter of monolith block)
    cpsi=400;       % cells per square inch
    % Details of porous media
    re=10e-9;       % equivalent pore radius in m
    por=0.4;        % Porosity, m3 gas/m3 total
    tau=4;          % Tortuosity
    % Active sites per catalyst/washcoat
    LPt=7.6e-6;     % mol(Pt)/g(cat) Pt loading data from CO chemisorption - on wc & cord
    aBET=100;       % m2(cat)/g(cat) (BET surface area from ASAP & Autopore)
    rho_wc=6e5;     % g/m3 loose bulk density from Disperal data sheet, sensitivity: 500, 700 g/l
    LPtc=0.005468;  % g(Pt)/g(cat) data from weighing, Degussa claim: 120g/ft3 = 0.006763
            % Degussa has 80 or 90 g/ft^3 of monolith
else
    % Dimension
    dm=0.118;       % m (diameter of monolith block)
    cpsi=400;       % cells per square inch
    % Active sites per catalyst/washcoat
    
end
% Common information of catalysts and further calculations
pi=3.141593;        % pi
zl=0.005;           % m (channel length) - 5 mm thin slice
Avo=6.022e23;       % mol-1 Avogadro number
MPt=195.08;         % g(Pt)/mol(Pt) molar mass of Pt
a_m=8.07e-20;       % m2 surface area occupied by a Pt atom on a polycrystalline surface (assume f.c.c.)from handbook
cpsm=cpsi/0.000645; % cells per square meter
Amon=pi*(dm/2)^2;   % m2, cross section area of monolith slice
ncell=cpsm*Amon;    % no. of cells
V_slice=Amon*zl;    % m3, volume of a slice including voidage 
H=LPt/aBET;         % mol(Pt) m-2(cat) (assume uniform active sites density)
Av=aBET*rho_wc;     % m2(cat) m-3(cat) (internal surface area per reactor vol)
                    % Mind Samolons define this as m2(Pt) m-3(cat)
r_gtc=por/(1-por);  % m3 (gas) m-3 (cat) ratio of volume of gas to solid
N_atom=LPtc/MPt*Avo;% no. of Pt atoms as deposited
Ns_atom=LPt*Avo;    % no. of Pt surface atoms (active for reaction)
D=Ns_atom/N_atom;   % Pt dispersion
%D=0.22;             % D=0.22, G=9.35e-5; D=1, G=2.06e-5; D=0.75, G=2.74e-5;
Ssp=a_m*Avo/MPt*D;  % m2(Pt) g-1(Pt) Specific surface area
%G=1/(Ssp*MPt);      % mol(Pt) m-2(Pt) from its reciprocal of m2(Pt) mol-1(Pt): 2.88e-4 when D=0.071
                    % Mind Salomons define this as 2.72e-5 mol(Pt) m-2(Pt): too high, mine is 2.2e-8
                    % Herz & Marin defined this as 6.835e-2mol(Pt) m-2(Pt) from CO chemisorption 7.5e2 cm2 g-1
G=2.0573e-5;
                    
%% Parameters for bulk diffusion coefficient (Stefan-Maxwell hard sphere model)
Mco=28.01;      % Molar mass of CO in kg/kmol
Mair=28.96;     % Molar mass of air in kg/kmol
Mo2=32.00;      % Molar mass of oxygen in kg/kmol
Mco2=44.01;     % Molar mass of carbon dioxide in kg/kmol
Mc3h8=44.1;     % Molar mass of propane in kg/kmol
Vco=18.9;       % Molecular volume of CO in m3/kmol in liq form at its b.p.
Vair=20.1;      % Molecular volume of air in m3/kmol in liq form at its b.p.
Vco2=22.262;    % Molecular volume of CO2 in m3/kmol in liq form at its b.p. (CHECK this in H&K)
Vc3h8=65.34;    % Molecular volume of C3H8 in m3/kmol (est. from H&K pg 229)

%% Initial fractional coverage of active sites
thetav=0;%0.00272; % vacant sites
thetaCO=1;%0.99228;% CO sites
thetaO=0;%0.005;   % O sites
thetaOO=0;      % OO sites

%% Setting up grids
    %% Grid in axial direction (gas and solid phase)
    nz=3;  %40
    dz=zl/nz;
    for i=1:nz
        z(i)=i*dz;
    end
    dzs=dz^2;

    %% Grid in radial direction (gas phase)
    nr=5;  % Change this and change RCOtot
    dr=r0/(nr-1);
    for j=1:nr
        r(j)=(j-1)*dr;
    end
    drs=dr^2;
    % jj=1 at boundary between solid and gas phase
    
    %% Grid in radial direction (solid phase)
    ns=3;  % Change this and change RCOtot
    ds=s0/(ns-1);
    for jj=1:ns
        s(jj)=(jj-1)*ds;
    end
    dss=ds^2;
    % jj=1 at boundary between solid and gas phase

    %% Grid in radial direction (solid cordierite phase)
    nu=3;  
    du=u0/(nu-1);
    for jjj=1:nu
        u(jjj)=(jjj-1)*du;
    end
    dus=du^2;
    % jjj=1 at boundary between washcoat and cordierite

%% Independent variable for ODE integration
tf=1800.0;
td=1.0;
tout=[0.0:td:tf]';
nout=tf;
ncall=0;

%% Initial condition
% Assuming flat inlet profile
for i=1:nz
    % Gas phase
    for j=1:nr    
        Ya(i,j)=Ya0;
        Yb(i,j)=Yb0;
        Yc(i,j)=Yc0;
        Tk(i,j)=Tk0;
        % Conversion of 1D ICs to 1D vector as ODE integrators only deal
        % with 1D dependent variables. This can be done in two ways:
        % 1) explicit programming of matrix subscripting
        % 2) the use of MATLAB's reshape utility
        y0((i-1)*nr+j)=Ya(i,j);
        y0((i-1)*nr+j+nz*nr)=Yb(i,j);
        y0((i-1)*nr+j+2*nz*nr)=Yc(i,j);
        y0((i-1)*nr+j+3*nz*nr)=Tk(i,j);
    end
    % Solid washcoat phase
    for jj=1:ns
        Yas(i,jj)=Ya0;
        Ybs(i,jj)=Yb0;
        Ycs(i,jj)=Yc0;
        Tks(i,jj)=Tk0;
        fcv(i,jj)=thetav;
        fcCO(i,jj)=thetaCO;
        fcO(i,jj)=thetaO;
        fcOO(i,jj)=thetaOO;
        % Conversion to 1D vector
        y0(4*nz*nr+(i-1)*ns+jj)=Yas(i,jj);
        y0(4*nz*nr+(i-1)*ns+jj+nz*ns)=Ybs(i,jj);
        y0(4*nz*nr+(i-1)*ns+jj+2*nz*ns)=Ycs(i,jj);
        y0(4*nz*nr+(i-1)*ns+jj+3*nz*ns)=Tks(i,jj);
        y0(4*nz*nr+(i-1)*ns+jj+4*nz*ns)=fcv(i,jj);
        y0(4*nz*nr+(i-1)*ns+jj+5*nz*ns)=fcCO(i,jj);
        y0(4*nz*nr+(i-1)*ns+jj+6*nz*ns)=fcO(i,jj);
        y0(4*nz*nr+(i-1)*ns+jj+7*nz*ns)=fcOO(i,jj);
    end
    % Solid cordierite phase
    for jjj=1:nu
        Tku(i,jjj)=Tk0;
        %Conversion to 1D vector
        y0(4*nz*nr+8*nz*ns+(i-1)*nu+jjj)=Tku(i,jjj);
    end
end

%% ODE integration
reltol=1.0e-5;     % default 1e-14
abstol=1.0e-5;     % default 1e-14
options=odeset('RelTol',reltol,'AbsTol',abstol);
mf=1;
if(mf==1)
    [t,y]=ode15s(@pde_1wc_thin_hc_y2D,tout,y0,options);
end

%% 1D to 2D matrices
for it=1:nout
    for i=1:nz
        % Gas phase
        for j=1:nr
            Ya(it,i,j)=y(it,(i-1)*nr+j);
            Yb(it,i,j)=y(it,(i-1)*nr+j+nz*nr);
            Yc(it,i,j)=y(it,(i-1)*nr+j+2*nz*nr);
            Tk(it,i,j)=y(it,(i-1)*nr+j+3*nz*nr);
            % This returns a 3D matrix
        end
        % Solid washcoat phase
        for jj=1:ns
            Yas(it,i,jj)=y(it,4*nz*nr+(i-1)*ns+jj);
            Ybs(it,i,jj)=y(it,4*nz*nr+(i-1)*ns+jj+nz*ns);
            Ycs(it,i,jj)=y(it,4*nz*nr+(i-1)*ns+jj+2*nz*ns);
            Tks(it,i,jj)=y(it,4*nz*nr+(i-1)*ns+jj+3*nz*ns);
            fcv(it,i,jj)=y(it,4*nz*nr+(i-1)*ns+jj+4*nz*ns);
            fcCO(it,i,jj)=y(it,4*nz*nr+(i-1)*ns+jj+5*nz*ns);
            fcO(it,i,jj)=y(it,4*nz*nr+(i-1)*ns+jj+6*nz*ns);
            fcOO(it,i,jj)=y(it,4*nz*nr+(i-1)*ns+jj+7*nz*ns);
            % This returns a 3D matrix
        end
        % Solid cordierite phase
        for jjj=1:nu
            Tku(it,i,jjj)=y(it,4*nz*nr+8*nz*ns+(i-1)*nu+jjj);
            % This returns a 3D matrix
        end
    end
end

%% Computation of reaction rate (from amount reacted)
for it=1:nout
    rCO(it)=Mflow*(Yae-Ya(it,nz,1));            % CO Reaction rate in mol/s, assume plug flow
    r2CO(it)=rCO(it)/V_slice;                   % CO Reaction rate per slice
    COconv(it)=(Yae-Ya(it,nz,1))/Yae;           % CO conversion
    rHC(it)=Mflow*(Yce-Yc(it,nz,1));            % HC Reaction rate in mol/s
    r2HC(it)=rHC(it)/V_slice;                   % HC Reaction rate per slice
    HCconv(it)=(Yce-Yc(it,nz,1))/Yce;           % HC conversion
end

%% Computation of inlet temperature profile
% Varying temperature to see light off and hysteresis
        for it=1:nout
            rt1=200;                % linear ramp time in second
            rt2=500;                % linear ramp time in second
            rt3=300;                % 300 s of constant temperature
            rt4=100;                % linear ramp time in second
            rt5=200;                % linear ramp time in second
            rt6=500;                % linear ramp time in second
            Tkef1=505;              % temperature point 1
            Tkef2=545;              % temperature point 2, max
            Tkef4=485;              % temperature point 3
            Tkef5=445;              % temperature point 4
            mp1=(Tkef1-Tke)/rt1;    % multiplier 1
            mp2=(Tkef2-Tkef1)/rt2;  % multiplier 2
            % mp3 is flat therefore no equation
            mp4=(Tkef4-Tkef2)/rt4;  % multiplier 4
            mp5=(Tkef5-Tkef4)/rt5;  % multiplier 5
            mp6=(Tke-Tkef5)/rt6;    % multiplier 6
            if it<rt1                
                Tkee(it)=Tke+mp1*it;
            elseif it<(rt1+rt2)      
                Tkee(it)=Tkef1+mp2*(it-rt1);
            elseif it<(rt1+rt2+rt3)       
                Tkee(it)=Tkef2;
            elseif it<(rt1+rt2+rt3+rt4)                
                Tkee(it)=Tkef2+mp4*(it-(rt1+rt2+rt3));
            elseif it<(rt1+rt2+rt3+rt4+rt5)
                Tkee(it)=Tkef4+mp5*(it-(rt1+rt2+rt3+rt4));
            elseif it<(rt1+rt2+rt3+rt4+rt5+rt6)
                Tkee(it)=Tkef5+mp6*(it-(rt1+rt2+rt3+rt4+rt5));
            else
                Tkee(it)=Tke;
            end
        end

   
%% Display a heading and centerline output
fprintf('\n nr = %2d  ns = %2d  nz = %2d\n',nr,ns,nz);

fprintf('\n ncall = %5d\n',ncall);
fprintf('\n nodes = %5d\n',nz*(1+ns));
fprintf('\n dt = %5d\n',td);

%% Parametric plots

figure(2);
plot(1:nout,Tk(1:nout,1,1),'g')
hold
plot(1:nout,Tks(1:nout,3,3),'r-')
plot(1:nout,Tku(1:nout,2,2),'m-');
plot(1:nout,Tkee(1:nout),'b');
axis tight
title('Tk & Tks');
xlabel('t');
ylabel('temperature');
legend('Tk(1,1)','Tks(3,3)','Tku(2,2)','entering T');

figure(4);
plot(1:nout,Yb(1:nout,nz,1),'g')
hold
plot(1:nout,Ya(1:nout,nz,1),'r')
plot(1:nout,Yc(1:nout,nz,1),'b');
axis tight
title('Y(r=0,z=L,t)');
xlabel('t');
ylabel('mol fraction');
legend('CO2 mol fraction','CO mol fraction','HC mol fraction');

%figure(5);
%plot(1:nout,Ya(1:nout,nz,1),'g')
%axis tight
%title('Ya(r=0,z=L,t)');
%xlabel('t');
%ylabel('CO mol fraction');

%figure(6);
%plot(z,Tk(nout,1:nz),'r-');
%hold
%plot(z,Tks(nout,1:nz,1),'b-');
%plot(z,Tks(nout,1:nz,5),'g-');
%axis tight;
%title('Tk & Tks along z at t=tf');
%xlabel('z');
%ylabel('temperature');
%legend('Tk','Tks1','Tks5');

%figure(7);
%plot(z,cas(1,1:nz,1),'r-');
%hold
%plot(z,cas(nout,1:nz,1),'b-');
%plot(z,ca(nout,1:nz),'g-');
%axis tight
%title('cas along z at different time');
%xlabel('z');
%ylabel('concentration');
%legend('cas-1','cas-tf','ca-tf');

%figure(8);
%plot(z,cbs(1,1:nz,1),'r-');
%hold
%plot(z,cbs(nout,1:nz,1),'co');
%plot(z,cb(nout,1:nz),'g-');
%axis tight
%title('cbs along z at different time');
%xlabel('z');
%ylabel('concentration');
%legend('cbs-1','cbs-tf','cb-tf');

%figure(9);
%plot(z,cbs(nout,1:nz,1),'r-');
%hold
%plot(z,cbs(nout,1:nz,3),'y-');
%plot(z,cbs(nout,1:nz,5),'b-');
%plot(z,cbs(nout,1:nz,7),'g-');
%plot(z,cbs(nout,1:nz,ns),'m-');
%axis tight
%title('cbs along z at tf');
%xlabel('z');
%ylabel('concentration');
%legend('1','3','5','7','ns');

figure(10);

subplot(2,2,1)
plot(1:nout,fcCO(1:nout,2,2),'r');
axis ([0 1800 0 1]);
title('fcCO');
xlabel('t');

subplot(2,2,2)
plot(1:nout,fcO(1:nout,2,2),'b');
axis ([0 1800 0 1]);
title('fcO');
xlabel('t');

subplot(2,2,3)
plot(1:nout,fcOO(1:nout,2,2),'m');
axis ([0 1800 0 1]);
title('fcOO');
xlabel('t');

subplot(2,2,4)
plot(1:nout,fcv(1:nout,2,2),'g');
axis ([0 1800 0 1]);
title('fcv');
xlabel('t');

%figure(12);
%plot(1:nout,rCO(1:nout));
%axis tight;
%title('Reaction rate from conversion');
%xlabel('t');
%ylabel('rate');

figure(13);
plot(1:nout,r2CO(1:nout),'r-');
hold
plot(1:nout,r2HC(1:nout),'b-');
axis tight;
title('Reaction rate per slice');
xlabel('t');
ylabel('rate mol/m3s');

figure(14);
plot(Tks(1:1000,2,2),r2CO(1:1000),'r-');
hold
plot(Tks(1001:nout,2,2),r2CO(1001:nout),'b-');
axis tight;
title('light-off hysteresis');
xlabel('solid temperature in K');
ylabel('CO reaction rate in mol/m3s of reactor volume');
legend('ignition','extinction');

figure(15);
plot(1:nout,COconv(1:nout),'r-');
hold
plot(1:nout,HCconv(1:nout),'b-');
axis ([0 1500 0 1]);
title('CO & HC conversion against time');
xlabel('time');
ylabel('CO & HC conversion');
legend('CO','HC');

figure(16);
plot(Tks(1:1000,2,2),COconv(1:1000),'r-');
hold
plot(Tkee(1:1000),COconv(1:1000),'rx');
plot(Tks(1001:nout,2,2),COconv(1001:nout),'b-');
plot(Tkee(1001:nout),COconv(1001:nout),'bx');
axis ([400 600 0 1.0]);
title('light-off hysteresis');
xlabel('temperature');
ylabel('CO conversion');
legend('ignition-solid','ignition','extinction-solid','extinction');

figure(17);

subplot(2,2,1)
plot(HisTks,HisRa,'r-');
axis ([400 600 0 15]);
title('HisRa vs Tks, CO ads');
xlabel('Tks');
ylabel('Ra');

subplot(2,2,2)
plot(HisTks,HisRb,'b-');
axis ([400 600 0 15]);
title('HisRb vs Tks, O2 ads');
xlabel('Tks');
ylabel('Rb');

subplot(2,2,3)
plot(HisTks,HisRc,'g-');
axis ([400 600 0 15]);
title('HisRc vs Tks, L-H');
xlabel('Tks');
ylabel('Rc');

subplot(2,2,4)
plot(HisTks,HisRd,'m-');
axis tight;
title('HisRd vs Tks, compress');
xlabel('Tks');
ylabel('Rd');

figure(18);

subplot(2,2,1)
plot(Tks(1:1000,2,2),fcCO(1:1000,2,2),'r-');
hold
plot(Tks(1001:nout,2,2),fcCO(1001:nout,2,2),'b-');
axis ([400 600 0 1]);
xlabel('solid temperature in K');
ylabel('CO coverage');
legend('ignition','extinction');

subplot(2,2,2)
plot(Tks(1:1000,2,2),fcO(1:1000,2,2),'r-');
hold
plot(Tks(1001:nout,2,2),fcO(1001:nout,2,2),'b-');
axis ([400 600 0 1]);
xlabel('solid temperature in K');
ylabel('O coverage');

subplot(2,2,3)
plot(Tks(1:1000,2,2),fcOO(1:1000,2,2),'r-');
hold
plot(Tks(1001:nout,2,2),fcOO(1001:nout,2,2),'b-');
axis ([400 600 0 1]);
xlabel('solid temperature in K');
ylabel('OO coverage');

subplot(2,2,4)
plot(Tks(1:1000,2,2),fcv(1:1000,2,2),'r-');
hold
plot(Tks(1001:nout,2,2),fcv(1001:nout,2,2),'b-');
axis tight;
xlabel('solid temperature in K');
ylabel('vacant coverage');

figure(19);
title('Individual reaction rate at (2,2)');

subplot(2,2,1)
plot(History_t,HisRa,'r-');
axis ([0 1800 0 15]);
title('CO ads-des');
xlabel('time');

subplot(2,2,2)
plot(History_t,HisRb,'b-');
axis ([0 1800 0 15]);
title('O2 ads');
xlabel('time');

subplot(2,2,3)
plot(History_t,HisRc,'g-');
axis ([0 1800 0 15]);
title('L-H');
xlabel('time');

subplot(2,2,4)
plot(History_t,HisRd,'m-');
axis ([0 1800 0 15]);
title('compress');
xlabel('time');

%subplot(2,3,5)
%plot(History_t,HisRe,'c-');
%axis ([0 1500 0 15]);
%title('LH CO-OO');
%xlabel('time');

%subplot(2,3,6)
%plot(History_t,HisRf,'y-');
%axis ([0 1500 0 15]);
%title('OO->O');
%xlabel('time');

figure(20);
plot(Tks(1:1000,2,2),Ya(1:1000,1,1),'r-');
hold
plot(Tks(1000:nout,2,2),Ya(1000:nout,1,1),'rx');
plot(Tks(1:1000,2,2),Ya(1:1000,2,1),'b-');
plot(Tks(1000:nout,2,2),Ya(1000:nout,2,1),'bx');
plot(Tks(1:1000,2,2),Ya(1:1000,nz,1),'g-');
plot(Tks(1000:nout,2,2),Ya(1000:nout,nz,1),'gx');
axis tight;
title('Ya during light off');
xlabel('solid temperature in K');
ylabel('CO mol fraction');
legend('ignition:1','extinction:1','ignition:2','extinction:2','ignition:nz','extinction:nz');

figure(21);
plot(Tks(1:1000,2,2),Yas(1:1000,1,3),'r-');
hold
plot(Tks(1000:nout,2,2),Yas(1000:nout,1,3),'rx');
plot(Tks(1:1000,2,2),Yas(1:1000,nz,3),'b-');
plot(Tks(1000:nout,2,2),Yas(1000:nout,nz,3),'bx');
plot(Tks(1:1000,2,2),Yas(1:1000,1,1),'g-');
plot(Tks(1000:nout,2,2),Yas(1000:nout,1,1),'gx');
plot(Tks(1:1000,2,2),Yas(1:1000,nz,1),'y-');
plot(Tks(1000:nout,2,2),Yas(1000:nout,nz,1),'yx');
plot(Tks(1:1000,2,2),Yas(1:1000,nz,3),'m-');
plot(Tks(1000:nout,2,2),Yas(1000:nout,nz,3),'mx');
axis tight;
title('Yas during light-off');
xlabel('temperature');
ylabel('CO mol fraction');
legend('ignition:1,3','extinction:1,3','ignition:nz,3','extinction:nz,3','ignition:1,1','extinction:1,1','ignition:nz,1','extinction:nz,1','ignition:nz,3','extinction:nz,3');

% Radial profiles (at t=1000)
for i=1:nz
    for j=1:nr
        Ya_rad(i,j)=Ya(1000,i,j);
        Tk_rad(i,j)=Tk(1000,i,j);
    end
    for jj=1:ns
        Yas_rad(i,jj)=Yas(1000,i,jj);
        Tks_rad(i,jj)=Tks(1000,i,jj);
    end
    for jjj=1:nu
        Tku_rad(i,jjj)=Tku(1000,i,jjj);
    end
end

figure(22);

subplot(2,2,1)
plot(z,Tk_rad(1:nz,2),'r-');
hold
plot(z,Tks_rad(1:nz,2),'g-');
plot(z,Tku_rad(1:nz,2),'b-');
xlabel('z - axial distance');
ylabel('Temperature');
title('Axial gas & solid phase temperature profile');
legend('Tk','Tks','Tku');

%subplot(2,2,2)
%plot(z,Tks_rad(1:nz,2),'r-');
%hold
%plot(z,Tku_rad(1:nz,2),'b-');
%xlabel('z - axial distance');
%ylabel('Tks & Tku');
%title('Axial solid phase temperature profile');
%legend('wc','cord');

subplot(2,2,3)
surf(s,z,Tks_rad)
axis([0 s0 0 zl 400 700]);
xlabel('s - radial distance');
ylabel('z - axial distance');
zlabel('Tks');
title(['surface plot of temperature at t=700s']);
view(-130,62);
rotate3d on

subplot(2,2,4)
surf(u,z,Tku_rad)
axis([0 u0 0 zl 400 700]);
xlabel('u - radial distance');
ylabel('z - axial distance');
zlabel('Tku');
title(['surface plot of temperature at t=700s']);
view(-130,62);
rotate3d on

figure(23);

subplot(2,2,1)
plot(z,Ya_rad(1:nz,1),'r-');
hold
plot(z,Ya_rad(1:nz,3),'g-');
plot(z,Ya_rad(1:nz,5),'m-');
plot(z,Yas_rad(1:nz,1),'c-');
plot(z,Yas_rad(1:nz,2),'b-');
xlabel('z - axial distance');
ylabel('mol fraction');
title('Axial gas phase concentration profile at time 1000s');
legend('Ya1','Ya3','Ya5','Yas1','Yas2');

%subplot(2,2,2)
%plot(z,Yas_rad(1:nz,2),'r-');
%xlabel('z - axial distance');
%ylabel('Yas');
%title('Axial solid phase concentration profile');

subplot(2,2,3)
surf(r,z,Ya_rad)
axis([0 r0 0 zl 0 0.003]);
xlabel('r - radial distance');
ylabel('z - axial distance');
zlabel('Ya');
title(['surface plot of concentration at t=1000s']);
view(-130,62);
rotate3d on

subplot(2,2,4)
surf(s,z,Yas_rad)
axis([0 s0 0 zl 0 0.003]);
xlabel('s - radial distance');
ylabel('z - axial distance');
zlabel('Yas');
title(['surface plot of concentration at t=1000s']);
view(-130,62);
rotate3d on