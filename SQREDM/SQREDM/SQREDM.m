function Out = SQREDM(D,dim,pars)
%
% This code aims to solve the model
%
%   min_Z  || sqrt(H).*(sqrt(Z)-D) ||^2 + (rho/2) ||Z+P_Kr(-Z)||^2
%    s.t.    L<=Z<=U
%
%
% INPUTS:
%
%	D   : n-by-n dissimilarities matrix                          [required]
%         diag(D) = 0
%         dissimilarities are UNSQUARED, i.e.,
%                          D_ij=||point_i-point_j||+noise 
%         
%   dim : the embedding dimension  (e.g. dim = 2 or 3)           [required]
%
%   pars: parameters and other information                       [optional]
%         pars.m   : m  -- the number of given points, m>=0  
%         pars.PP  : PP -- dim-by-n matrix of coordinates of n points with
%                          first m(>=0) columns being given
%                    For sensor network localization (SNL)
%                    PP = [PA, PS]
%                    PA -- dim-by-m matrix of coordinates of m anchors
%                    PS -- dim-by-(n-m) matrix of coordinates of (n-m) sensors
%        pars.rho  : initial value of rho; (default rho=sqrt(n))  
%        pars.update:0  -- fix rho (default); 1 -- update rho during process       
%        pars.LOWBD: upper bound i.e., L=pars.LOWBD.^2, Z_{ij}>=L_{ij}^2 
%        pars.UPPBD: lower bound i.e., U=pars.LOWBD.^2, Z_{ij}<=U_{ij}^2
%                    Note: elements of pars.LOWBD and pars.UPPBD are UNSQUARED distances                          
%        pars.range: the communication range of two points, which means
%                    upper bound for Z_{ij}<=pars.range^2 if (D_{ij}>0  & i~=j)
%                    lower bound for Z_{ij}>=pars.range^2 if (D_{ij}==0 & i~=j)
%                    Note: pars.range is particular for SNL problem. If pars.range
%                    exists, no need pars.LOWBD and pars.UPPBD
%        pars.Otol : tolerance for objective, default Otol=sqrt(n)*1e-3 
%        pars.Etol : tolerance for eigenvalue, default Etol=1e-3  
%                    Note: If the noise is relatively large, change Etol=1e-2
%        pars.draw : 1--plot localizations in Re^dim (default); 0--no plot 
%
%
% OUTPUTS:
%
% If NO pars.PP exists 
%       Out.X:      dim-by-n matrix,  final coordinates 
%       Out.Time:   total time
%       Out.stress: relative stress
%
% If pars.PP exists 
%       Out.X:      dim-by-(n-m) matrix, coordinates before refinement 
%       Out.rX:     dim-by-(n-m) matrix, coordinates after refinement 
%       Out.Time:   total time including time for refinement
%       Out.rTime:  time for refinement
%       Out.RMSD:   Root Mean Square Distance (RMSD) before refinement 
%       Out.rRMSD:  Root Mean Square Distance (RMSD) after refinement 
%
% Refinement step is taken from Kim-Chaun Toh SNLSDP solver
%
% Send your comments and suggestions to   [ sz3g14@soton.ac.uk ]                       
%
% Warning: Accuracy may not be guaranteed!!!!!   
%
% This version: March 1st, 2018,   written by Shenglong Zhou   


t0=tic;

% parameters design
n  = size(D,1);
if nargin==2; pars=[]; end
[m,itmax,Eigtol,Objtol] = getparameters(n,pars);

if m>0; 
    fprintf('\nNumber of given points : %3d\n',m);
    fprintf('Number of unknown points: %3d\n',n-m);
    fprintf('Procrustes analysis and refinements step will be done!\n');
else
    fprintf('\nNo points are given!\n'); 
    fprintf('Number of unknown points: %3d\n',n);
end

Do    = D;
nD    = nnz(D);
rate  = nD/n/n;
% shortest path to complete the missing elements of D
    fprintf('Available dissimilarities rate: %1.2f \n',rate);
if  rate<0.05
    fprintf('Suggest providing more available dissimilarities!\n');
end
if rate<0.9
    D2 = D.*D;    
    disp('Contruct the shortest paths...');
    ts = tic;
    SD = graphallshortestpaths(sparse(sqrt(D2)));
    fprintf('Done the shortest paths by using %1.2f seconds\n',toc(ts));
    if any(SD == Inf)
    error('The neighborhood graph is not connected, increase range!');
    end
    SD = max(SD, D);  % to only replace those missing distances  
else
    SD = D;
    disp('No shortest paths calculated!');
end

fSD   = full(SD);
scale = max(fSD(:));
if scale<=10; scale=1; else D=D./scale; fSD=fSD./scale; end

H = full(spones(D)) ;
D  = full(D);
r  = dim;
T  = [];
if m>0; T = 1:m; H(T,T)=0;  end

% ilitialize ( L  U  Z rho) -----------------------------------------------
Z  = fSD.*fSD; 
UB = n*max(Z(:));
L  = zeros(n); 
U  = UB*ones(n);

if isfield(pars,'LOWBD');  
    L  = (pars.LOWBD/scale).^2;  
end
if isfield(pars,'UPPBD');  
    U  = (pars.UPPBD/scale).^2;
    HU = spones(tril(U,-1));
    if nnz(HU)<(n^2-n)/2; U=U+(1-HU-HU')*UB; end   % if U_{ij}=0 set U_{ij}=UB
    U(U==inf)=UB;
end 
if isfield(pars,'range');  
    H1 = 1-H; 
    rs = (pars.range/scale)^2; 
    L  = L.*H  + H1*rs;                           
    U  = U.*H1 + H*rs;
end

L(T,T)       = Z(T,T);      
U(T,T)       = Z(T,T); 
L(1:n+1:end) = 0;               
U(1:n+1:end) = 0;  
Z      = min( max(Z, L ), U );

rho    = sqrt(n);
update = 0; 
draw   = 1;
if isfield(pars,'rho');    rho    = pars.rho;    end
if isfield(pars,'update'); update = pars.update; end 
if isfield(pars,'draw');   draw   = pars.draw;   end

H2    = H.*H;  
H2D   = H2.*D;
H2r   = H2/rho;
H2Dr  = H2D/rho;
TH    = find(H2Dr>0);  
PZ    = ProjKr(-Z,r);

fprintf('Start to run ... \n');
fprintf('\n------------------------------------------------\n');
ErrEig  = Inf; 
FZr     = FNorm(sqrt(Z(TH))-D(TH))+rho*FNorm(Z+PZ)/2; 

for iter= 1:itmax  
    
    Z   = min( max( Lhalf(-PZ-H2r,H2Dr,TH), L ), U );
    PZ  = ProjKr(-Z,r);
    
    % stop criteria    
    ErrEig0 = ErrEig;
    gZ      = FNorm(Z+PZ)/2;
    ErrEig  = 2*gZ/FNorm(JXJ(Z));
    FZro    = FZr;
    FZr     = FNorm(sqrt(Z(TH))-D(TH))+rho*gZ;
    ErrObj  = abs(FZro-FZr)/(rho+FZro);
    fprintf('Iter: %3d  ErrEig: %.3e  ErrObj: %.3e\n',iter, ErrEig, ErrObj);
    
    if (ErrEig<Eigtol | abs(ErrEig0-ErrEig)<1e-8) & ...
       ErrObj<Objtol & iter>9; break; end 
    
    % update rho if update==1
    if update==1 & mod(iter,10)==0
        rho0 = rho;
        if ErrEig>Eigtol & ErrObj<Objtol/10; rho= 1.1*rho; end
        if ErrEig<Eigtol/10 & ErrObj>Objtol; rho= 0.9*rho; end  
        if rho0~=rho; H2r = H2/rho; H2Dr  = H2D/rho;       end
    end
    
end


% Results and graph output ------------------------------------------------
if isfield(pars, 'PP')  
    Ts        = 1+m:n;
    PPs       = pars.PP/scale;
    Xs        = procrustes_zhou(m, PPs, Z);            %procrustes analysis
    Out.Time  = toc(t0); tref=tic;   
    Xref      = refinepositions(Xs,PPs(:,T),[D(Ts,Ts) D(Ts,T)]);%refinement
    Out.rTime = toc(tref);  
    Out.Z     = Z*scale^2;
    Out.X     = Xs*scale;
    Out.rX    = Xref*scale;
    Out.Time  = Out.Time+Out.rTime;
    Out.RMSD  = sqrt(FNorm(pars.PP(:,Ts)-Out.X)/(n-m));
    Out.rRMSD = sqrt(FNorm(pars.PP(:,Ts)-Out.rX)/(n-m));   
    if  draw
        plot_SNL(Out.X,Out.rX,Out.RMSD,Out.rRMSD,pars);      
    end
    fprintf('------------------------------------------------\n');
    fprintf('rTime:   %1.3fsec\n',  Out.rTime)
    fprintf('Time:    %1.3fsec\n',  Out.Time)
    fprintf('RMSD:    %1.2e \n',    Out.RMSD )
    fprintf('rRMSD:   %1.2e \n',    Out.rRMSD)
else
    Out.Time   = toc(t0);
    [U,E]      = eig(JXJ(-Z)/2);
    Eig        = sort(real(diag(E)),'descend');
    Er         = real((Eig(1:r)).^0.5); 
    Out.Eigs   = Eig*(scale^2);
    Out.X      = (sparse(diag(Er))*real(U(:,1:r))')*scale;
    Z          = sqrt(Z)*scale;
    Out.stress = sqrt(FNorm(sqrt(Z(H~=0))-Do(H~=0))/FNorm(Do(H~=0)));
    fprintf('------------------------------------------------\n');
    fprintf('Time:      %1.3fsec\n',  Out.Time)
    fprintf('Stress:    %1.2e \n',    Out.stress)    
end

end

% ------------------------------------------------------------------------
function  [m,itmax,Eigtol,Objtol] = getparameters(n,pars)
    itmax  = 2000; 
    m      = 0;
    Objtol = sqrt(n)*1e-3;
    Eigtol = 1e-2;   
    if isfield(pars, 'Otol');   Objtol=pars.Otol;  end
    if isfield(pars, 'm');      m=pars.m;          end
    if isfield(pars, 'Etol');   Eigtol=pars.Etol;  end
end
 
% ------------------------------------------------------------------------
function  JXJ = JXJ( X )
% Compute J*X*J where J = I-ee'/n;
    nX   = size(X,1);
    Xe   = sum(X, 2);  
    eXe  = sum(Xe);    
    JXJ  = repmat(Xe,1,nX);
    JXJt = repmat(Xe',nX,1);
    JXJ  = -(JXJ + JXJt)/nX;
    JXJ  = JXJ + X + eXe/nX^2;
end

% ------------------------------------------------------------------------
function x =Lhalf( w, a, In )
% min 0.5*(x-w)^2-2*a*x^(0.5) s.t. x>=0. (where a>=0)
[mw,nw]    = size(w);
if nnz(In)/mw/nw < 0.4
	x      = w; 
    a2     = zeros(mw,nw); 
    w3     = zeros(mw,nw) ;
    d      = zeros(mw,nw) ;
    a2(In) = a(In).^2/4; 
    w3(In) = w(In).^3/27; 
    d(In)  = a2(In)-w3(In);

    I1=In(find(d(In)<0));
    if ~isempty(I1) 
        x(I1)=(2*w(I1)/3).*(1+cos((2/3)*acos(sqrt(a2(I1)./w3(I1)))));
    end
    
    I2=In(find( d(In)>=0 & w(In)>=0));
    if ~isempty(I2)
        st2=sqrt(d(I2));
        x(I2)= ((a(I2)/2+st2).^(1/3)+(a(I2)/2-st2).^(1/3)).^2;
    end
    
    I3=In(find( d(In)>=0 & w(In)<0));
    if ~isempty(I3)
        st3=sqrt(d(I3));
        x(I3)= ((a(I3)/2+st3).^(1/3)-(st3-a(I3)/2).^(1/3)).^2;
    end
else
    x=zeros(mw,nw);
    a2=a.^2/4; w3=w.^3/27; d=a2-w3;
    I1=find(d<0); 
    if ~isempty(I1) 
        x(I1)=(2*w(I1)/3).*(1+cos((2/3)*acos(sqrt(a2(I1)./w3(I1)))));
    end
    
    I2=find( d>=0 & w>=0);
    if ~isempty(I2)
        st2=sqrt(d(I2));
        x(I2)= ((a(I2)/2+st2).^(1/3)+(a(I2)/2-st2).^(1/3)).^2;
    end
    
    I3=find( d>=0 & w<0);
    if ~isempty(I3)
        st3=sqrt(d(I3));
        x(I3)= ((a(I3)/2+st3).^(1/3)-(st3-a(I3)/2).^(1/3)).^2;
    end
 
end
x = real(x);
end
% ------------------------------------------------------------------------
function fn= FNorm(A)
% Compute the Frobenius norm of A, i.e., ||A||_F^2
    fn=sum(sum(A.*A));
end

% ------------------------------------------------------------------------
function Z0= ProjKr(A,r)
% The projection of A on cone K_+^n(r)  
	JAJ     = JXJ(A);
	JAJ     = (JAJ+JAJ')/2
	[V0,P0] = eigs(JAJ,r,'LA');
	Z0      = real(V0*max(0,P0)*V0'+A-JAJ); 
end

% ------------------------------------------------------------------------
function Xs= procrustes_zhou(m0, PP, D0)
% Input:
% m0 -- the number of given anchors in Re^r0
% PP -- r0-by-n0 matrix whose first m0 columns are  given anchors 
%       and rest (n0-m0) columns are given sensors
% D0--  n0-by-n0 distance matrix 
%      [anchors-anchors (squred)distance, anchors-sensors (squred)distance
%       sensors-anchors (squred)distance, sensors-sensors (squred)distance]
%       containing m0 computed anchors and (n0-m0) sensors in Re^r0
% Output:
% Xs -- r0-by-(n0-m0) matrix containing (n0-m0) sensors in Re^r0
    r0     = size(PP,1);
    n0     = size(D0,2);
    JDJ    = JXJ(D0);
    JDJ    = -(JDJ+JDJ')/4;
    [U1,D1]= eigs(JDJ,r0,'LA');   
    X0     = (D1.^(1/2))*(U1');   
    if m0>0
    A      = PP(:,1:m0);
    [Q,~,a0,p0] = procrustes_qi(A,X0(:,1:m0));	
	Z0     = Q'*(X0-p0(:, ones(n0, 1))) + a0(:, ones(n0, 1)); 
    Xs     = Z0(:,m0+1:n0);
    Xa     = Z0(:, 1:m0);
    %Xs    = Xs*mean(sum(Xa.*A)./sum(Xa.^2));
    Xs     = Xs*max(1,sum(sum(Xa.*A))/sum(sum(Xa.^2)));
    else        
    [~,~,Map]= procrustes(PP',X0'); 
    Xs       = (Map.b*Map.T')*X0 + diag(Map.c(1,:))*ones(size(X0));
    end
end

% ------------------------------------------------------------------------
function [Q, P, a0, p0] = procrustes_qi(A, P)
% Procrustes analysis for rigid localization of anchors in A
% Input:
% A -- r-by-m0 matrix containing m0 anchors in Re^r
% P -- r-by-m0 matrix containing m0 computed locations of m0 anchors in Re^r
% 
% Output:
% Q  -- the orthogonal transformation
% P  -- the resulting coordinators of P after transformation
% a0 -- the tranlation vector that simply centrizes A
% p0 -- the tranlation vector that simply centrizes P
%
%
    m0 = size(A,2);
    % centerizing A and P
    a0 = sum(A, 2)/m0;
    p0 = sum(P,2)/m0;
    A  = A - a0(:, ones(m0,1));
    P1 = P - p0(:, ones(m0,1));
    P  = P1*A';
    [U0, ~, V] = svd(P);
    Q  = U0*V';
    P  = Q'*P1 + a0(:, ones(m0,1));
end
% ------------------------------------------------------------------------
function  plot_SNL(A,B,a,b,pars)
% Graph embedding ponits
figure,
if ~isfield(pars,'E')   
    PP=pars.PP;
    if isfield(pars, 'm'); m=pars.m; else m=0;  end
    [r,n]= size(PP);
    T   = 1:m;
    T1  = m+1:n;
    if r==2
        Af=[PP(:,T) A]; 
        Bf=[PP(:,T) B];
        for j=1:2
            if j==2; A=B; Af=Bf; a=b; end
            subplot(1,2,j),
            set(gca,'FontName','Times','FontSize',8)
            plot(PP(1,T),PP(2,T),'gs','markersize',4,'linewidth',2);     
            for i=1:n-m; 
                line([A(1,i) PP(1,m+i)], [A(2,i) PP(2,m+i)],...
                    'linewidth',.1,'color','b'); hold on
            end
            plot(PP(1,T1),PP(2,T1),'bo','markersize',4.5);hold on
            plot(A(1,:),A(2,:),'m*','markersize',3);hold on          
            plot(PP(1,T),PP(2,T),'gs','markersize',4,'linewidth',2);
            ZZ = [PP Af]';
            if j==1;
                xlabel(['Before refinement: RMSD = ', sprintf('%4.2e', a)],...
                    'FontName','Times','FontSize',8);
            else
                xlabel(['After refinement: rRMSD = ', sprintf('%4.2e', a)],...
                    'FontName','Times','FontSize',8);
            end
            axis([min(min(ZZ(:,1))) max(max(ZZ(:,1))) min(min(ZZ(:,2))) max(max(ZZ(:,2)))]) 
            hold on
        end
    else
        Af=[PP(:,T) A];
        Bf=[PP(:,T) B];  
        for j=1:2
            if j==2; A=B; Af=Bf; a=b; end
            subplot(1,2,j),       
            set(gca,'FontName','Times','FontSize',8)
           
            plot3(PP(1,T1),PP(2,T1),PP(3,T1),'bo','markersize',4.5);hold on
            plot3(A(1,:),A(2,:),A(3,:),'m*','markersize',3);hold on   
            plot3(PP(1,T),PP(2,T),PP(3,T),'gs','markersize',4,'linewidth',2); hold on            
            for i=1:n-m; 
                line([A(1,i) PP(1,m+i)], [A(2,i) PP(2,m+i)], [A(3,i) PP(3,m+i)],...
                    'linewidth',.1,'color','b'); hold on
            end
            ZZ = [PP Af]';
            plot3(PP(1,T),PP(2,T),PP(3,T),'gs','markersize',4,'linewidth',2); hold on        
            if j==1;
                title(['Before refinement: RMSD = ', sprintf('%4.2e', a)],...
                    'FontName','Times','FontSize',8);
            else
                title(['After refinement: rRMSD = ', sprintf('%4.2e', a)],...
                    'FontName','Times','FontSize',8);
            end
            grid on;
            axis([min(min(ZZ(:,1))) max(max(ZZ(:,1))) min(min(ZZ(:,2))) max(max(ZZ(:,2)))...
                min(min(ZZ(:,3))) max(max(ZZ(:,3)))]) 
            hold on             
        end
    end    
else
    if a <b; B=A; end
    Xfrefine=[pars.PP(:,1:pars.m) B];
    scatter(Xfrefine(1,pars.E),Xfrefine(2,pars.E),12,'filled');
	hold on
	scatter(Xfrefine(1,pars.D),Xfrefine(2,pars.D),12,'filled');
	scatter(Xfrefine(1,pars.M),Xfrefine(2,pars.M),12,'filled');
 	scatter(Xfrefine(1,1:pars.m),Xfrefine(2,1:pars.m),12,'filled','k');
    set(gca,'FontName','Times','FontSize',8);
    if a <b;          
        title(['Before refinement: RMSD = ',  sprintf('%4.2e', a)],...
                    'FontName','Times','FontSize',9);  
    else
        title(['After refinement: rRMSD = ', sprintf('%4.2e', b)],...
                    'FontName','Times','FontSize',9);
    end
	axis equal;
	axis([-5 120 -5 55]);
end
    
end

