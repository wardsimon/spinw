function spectra = spinwave_structureFact(obj, hkl, omega, V, varargin)
% calculates spin correlation function using linear spin wave theory
%
% ### Syntax
%
% `spectra = spinwave(obj,Q)`
%
% `spectra = spinwave(___,Name,Value)`
%
% ### Description
%
% `spinwave(obj,Q,Name,Value)` calculates spin wave dispersion and
% spin-spin correlation function at the reciprocal space points $Q$. The
% function can solve any single-k magnetic structure exactly and any
% multi-k magnetic structure appoximately and quadratic spinw-spin
% interactions as well as single ion anisotropy and magnetic field.
% Biquadratic exchange interactions are also implemented, however only for
% $k_m=0$ magnetic structures.
%
% If the magnetic ordering wavevector is non-integer, the dispersion is
% calculated using a coordinate system rotating from unit cell to unit
% cell. In this case the spin Hamiltonian has to fulfill this extra
% rotational symmetry which is not checked programatically.
%
% Some of the code of the function can run faster if mex files are used. To
% switch on mex files, use the `swpref.setpref('usemex',true)` command. For
% details see the [sw_mex] and [swpref.setpref] functions.
%
% ### Examples
%
% To calculate and plot the spin wave dispersion of the
% triangular lattice antiferromagnet ($S=1$, $J=1$) along the $(h,h,0)$
% direction in reciprocal space we create the built in triangular lattice
% model using `sw_model`.
%
% ```
% >>tri = sw_model('triAF',1)
% >>spec = tri.spinwave({[0 0 0] [1 1 0]})
% >>sw_plotspec(spec)
% >>snapnow
% ```
%
% ### Input Arguments
%
% `obj`
% : [spinw] object.
%
% `Q`
% : Defines the $Q$ points where the spectra is calculated, in reciprocal
%   lattice units, size is $[3\times n_{Q}]$. $Q$ can be also defined by
%   several linear scan in reciprocal space. In this case `Q` is cell type,
%   where each element of the cell defines a point in $Q$ space. Linear scans
%   are assumed between consecutive points. Also the number of $Q$ points can
%   be specified as a last element, it is 100 by defaults.
%
%   For example to define a scan along $(h,0,0)$ from $h=0$ to $h=1$ using
%   200 $Q$ points the following input should be used:
%   ```
%   Q = {[0 0 0] [1 0 0]  200}
%   ```
%
%   For symbolic calculation at a general reciprocal space point use `sym`
%   type input.
%
%   For example to calculate the spectrum along $(h,0,0)$ use:
%   ```
%   Q = [sym('h') 0 0]
%   ```
%   To calculate spectrum at a specific $Q$ point symbolically, e.g. at
%   $(0,1,0)$ use:
%   ```
%   Q = sym([0 1 0])
%   ```
%
% ### Name-Value Pair Arguments
%
% `'formfact'`
% : If true, the magnetic form factor is included in the spin-spin
%   correlation function calculation. The form factor coefficients are
%   stored in `obj.unit_cell.ff(1,:,atomIndex)`. Default value is `false`.
%
% `'formfactfun'`
% : Function that calculates the magnetic form factor for given $Q$ value.
%   value. Default value is `@sw_mff`, that uses a tabulated coefficients
%   for the form factor calculation. For anisotropic form factors a user
%   defined function can be written that has the following header:
%   ```
%   F = formfactfun(atomLabel,Q)
%   ```
%   where the parameters are:
%   * `F`           row vector containing the form factor for every input
%                   $Q$ value
%   * `atomLabel`   string, label of the selected magnetic atom
%   * `Q`           matrix with dimensions of $[3\times n_Q]$, where each
%                   column contains a $Q$ vector in $\\ang^{-1}$ units.
%
% `'gtensor'`
% : If true, the g-tensor will be included in the spin-spin correlation
%   function. Including anisotropic g-tensor or different
%   g-tensor for different ions is only possible here. Including a simple
%   isotropic g-tensor is possible afterwards using the [sw_instrument]
%   function.
%
% `'fitmode'`
% : If `true`, function is optimized for multiple consecutive calls (e.g.
%   the output spectrum won't contain the copy of `obj`), default is
%   `false`.
%
% `'notwin'`
% : If `true`, the spectra of the twins won't be calculated. Default is
% `false`.
%
% `'sortMode'`
% : If `true`, the spin wave modes will be sorted by continuity. Default is
%   `true`.
%
% `'optmem'`
% : Parameter to optimise memory usage. The list of Q values will be cut
%   into `optmem` number of pieces and will be calculated piece by piece to
%   decrease peak memory usage. Default value is 0, when the number
%   of slices are determined automatically from the available free memory.
%
% `'tol'`
% : Tolerance of the incommensurability of the magnetic ordering wavevector.
%   Deviations from integer values of the ordering wavevector smaller than
%   the tolerance are considered to be commensurate. Default value is
%   $10^{-4}$.
%
% `'omega_tol'`
% : Tolerance on the energy difference of degenerate modes when
%   diagonalising the quadratic form, default value is $10^{-5}$.
%
% `'hermit'`
% : Method for matrix diagonalization with the following logical values:
%
%   * `true`    using Colpa's method (for details see [J.H.P. Colpa, Physica 93A (1978) 327](http://www.sciencedirect.com/science/article/pii/0378437178901607)),
%               the dynamical matrix is converted into another Hermitian
%               matrix, that will give the real eigenvalues.
%   * `false`   using the standard method (for details see [R.M. White, PR 139 (1965) A450](https://journals.aps.org/pr/abstract/10.1103/PhysRev.139.A450))
%               the non-Hermitian $\mathcal{g}\times \mathcal{H}$ matrix
%               will be diagonalised, which is computationally less
%               efficient. Default value is `true`.
%
% {{note Always use Colpa's method, except when imaginary eigenvalues are
%   expected. In this case only White's method work. The solution in this
%   case is wrong, however by examining the eigenvalues it can give a hint
%   where the problem is.}}
%
% `'saveH'`
% : If true, the quadratic form of the Hamiltonian is also saved in the
%   output. Be carefull, it can take up lots of memory. Default value is
%   `false`.
%
% `'saveV'`
% : If true, the matrices that transform the normal magnon modes into the
%   magnon modes localized on the spins are also saved into the output. Be
%   carefull, it can take up lots of memory. Default value is `false`.
%
% `'saveSabp'`
% : If true, the dynamical structure factor in the rotating frame
%   $S'(k,\omega)$ is saved. Default value is `false`.
%
% `'title'`
% : Gives a title string to the simulation that is saved in the output.
%
% `'fid'`
% : Defines whether to provide text output. The default value is determined
%   by the `fid` preference stored in [swpref]. The possible values are:
%   * `0`   No text output is generated.
%   * `1`   Text output in the MATLAB Command Window.
%   * `fid` File ID provided by the `fopen` command, the output is written
%           into the opened file stream.
%
% `'tid'`
% : Determines if the elapsed and required time for the calculation is
%   displayed. The default value is determined by the `tid` preference
%   stored in [swpref]. The following values are allowed (for more details
%   see [sw_timeit]):
%   * `0` No timing is executed.
%   * `1` Display the timing in the Command Window.
%   * `2` Show the timing in a separat pup-up window.
%
% ### Output Arguments
%
% `spectra`
% : structure, with the following fields:
%   * `omega`   Calculated spin wave dispersion with dimensions of
%               $[n_{mode}\times n_{Q}]$.
%   * `Sab`     Dynamical structure factor with dimensins of
%               $[3\times 3\times n_{mode}\times n_{Q}]$. Each
%               `(:,:,i,j)` submatrix contains the 9 correlation functions
%               $S^{xx}$, $S^{xy}$, $S^{xz}$, etc. If given, magnetic form
%               factor is included. Intensity is in \\hbar units, normalized
%               to the crystallographic unit cell.
%   * `H`       Quadratic form of the Hamiltonian. Only saved if `saveH` is
%               true.
%   * `V`       Transformation matrix from the normal magnon modes to the
%               magnons localized on spins using the following:
%               $x_i = \sum_j V_{ij} \times x_j'$
%               Only saved if `saveV` is true.
%   * `Sabp`    Dynamical structure factor in the rotating frame,
%               dimensions are $[3\times 3\times n_{mode}\times n_{Q}]$,
%               but the number of modes are equal to twice the number of
%               magnetic atoms.
%   * `formfact`  Cell containing the labels of the magnetic ions if form
%               factor in included in the spin-spin correlation function.
%   * `cmplxBase` The local coordinate system on each magnetic moment is
%               defined by the complex magnetic moments:
%               $\begin{align}  e_1 &= \Im(\hat{M})\\
%                               e_3 &= Re(\hat{M})\\
%                               e_2 &= e_3\times e_1
%               \end{align}$
%
%   * `hkl`     Contains the input $Q$ values, dimensions are $[3\times n_{Q}]$.
%   * `hklA`    Same $Q$ values, but in $\\ang^{-1}$ unit, in the
%               lab coordinate system, dimensins are $[3\times n_{Q}]$.
%   * `incomm`  Logical value, tells whether the calculated spectra is
%               incommensurate or not.
%   * `obj`     The copy (clone) of the input `obj`, see [spinw.copy].
%
% The number of magnetic modes (labeled by `nMode`) for commensurate
% structures is double the number of magnetic atoms in the magnetic cell.
% For incommensurate structures this number is tripled due to the
% appearance of the $(Q\pm k_m)$ Fourier components in the correlation
% functions. For every $Q$ points in the following order:
% $(Q-k_m,Q,Q+k_m)$.
%
% If several twins exist in the sample, `omega` and `Sab` are packaged into
% a cell, that contains $n_{twin}$ number of matrices.
%
% ### See Also
%
% [spinw] \| [spinw.spinwavesym] \| [sw_mex] \| [spinw.powspec] \| [sortmode]
%

% $Name: SpinW$ ($Version: 3.1$)
% $Author: S. Tóth and S. Ward$ ($Contact: admin@spinw.org, @spinw4 on Twitter$)
% $Revision: 1591$ ($Date: 25-Apr-2019$)
% $License: GNU GENERAL PUBLIC LICENSE$

pref = swpref;

% for linear scans create the Q line(s)
if nargin > 1
    hkl = sw_qscan(hkl);
else
    hkl = [];
end

% save warning of eigorth
orthWarn0 = false;

% save warning for singular matrix
singWarn0 = warning('off','MATLAB:nearlySingularMatrix');

% calculate symbolic spectrum if obj is in symbolic mode
if obj.symbolic
    if numel(hkl) == 3
        hkl = sym(hkl);
    end
    
    if ~isa(hkl,'sym')
        inpForm.fname  = {'fitmode'};
        inpForm.defval = {false    };
        inpForm.size   = {[1 1]    };
        param0 = sw_readparam(inpForm, varargin{:});
        
        if ~param0.fitmode
            warning('spinw:spinwave:MissingInput','No hkl value was given, spin wave spectrum for general Q (h,k,l) will be calculated!');
        end
        spectra = obj.spinwavesym(varargin{:});
    else
        spectra = obj.spinwavesym(varargin{:},'hkl',hkl);
    end
    return
end

% help when executed without argument
if nargin==1
    swhelp spinw.spinwave
    spectra = [];
    return
end

title0 = 'Numerical LSWT spectrum';

inpForm.fname  = {'fitmode' 'notwin' 'sortMode' 'optmem' 'tol' 'hermit'};
inpForm.defval = {false     false    true       0        1e-4  true    };
inpForm.size   = {[1 1]     [1 1]    [1 1]      [1 1]    [1 1] [1 1]   };

inpForm.fname  = [inpForm.fname  {'omega_tol' 'saveSabp' 'saveV' 'saveH'}];
inpForm.defval = [inpForm.defval {1e-5        false      false   false  }];
inpForm.size   = [inpForm.size   {[1 1]       [1 1]      [1 1]   [1 1]  }];

inpForm.fname  = [inpForm.fname  {'formfact' 'formfactfun' 'title' 'gtensor'}];
inpForm.defval = [inpForm.defval {false       @sw_mff      title0  false    }];
inpForm.size   = [inpForm.size   {[1 -1]      [1 1]        [1 -2]  [1 1]    }];

inpForm.fname  = [inpForm.fname  {'cmplxBase' 'tid' 'fid' }];
inpForm.defval = [inpForm.defval {false       -1    -1    }];
inpForm.size   = [inpForm.size   {[1 1]       [1 1] [1 1] }];

param = sw_readparam(inpForm, varargin{:});

if ~param.fitmode
    % save the time of the beginning of the calculation
    spectra.datestart = datestr(now);
end

if param.fitmode
    param.sortMode = false;
    param.tid = 0;
end

if param.tid == -1
    param.tid = pref.tid;
end

if param.fid == -1
    param.fid = pref.fid;
end
fid = param.fid;

% generate magnetic structure in the rotating noation
magStr = obj.magstr;

% size of the extended magnetic unit cell
nExt    = magStr.N_ext;
% magnetic ordering wavevector in the extended magnetic unit cell
km = magStr.k.*nExt;

% whether the structure is incommensurate
incomm = any(abs(km-round(km)) > param.tol);

% Transform the momentum values to the new lattice coordinate system
hkl = obj.unit.qmat*hkl;

% Calculates momentum transfer in A^-1 units.
hklA = 2*pi*(hkl'/obj.basisvector)';

% Check for 2*km
tol = param.tol*2;
helical =  sum(abs(mod(abs(2*km)+tol,1)-tol).^2) > tol;

% number of Q points
nHkl0 = size(hkl,2);

% define Q scans for the twins
nTwin = size(obj.twin.vol,2);
if param.notwin
    nTwin = 1;
end

% if the single twin has no rotation set param.notwin true
rotc1 = obj.twin.rotc(:,:,1)-eye(3);
if (nTwin == 1) && norm(rotc1(:))==0
    param.notwin = true;
end

if ~param.notwin
    % In the abc coordinate system of the selected twin the scan is
    % rotated opposite direction to rotC.
    hkl  = obj.twinq(hkl);
    nHkl = nHkl0*nTwin;
else
    nHkl = nHkl0;
    hkl  = {hkl};
end

if incomm
    % TODO
    if ~helical && ~param.fitmode
        warning('spinw:spinwave:Twokm',['The two times the magnetic ordering '...
            'wavevector 2*km = G, reciproc lattice vector, use magnetic supercell to calculate spectrum!']);
    end
    
    hkl0 = cell(1,nTwin);
    hklExt = cell(1,nTwin);
    
    for tt = 1:nTwin
        % without the k_m: (k, k, k)
        hkl0{tt} = repmat(hkl{tt},[1 3]);
        
        % for wavevectors in the extended unit cell km won't be multiplied by
        % nExt (we devide here to cancel the multiplication later)
        kme = km./nExt;
        hklExt{tt}  = [bsxfun(@minus,hkl{tt},kme') hkl{tt} bsxfun(@plus,hkl{tt},kme')];
        
        % calculate dispersion for (k-km, k, k+km)
        hkl{tt}  = [bsxfun(@minus,hkl{tt},km') hkl{tt} bsxfun(@plus,hkl{tt},km')];
    end
    nHkl0 = nHkl0*3;
else
    hkl0   = hkl;
    hklExt = hkl;
end

hkl    = cell2mat(hkl);
hkl0   = cell2mat(hkl0);
hklExt = cell2mat(hklExt);

% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
[~, SI, RR] = obj.intmatrix('fitmode',true,'conjugate',true);


hklExt0 = bsxfun(@times,hkl0,nExt')*2*pi;

M0 = magStr.S;
S0 = sqrt(sum(M0.^2,1));
% normal to rotation of the magnetic moments
n  = magStr.n;
nMagExt = size(M0,2);


if ~param.cmplxBase
    if obj.symbolic
        e3 = simplify(M0./[S0; S0; S0]);
        % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
        e2  = [zeros(1,nMagExt); e3(3,:); -e3(2,:)];
        % select zero vector and make them parallel to [0,0,1]
        selidx = abs(e2)>0;
        if isa(selidx,'sym')
            e2(3,~any(~sw_always(abs(e2)==0))) = 1;
        else
            e2(3,~any(abs(e2)>0)) = 1;
        end
        E0 = sqrt(sum(e2.^2,1));
        e2  = simplify(e2./[E0; E0; E0]);
        % e1 = e2 x e3
        e1  = simplify(cross(e2,e3));
    else
        % e3 || Si
        e3 = bsxfun(@rdivide,M0,S0);
        % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
        e2  = [zeros(1,nMagExt); e3(3,:); -e3(2,:)];
        e2(3,~any(abs(e2)>1e-10)) = 1;
        e2  = bsxfun(@rdivide,e2,sqrt(sum(e2.^2,1)));
        % e1 = e2 x e3
        e1  = cross(e2,e3);
    end
else
    F0  = obj.mag_str.F;
    RF0 = sqrt(sum(real(F0).^2,1));
    IF0 = sqrt(sum(imag(F0).^2,1));
    % e3 = real(M)
    e3  = real(F0)./repmat(RF0,[3 1]);
    % e1 = imag(M) perpendicular to e3
    e1  = imag(F0)./repmat(IF0,[3 1]);
    e1  = e1-bsxfun(@times,sum(e1.*e3,1),e3);
    e1  = e1./repmat(sqrt(sum(e1.^2,1)),[3 1]);
    % e2 = cross(e3,e1)
    e2  = cross(e3,e1);
end
% assign complex vectors that define the rotating coordinate system on
% every magnetic atom
zed = e1 + 1i*e2;

if param.gtensor
    
    gtensor = SI.g;
    
    if incomm
        % keep the rotation invariant part of g-tensor
        nx  = [0 -n(3) n(2);n(3) 0 -n(1);-n(2) n(1) 0];
        nxn = n'*n;
        m1  = eye(3);
        gtensor = 1/2*gtensor - 1/2*mmat(mmat(nx,gtensor),nx) + 1/2*mmat(mmat(nxn-m1,gtensor),nxn) + 1/2*mmat(mmat(nxn,gtensor),2*nxn-m1);
    end
end

% empty Sab
Sab = zeros(3,3,2*nMagExt,0);

warn1 = false;

% calculate all magnetic form factors
if param.formfact
    spectra.formfact = true;
    % Angstrom^-1 units for Q
    hklA0 = 2*pi*(hkl0'/obj.basisvector)';
    % store form factor per Q point for each atom in the magnetic supercell
    % TODO check prod(nExt)? instead of nExt
    %FF = repmat(param.formfactfun(permute(obj.unit_cell.ff(1,:,obj.matom.idx),[3 2 1]),hklA0),[1 nExt]);
    FF = repmat(param.formfactfun(permute(obj.unit_cell.ff(1,:,obj.matom.idx),[3 2 1]),hklA0),[prod(nExt) 1]);
else
    spectra.formfact = false;
end

nHklMEM = size(hklExt,2);

% Calculates correlation functions.
% V right
VExtR = repmat(permute(V  ,[4 5 1 2 3]),[3 3 1 1 1]);
% V left: conjugate transpose of V
VExtL = conj(permute(VExtR,[1 2 4 3 5]));

% Introduces the exp(-ikR) exponential factor.
ExpF =  exp(-1i*sum(repmat(permute(hklExt0,[1 3 2]),[1 nMagExt 1]).*repmat(RR,[1 1 nHklMEM]),1));
% Includes the sqrt(Si/2) prefactor.
ExpF = ExpF.*repmat(sqrt(S0/2),[1 1 nHklMEM]);

ExpFL =      repmat(permute(ExpF,[1 4 5 2 3]),[3 3 2*nMagExt 2]);
% conj transpose of ExpFL
ExpFR = conj(permute(ExpFL,[1 2 4 3 5]));

zeda = repmat(permute([zed conj(zed)],[1 3 4 2]),[1 3 2*nMagExt 1 nHklMEM]);
% conj transpose of zeda
zedb = conj(permute(zeda,[2 1 4 3 5]));

% calculate magnetic structure factor using the hklExt0 Q-values
% since the S(Q+/-k,omega) correlation functions also belong to the
% F(Q)^2 form factor

if param.formfact
    % include the form factor in the z^alpha, z^beta matrices
    zeda = zeda.*repmat(permute(FF(:,hklIdxMEM),[3 4 5 1 2]),[3 3 2*nMagExt 2 1]);
    zedb = zedb.*repmat(permute(FF(:,hklIdxMEM),[3 4 1 5 2]),[3 3 2 2*nMagExt 1]);
end

if param.gtensor
    % include the g-tensor
    zeda = mmat(repmat(permute(gtensor,[1 2 4 3]),[1 1 1 2]),zeda);
    zedb = mmat(zedb,repmat(gtensor,[1 1 2]));
end
% Dynamical structure factor from S^alpha^beta(k) correlation function.
% Sab(alpha,beta,iMode,iHkl), size: 3 x 3 x 2*nMagExt x nHkl.
% Normalizes the intensity to single unit cell.
Sab = cat(4,Sab,squeeze(sum(zeda.*ExpFL.*VExtL,4)).*squeeze(sum(zedb.*ExpFR.*VExtR,3))/prod(nExt));

if param.sortMode
    % sort the spin wave modes
    [omega, Sab] = sortmode(omega,reshape(Sab,9,size(Sab,3),[]));
    Sab          = reshape(Sab,3,3,size(Sab,2),[]);
end

[~,singWarn] = lastwarn;
% restore warning for singular matrix
warning(singWarn0.state,'MATLAB:nearlySingularMatrix');

% If number of formula units are given per cell normalize to formula
% unit
if obj.unit.nformula > 0
    Sab = Sab/double(obj.unit.nformula);
end

sw_timeit(100,2,param.tid);

fprintf0(fid,'Calculation finished.\n');

if warn1 && ~param.fitmode
    warning('spinw:spinwave:NonPosDefHamiltonian',['To make the Hamiltonian '...
        'positive definite, a small omega_tol value was added to its diagonal!'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END MEMORY MANAGEMENT LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if incomm
    % resize matrices due to the incommensurability (k-km,k,k+km) multiplicity
    kmIdx = repmat(sort(repmat([1 2 3],1,nHkl0/3)),1,nTwin);
    % Rodrigues' rotation formula.
    nx  = [0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];
    nxn = n'*n;
    K1 = 1/2*(eye(3) - nxn - 1i*nx);
    K2 = nxn;
    
    % keep the rotation invariant part of Sab
    %nx  = [0 -n(3) n(2);n(3) 0 -n(1);-n(2) n(1) 0];
    %nxn = n'*n;
    m1  = eye(3);
    
    % if the 2*km vector is integer, the magnetic structure is not a true
    % helix
    %tol = param.tol*2;
    %helical =  sum(abs(mod(abs(2*km)+tol,1)-tol).^2) > tol;
    
    if helical
        % integrating out the arbitrary initial phase of the helix
        Sab = 1/2*Sab - 1/2*mmat(mmat(nx,Sab),nx) + 1/2*mmat(mmat(nxn-m1,Sab),nxn) + 1/2*mmat(mmat(nxn,Sab),2*nxn-m1);
    end
    
    % Save the structure factor in the rotating frame
    if param.saveSabp
        Sabp = Sab(:,:,:,kmIdx==2);
        omegap = omega(:,kmIdx==2);
    end
    
    % dispersion
    omega = [omega(:,kmIdx==1); omega(:,kmIdx==2); omega(:,kmIdx==3)];
    % exchange matrices
    Sab   = cat(3,mmat(Sab(:,:,:,kmIdx==1),K1), mmat(Sab(:,:,:,kmIdx==2),K2), ...
        mmat(Sab(:,:,:,kmIdx==3),conj(K1)));
    
    hkl   = hkl(:,kmIdx==2);
    nHkl0 = nHkl0/3;
else
    helical = false;
end

if ~param.notwin
    % Rotate the calculated correlation function into the twin coordinate
    % system using rotC
    SabAll = cell(1,nTwin);
    for ii = 1:nTwin
        % select the ii-th twin from the Q points
        idx    = (1:nHkl0) + (ii-1)*nHkl0;
        % select correlation function of twin ii
        SabT   = Sab(:,:,:,idx);
        % size of the correlation function matrix
        sSabT  = size(SabT);
        % convert the matrix into cell of 3x3 matrices
        SabT   = reshape(SabT,3,3,[]);
        % select the rotation matrix of twin ii
        rotC   = obj.twin.rotc(:,:,ii);
        % rotate correlation function using arrayfun
        SabRot = arrayfun(@(idx)(rotC*SabT(:,:,idx)*(rotC')),1:size(SabT,3),'UniformOutput',false);
        SabRot = cat(3,SabRot{:});
        % resize back the correlation matrix
        SabAll{ii} = reshape(SabRot,sSabT);
    end
    Sab = SabAll;
    
    if nTwin == 1
        Sab = Sab{1};
    else
        omega = mat2cell(omega,size(omega,1),repmat(nHkl0,[1 nTwin]));
    end
    
end

% Creates output structure with the calculated values.
spectra.omega    = omega;
spectra.Sab      = Sab;
spectra.hkl      = obj.unit.qmat\hkl(:,1:nHkl0);
spectra.hklA     = hklA;
spectra.incomm   = incomm;
spectra.helical  = helical;
spectra.norm     = false;
spectra.nformula = double(obj.unit.nformula);

% Save different intermediate results.
if param.saveV
    spectra.V = Vsave;
end
if param.saveH
    spectra.H = Hsave;
end
if param.saveSabp && incomm
    spectra.Sabp = Sabp;
    spectra.omegap = omegap;
end

% save the important parameters
spectra.param.notwin    = param.notwin;
spectra.param.sortMode  = param.sortMode;
spectra.param.tol       = param.tol;
spectra.param.omega_tol = param.omega_tol;
spectra.param.hermit    = param.hermit;
spectra.title           = param.title;
spectra.gtensor         = param.gtensor;

if ~param.fitmode
    spectra.dateend = datestr(now);
    spectra.obj = copy(obj);
end

if ~param.gtensor && any(obj.single_ion.g)
    warning('spinw:spinwave:NonZerogTensor',['The SpinW model defines a '...
        'g-tensor that is not included in the calculation. Anisotropic '...
        'g-tensor values cannot be applied afterwards as they change relative'...
        'spin wave intensities!'])
end

% issue eigorth warning
if orthWarn0
    warning('spinw:spinwave:NoOrth','Eigenvectors of defective eigenvalues cannot be orthogonalised at some q-point!');
end

if strcmp(singWarn,'MATLAB:nearlySingularMatrix')
    lineLink = 'line 846';
    if feature('HotLinks')
        lineLink = ['<a href="matlab:opentoline([''' sw_rootdir 'swfiles' filesep '@spinw' filesep 'spinwave.m''' '],846,0)">' lineLink '</a>'];
    end
    warning('spinw:spinwave:nearlySingularMatrix',['Matrix is close '...
        'to singular or badly scaled. Results may be inaccurate.\n> In spinw/spinwave (' lineLink ')']);
    %fprintf(repmat('\b',[1 30]));
end

end