function gencoupling2(obj, varargin)
% generates the COUPLING property of sw object
%
% GENCOUPLING(obj, 'option1', value, ...)
%
% The function calculates equivalent bonds between magnetic atoms. These
% are determined either based on crystal symmetry or bond length (with
% tolDist tolerance). If the space group index of 0 is defined
% (obj.lattice.sym=0), the equivalent bonds will be based on bond length.
% For space group index larger than 0, the symmetry equivalent bonds will
% be determined. This can ve overwritten by the forceNoSym parameter to
% consider bond length.
%
% IMPORTANT!
%   This function has to be used after the crystal structure is defined.
%   The SW.ADDCOUPLING, SW.COUPLINGTABLE functions will only work
%   afterwards.
%
% Input:
%
% obj           spinw class object.
%
% Options:
%
% forceNoSym    If true, equivalent bonds are generated based on
%               bond length with .tolDist tolerance. If false symmetry
%               operators will be used if they are given
%               (obj.lattice.sym>0).
% maxDistance   Maximum bond length that will be stored in the
%               obj.coupling property in units of Angstrom. Default is 8.
% maxSym        Maximum bond length until the symmetry equivalent bonds are
%               generated. It is usefull if long bonds have to be generated
%               for the dipolar interaction, but the symmetry analysis of
%               them is not necessary. Default value is equal to
%               maxDistance.
% tolDist       Tolerance of distance, within two bonds are regarded
%               equivalent, default is 1e-3 Angstrom. Only used, when no
%               space group is defined.
% dMin          Minimum bond length, below which an error is triggered.
%               Default value is 0.5 Angstrom.
%
% Output:
%
% The obj.coupling field will be filled with values, depending on the
% crystal geometry.
%
% Example:
%
% cryst = sw;
% cryst.genlattice('lat_const',[3 3 5],'angled',[90 90 120])
% cryst.addatom('r',[0 0 0])
% cryst.gencoupling
% cryst.couplingtable(1:3)
%
% A triangular lattice is created in cryst and after using gencoupling()
% the couplingtable() function lists the 1st, 2nd and 3rd neighbor bonds.
%
% See also SPINW, SPINW.SYMMETRY, SPINW.NOSYM.
%

isSym = obj.lattice.sym > 0;

inpForm.fname  = {'forceNoSym' 'maxDistance' 'tol' 'tolDist' 'dMin' 'maxSym'};
inpForm.defval = {false         8             1e-5  1e-3      0.5   []      };
inpForm.size   = {[1 1]        [1 1]         [1 1] [1 1]     [1 1]  [1 1]   };
inpForm.soft   = {false        false          false false    false  true    };


param = sw_readparam(inpForm, varargin{:});
tol   = param.tol;
tolD  = param.tolDist;

if isempty(param.maxSym)
    param.maxSym = param.maxDistance;
end

% to avoid some problem with symmetry if maxDistance equal to a lattice
% constant
param.maxDistance = param.maxDistance + tol;

fid = obj.fileid;

% calculate the height of the paralellepiped of a unit cell
hMax = 1./sqrt(sum(inv(obj.basisvector).^2,2));

% gives the number of extra unit cells along all 3 axes
% that are necessary to cover the minimum bond distance
nC = ceil(param.maxDistance./hMax);

% force no symmetry operator mode
if param.forceNoSym
    isSym = false;
end

if fid ~= 0
    fprintf0(fid,['Creating the bond list (maxDistance = %g ' symbol('ang')...
        ', nCell = %dx%dx%d) ...\n'],param.maxDistance-tol,nC);
end

% save the sym/nosym method into obj
obj.sym = isSym;

% Generate all atomic positions in the unit cell.
mAtom = obj.matom;

% Number of magnetic atoms in the unit cell.
nMagAtom = size(mAtom.r,2);

if nMagAtom == 0
    error('sw:gencoupling:NoMagAtom','There is no magnetic atom in the unit cell with S>0!');
end

% Use half 'cube' around the center unit cell and remove the identical
% bonds using unique() to simplify the calculation

% number of unit cells in the half 'cube'
cDim = [nC(1)+1 2*nC(2)+1 2*nC(3)+1];
nHalfCube = prod(cDim);

% generate all cell translations
cTr{1} = 0:nC(1);
cTr{2} = -nC(2):nC(2);
cTr{3} = -nC(3):nC(3);

% generate cell origin positions in l.u.
[cTr{1}, cTr{2}, cTr{3}] = ndgrid(cTr{:});
% cell origin translations: Na x Nb x Nc x 1 x 1 x3
cTr = cat(6,cTr{:});
% remove unit cells that would produce duplicate bonds
% mark them with NaN (enough to do along a-axis values)
cTr(1,:,(1:end)<=nC(3),1,1,1)       = nan;
cTr(1,(1:end)<=nC(2),nC(3)+1,1,1,1) = nan;
% positions of mgnetic atoms in the half 'cube' in l.u.
% Na x Nb x Nc x 1 x nMagAtom x 3
% atom2
aPos = bsxfun(@plus,permute(mAtom.r,[3:6 2 1]),cTr);
% generate all distances from the atoms in the (0,0,0) cell in l.u.
% r(atom2) - r(atom1)
% Na x Nb x Nc x nMagAtom x nMagAtom x 3
dR  = bsxfun(@minus,aPos,permute(mAtom.r,[3:5 2 6 1]));
% mark duplicate bonds bonds within the (0,0,0) cell with nan
R0 = permute(dR(1,nC(2)+1,nC(3)+1,:,:,:),[4:6 1:3]);
dR(1,nC(2)+1,nC(3)+1,:,:,:) = bsxfun(@times,tril(nan(nMagAtom))+1,R0);
% calculate the absolute value of the distances in Angstrom
dRA = sqrt(sum(mmat(permute(obj.basisvector,[3:7 1 2]),dR,[6 7]).^2,6));
% reshape the numbers into a column list of bonds
% 3 x Na x Nb x Nc x nMagAtom x nMagAtom
dl    = permute(repmat(cTr,[1 1 1 nMagAtom nMagAtom 1]),[6 1:5]);
atom1 = repmat(permute(1:nMagAtom,[1 3:5 2 6]),[1 cDim 1 nMagAtom]);
atom2 = repmat(permute(1:nMagAtom,[1 3:5 6 2]),[1 cDim nMagAtom 1]);
dRA   = permute(dRA,[6 1:5]);
% store everything in a single matrix
cMat  = [dl(:,:);atom1(1,:);atom2(1,:);dRA(1,:)];
% remove nan bonds
cMat = cMat(:,~isnan(cMat(1,:)));
% sort according to increasing distance
% cMat = sortrows(cMat',6)'; too slow
[~,cIdx] = sort(cMat(6,:));
cMat     = cMat(:,cIdx);
% cutoff at maximum distance
cMat = cMat(:,cMat(6,:)<=param.maxDistance);
% index the bonds
cMat(7,:) = cumsum([1 diff(cMat(6,:))>tolD]);

% check whether some atoms are too close
if cMat(6,1) < param.dMin
    error('sw:gencoupling:AtomPos',['Some atoms are too close (Dmin=' ...
        num2str(cMat(6,1)) '<' num2str(param.dMin) '), check your crystal structure!']);
end

dRA  = cMat(6,:);
% keep the bond length
cMat = cMat([1:5 7],:);
% cMat rows: [la, lb, lc, atom1, atom2, idx]

% symmetry equivalent couplings
if isSym
    % get the symmetry operators
    [symOp, symTr] = sw_gencoord(obj.lattice.sym);
    % store the final sorted couoplings in nMat
    nMat = zeros(6,0);
    ii  = 1;
    idx = 1;
    % maximum bond index for symmetry operations
    maxidxSym = max(cMat(6,dRA<=param.maxSym));
    
    while ii <= maxidxSym
        % select columns from sorM with a certain idx value
        sortMs = cMat(:,cMat(6,:) == ii);
        while size(sortMs,2)>0
            [genC, unC] = sw_gensymcoupling(obj, sortMs(:,1), {symOp, symTr}, tol);
            genCAll = [genC [-genC(1:3,:);genC([5 4],:)]];
            % remove from sortMs the identical couplings
            iNew = isnew(genCAll(1:5,:),sortMs(1:5,:),tol);
            sortMs(:,~iNew) = [];
            % remove identical couplings from the symmetry generated
            % list
            genC(:,~unC) = [];
            if sum(~iNew) ~= sum(unC)
                error('sw:gencoupling:SymProblem','Symmetry error! ii=%d, idx=%d. Try to change ''tol'' parameter.',ii,idx);
            end
            % move the non-unique (not new) couplings (symmetry equivalent ones)
            nMat = [nMat [genC;ones(1,size(genC,2))*idx]]; %#ok<AGROW>
            idx  = idx + 1;
        end
        ii = ii + 1;
    end
    
    % include the increase of bond index in case bonds are splitted due to
    % symmetry inequivalency
    cMat = cMat(:,cMat(6,:)>maxidxSym);
    cMat(6,:) = cMat(6,:) + nMat(6,end)-maxidxSym;
    cMat = [nMat cMat];
end

% default anisotropy and g-tensor values
aniso = int32(zeros(1,nMagAtom));
g     = int32(zeros(1,nMagAtom));

if fid ~= 0
    fprintf0(fid,'... %d bonds are retained out of %d generated!\n',size(cMat,2),nHalfCube*nMagAtom^2);
end

% save output structure
coupling.dl      = int32(cMat(1:3,:));
coupling.atom1   = int32(cMat(4,:));
coupling.atom2   = int32(cMat(5,:));
coupling.idx     = int32(cMat(6,:));
coupling.mat_idx = int32(zeros(3,size(coupling.idx,2)));
coupling.type    = int32(zeros(3,size(coupling.idx,2)));
coupling.rdip    = obj.coupling.rdip;

obj.coupling         = coupling;
obj.single_ion.aniso = aniso;
obj.single_ion.g     = g;

validate(obj);

end

function [isnew, symIdx] = isnew(A,B, tol)
% selects the new vectors from B. Dimensions of A and B have to be [3 nA]
% and [3 nB] respectively.
%
% A vector in B is considered new, if d(mod(vA-vB,1))<tol.
%

nA = size(A,2);
nB = size(B,2);

notequal = sum(abs(repmat(permute(A,[2 3 1]),[1 nB 1]) - repmat(permute(B,[3 2 1]),[nA 1 1])).^2,3) > tol;

isnew = all(notequal,1);

idx = 1:nB;

symIdx = arrayfun(@(idx)find(~notequal(:,idx),1,'first'),idx(~isnew));

end