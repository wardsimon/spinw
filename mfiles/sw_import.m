function obj = sw_import(fName, toPlot, obj0)
% create SpinW object from .cif and FullProf Studio .fst files
%
% At present the function can import Crystallographic Information Framework
% (.cif) files or FullProf Studio (.fst) files.
%
% obj = SW_IMPORT(fName, {toPlot})
%
% Input:
%
% fName     String, contain the file location and name of the .fst file.
% toPlot    If true the structure will be plotted, default is false.
%

if nargin == 0
    help sw_import
    return
end

if nargin < 2
    toPlot = false;
end

if nargin < 3
    obj0 = spinw;
end

[~,~,fExt] = fileparts(fName);

switch fExt
    case '.fst'
        % read the file into a string
        str = fileread(fName);
        % split string up into lines
        
        str = strsplit(str,'\n');
        % remove empty lines
        str = str(~cellfun(@(C)isempty(C),str));
        % store all imported data
        dat = struct('atom',struct('label',{},'r',{}),'cell',[],'box',[],'matom',struct('label',{},'r',{},'M',{}),'k',zeros(3,0));
        
        for ii = 1:numel(str)
            sstr = str{ii};
            
            if ismember(sstr(1),'!{}')
                % !COMMENT and grouping
                continue
            end
            
            sSplit = strsplit(sstr,' ');
            switch sSplit{1}
                case 'SPACEG'
                    dat.spgr = sstr(8:end);
                case 'CELL'
                    dat.cell = str2double(sSplit(2:7));
                case 'BOX'
                    dat.box = str2double(sSplit(2:7));
                case 'ATOM'
                    dat.atom(end+1).label = [sSplit{2} ' ' sSplit{3}];
                    dat.atom(end).r       = str2double(sSplit(4:6))';
                case 'K'
                    dat.k = [dat.k str2double(sSplit(2:4))'];
                case 'BKG'
                    % [RGBT]
                    dat.bkg = str2double(sSplit(2:5))*255;
                case 'MATOM'
                    if numel(sSplit{3})>1
                        sSplit{3}(2) = lower(sSplit{3}(2));
                    end
                    dat.matom(end+1).label = [sSplit{2} ' ' sSplit{3}];
                    dat.matom(end).r       = str2double(sSplit(4:6))';
                case 'SKP'
                    % M: 3 x 1 x nK
                    kIdx = str2double(sSplit{2});
                    Mk   = 0.5*(str2double(sSplit(4:6))+1i*str2double(sSplit(7:9)))*exp(-2*pi*1i*str2double(sSplit{10}));
                    % opposite phi
                    %Mk   = str2double(sSplit(4:6))+1i*str2double(sSplit(7:9))*exp(2*pi*1i*str2double(sSplit{10}));
                    dat.matom(end).M(:,1,kIdx) = transpose(Mk);
            end
        end
        
        % create the SpinW model
        obj0.genlattice('lat_const',dat.cell(1:3),'angled',dat.cell(4:6),'spgr',dat.spgr);
        if ~isempty(dat.atom)
            obj0.addatom('r',[dat.atom(:).r],'label',{dat.atom(:).label});
        elseif ~isempty(dat.matom)
            % add fictious spin to the atom
            obj0.addatom('r',[dat.matom(:).r],'label',{dat.matom(:).label},'S',ones(1,numel(dat.matom)));
        end
        
        if ~isempty(dat.k)
            
            % generate the single-k magnetic structure
            obj0.mag_str.k = dat.k;
            
            % determine real magnetic moment components [3 x nMagAtom x nKm]
            M0 = [dat.matom(:).M];
            
            % add extra phase due to km*r
            %M0 = bsxfun(@times,M0,exp(sum(bsxfun(@times,dat.matom.r,1i*2*pi*permute(dat.k,[1 3 2])),1)));
            
            
            % number of k-vectors
            nK = size(dat.k,2);
            % convert from the lattice component coordinate system to Descartes
            % coordinate system
            for ii = 1:nK
                M0(:,:,ii) = obj0.basisvector(1)*M0(:,:,ii);
            end
            
            % keep only the first wave vector
            % the given magnetic moments are in Bohr magneton, SpinW assumes g=2, thus
            % the moments are divided by this number
            g0 = 1;
            % calculate the real part of the magnetic structure
            obj0.mag_str.F = M0/g0;
        end
        
    case '.cif'
        
        cif0 = cif(fName);
        abc0 = [cif0.cell_length_a cif0.cell_length_b cif0.cell_length_c];
        ang0 = [cif0.cell_angle_alpha cif0.cell_angle_beta cif0.cell_angle_gamma];
        sym0 = cif0.('symmetry_space_group_name_H-M');
        %symi0 = cif0.symmetry_Int_Tables_number;
        
        if ~isempty(cif0.symmetry_equiv_pos_as_xyz)
            xyz0 = cif0.symmetry_equiv_pos_as_xyz;
        elseif ~isempty(cif0.space_group_symop_operation_xyz)
            xyz0 = cif0.space_group_symop_operation_xyz;
        else
            warning('spinw:WrongFormat','Missing symmetry operators, using P1')
            xyz0 = {'x,y,z'};
        end
        
        xyz0 = sprintf('%s; ',xyz0{:}); xyz0 = xyz0(1:end-2);
        %name0 = cif0.atom_site_type_symbol';
        cell0 = [cif0.atom_site_label cif0.atom_site_type_symbol];
        name0 = cellfun(@(x,y)strjoin({x y}),cell0(:,1),cell0(:,2),'UniformOutput',false)';
        r0    = mod([cif0.atom_site_fract_x cif0.atom_site_fract_y cif0.atom_site_fract_z]',1);
        
        % save formula units
        if ~isempty(cif0.cell_formula_units_Z)
            obj0.unit.nformula = int32(cif0.cell_formula_units_Z);
        end
        
        if numel(abc0)==3
            obj0.lattice.lat_const = abc0;
        end
        if numel(ang0) == 3
            obj0.lattice.angle = ang0*pi/180;
        end
        if numel(xyz0) > 3
            % determine the symmetry generators
            [symOp, symTr] = sw_gensym(sym0, xyz0);
            [symOp, symTr] = sw_symgetgen(symOp, symTr);
            % save generators into spinw pbject
            obj0.lattice.label = sym0;
            obj0.lattice.sym   = [symOp permute(symTr,[1 3 2])];
        end
        
        if size(name0,2) == size(r0,2)
            % add atoms to the crystal structure
            obj0.addatom('r',r0,'label',name0,'occ',cif0.atom_site_occupancy','biso',cif0.atom_site_B_iso_or_equiv)
        else
            error('spinw:WrongInput','The .cif file contains inconsistent information!')
        end
end

fprintf0(obj0.fileid,'Crystal structure is imported from %s.\n',fName);

if nargout > 0
    obj = obj0;
end

if toPlot
    plot(obj0,'range',reshape(dat.box',[2 3])')
end

end