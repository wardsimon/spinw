classdef brille < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lattice = struct('dlat',[], 'rlat', [])
        brillouinzone
        interpolant
        max_volume = 1E-4
    end
    properties(Hidden=true)
       extent = [1, 1, 1]
       spinwObj
    end
    
    methods
        function obj = brille(varargin)
            
            if nargin == 0
                return
            else
                obj.spinwObj = varargin{1};
            end
            
            obj.extent = obj.spinwObj.magstr.N_ext;
            
            obj.createLattice(obj.spinwObj.lattice.lat_const.*obj.extent,...
                obj.spinwObj.lattice.angle,...
                obj.spinwObj.lattice.label)
        end
        
        function createLattice(obj, alatt, angdeg, spgrp)
            pylens = m2p( alatt(1:3) );
            pyangs = m2p( angdeg(1:3) );
            obj.lattice.dlat = py.brille.Direct(pylens, pyangs, spgrp);
            obj.lattice.rlat = obj.lattice.dlat.star;
            % Check for the pesky P0 symmetry
            if (p2m(obj.lattice.rlat.hall) == 0)
                warning('brille:spinwaveBrille', 'The symmetry is P0, falling back to P1')
                obj.lattice.rlat.hall = uint16(1); % Force P1 symmetry
            end
            obj.createBrillouinZone()
            obj.createInterpolant()
        end
        
        function updateLattice(obj, varargin)
            if ~isempty(varargin)
                obj.createLattice(varargin{:})
            else
                % This is using the built in spinW object.
                if isempty(obj.spinwObj)
                    error('brille:NoSpinW', 'The brille object doesnt have a spinW object')
                end
                obj.createLattice(obj.spinwObj.lattice.lat_const.*obj.extent,...
                obj.spinwObj.lattice.angle,...
                obj.spinwObj.lattice.label)
            end
            obj.createBrillouinZone()
            obj.createInterpolant()
        end
        
        function points = getBrillouinZoneRLU(obj)
            points = p2m(obj.interpolant.rlu);
        end
        
        function setBrillouinZoneEigs(obj, V, type)
            obj.interpolant.fill(m2p(V), int32(type))
        end
        
        function V = queryInterpolant(obj, hkl)
            V = p2m(obj.interpolant.ir_interpolate_at(m2p(hkl')));
        end
        
        function fillFromSpinW(obj)
            % Calculate the eigen values/vectors of the Brillouin zone points
            [omega, V] = obj.spinwObj.spinwaveEigVecs(obj.getBrillouinZoneRLU()');

            % Reshape into what Brille requires
            VV = cat(3, permute(omega, [2, 1, 3]), permute(V, [3, 1, 2]));
            if incomm
                VV = reshape(VV, [], 3, size(VV, 2), size(VV, 3));
            end
            size_VV = size(VV);
            num_el = prod(size_VV(3:end));
            obj.setBrillouinZoneEigs(VV, [num_el,0,0]);
        end
        
        function spectra = spinwave(obj, hkl, varargin)
            obj.fillFromSpinW()
            % Call the interpolation object with the required HKL
            hkl = sw_qscan(hkl);
            VVnew = p2m(obj.interpolant.ir_interpolate_at(m2p(hkl')));

            % Do the inverse of Brille formatting to get spinW formated eigen
            % values/vectors
            if incomm
                VVnew = reshape(VVnew, [], size(VV, 3), size(VV, 4));
            end
            Vnew = permute(VVnew(:, :, 2:end), [2, 3, 1]);
            omegaNew = squeeze(VVnew(:, :, 1))';

            % Use the Brille eigen values/vectors to create a spectra object.
            spectra = obj.spinwObj.spinwave_structureFact(obj, hkl, omegaNew, Vnew, varargin{:});    
        end
    end
    
    methods(Hidden=true)
        
        function createBrillouinZone(obj)
            obj.brillouinzone = py.brille.BrillouinZone(obj.lattice.rlat);
        end
        
        function createInterpolant(obj)
            obj.interpolant = brille.BZTrellisQ(obj.brillouinzone, ...
                'max_volume', obj.max_volume, ...
                'complex', true);
        end
    end
end

