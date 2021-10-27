classdef Optimizer
    properties
        wavelength % Wavelength of incident light:
        period % Period of the device
        thickness % Device layer thickness
        angleInc % Incident angle in degrees
        targetOrder % Target diffraction order
    end
    properties (Dependent = true)        
        polarization % Polarization: can be 'TM', 'TE' or 'Both'
        tgtPhiCoeffs 
    end
    properties
        indexData 
        substrateIndex 
        superstrateIndex 
        gapIndex 
        symmetry 
        simGridPtsPerWL 
    end
    
    properties (Dependent = true, Hidden = true)
        nwl % number of wavelength values
        numPol % number of polarizations: 1 or 2
    end
    
    properties (GetAccess = protected, SetAccess = protected, Hidden = true)
        privatePolarization 
        privateTgtPhiCoeffs 
    end
    
    methods
        function self = rmres(self)
        end
        
        function tf = eq(self, another, varargin)
            tf = true;
            pr = properties(self);
            for i = 1:length(pr)
                tf = tf && self.(pr{i}) == another.(pr{i});
            end
        end
        
        function tf = ne(self, another, varargin)
            tf = ~eq(self, another, varargin{:});
        end
        
        function self = set(self, varargin)
            if numel(self) > 1
                self = arrayfun(@(x) set(x, varargin{:}), self);
                return
            end
            fields_given = varargin(1:2:end);
            vals_given = varargin(2:2:end);
            assert(length(fields_given) == length(vals_given), 'Unequal number of Name/Value arguments')
            for ii = 1:length(fields_given)
                self.(fields_given{ii}) = vals_given{ii};
            end
        end
        
        function p = gen_patt(~, p)
        end
        
        function numpol = get.numPol(self)
            if strcmpi(self.polarization, 'both')
                numpol = 2;
            else
                numpol = 1;
            end
        end
        
        function strct = optimizer2struct(self, varargin)
            warning('off', 'MATLAB:structOnObject')
            strct = arrayfun(@struct, self);
        end
        
        function varargout = imsize(self, n)
        % returns size of the pattern
            sz = ceil(round(self.simGridPtsPerWL * self.period / min(self.wavelength),1));
            if isscalar(sz), sz = [sz 1]; end
            if nargin > 1
                sz = [sz, ones(1, max(n)-length(sz))];
                sz = sz(n);
            end
            if nargout  < 2
                varargout = {sz};
            else
                varargout = num2cell(sz);
            end
        end
        
        function pol = npol(self, n)
            switch lower(self.polarization)
                case 'te'
                    pol = 1;
                case 'tm'
                    pol = 2;
                case 'both'
                    pol = 1:2;
            end
            if nargin == 2
                pol = pol(n);
            end
        end
        
        function n = get.nwl(self)
            n = length(self.wavelength);
        end
        
        function plot(~)
        end
        
        function val = get.polarization(self)
            val = self.privatePolarization;
        end
        
        function self = set.polarization(self, val)
            self.privatePolarization = val;
        end
        
        function val = get.tgtPhiCoeffs(self)
            val = self.privateTgtPhiCoeffs;
        end
        
        function self = set.tgtPhiCoeffs(self, val)
            self.privateTgtPhiCoeffs = val;
        end
        
    end
    
    methods (Static)
        
        function s = substituteShortNames(s)
        end
        
    end
end