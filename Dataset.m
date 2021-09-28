classdef Dataset
    properties (SetAccess = public, GetAccess = public)
        phase % phase [length(wavelength) x polarization x N]
        absEff % absulute efficiency [length(wavelength) x polarization x N]
        relEff % relative efficiency [length(wavelength) x polarization x N]
        pattern % 1xN cell array of pattern geometries, generator functions 
        % that take a two integer (or a scalal) value as argument for the
        % size of the image
        object (1,1) Optimizer % Optimizer object for record. It's assumed (and verified) that 
               % devices in a given dataset are identical (at least the
               % properties that Optimizer's eq method compares)
    end

    properties (Dependent)
        wavelength % returns wavelength property from the Optimizer in the 'object' property
    end
    
    properties (Dependent, Hidden = true)
        property % returns all (even hidden and private), but not Dependent properties
    end
    
    properties (SetAccess = public, GetAccess = public, Hidden = false)
        tgtPhiCoeffs % 3xnumPol vector, where rows are delta2_phi, delta_phi, phi_zero
        phiCoeffs % for the output of getPhi
        start StartDataset % stores start patterns and their characteristics used for the optimizations
        hash % SHA-1 hash cell array of devices (some of the methods use it)
    end
    
    properties (Access = public)
        info % comment and information string
    end
    
    methods
        
        %------------------------------------------
        %%% Constructor
        %------------------------------------------
        
        function obj = Dataset(data)
        % data can be:
        %   - instance or array of class Optimizer
        %   - a struct with fields of the same names as Dataset properties
        %   arrays of struct give arrays of Dataset, while arrays of
        %   Optimizer are gathered into instances of Dataset according to
        %   their eq method
            if nargin > 0
                if isa(data, 'Optimizer')
                    if numel(data) > 1
                        test = eq(data(1), data, false);
                        % if there are objects that correspond to different
                        % geometries, group them together 
                        if sum(test) ~= numel(data)
                            obj = [Dataset(data(test)) Dataset(data(~test))];
                            return
                        end
                    end
                    obj.object = rmres(data(1));
                    nwl = length(data(1).wavelength);
                    numPol = data(1).numPol;
                    robustInd = data(1).robustStartDeviation == 0;

                    % avoid squeeze for the fear of dimensions of length 1 that we'd like
                    % to keep (what if nwl == 1 and numPol == 1)
                    fields = {'phase', 'absEff', 'relEff'};
                    for f = 1:length(fields)
                        val = cat(5, data.(fields{f}));
                        if ~isempty(val)
                            obj.(fields{f}) = reshape(val(end, :, :, robustInd, :), nwl, numPol, []);
                        end
                    end

                    obj.tgtPhiCoeffs = cat(3, data.tgtPhiCoeffs);
                    if size(obj.tgtPhiCoeffs,1) < 3 
                        % exact values given, return to order:
                        % [wavelength, polarization, device] and get
                        % the coefficients
                        tgtphase = permute(obj.tgtPhiCoeffs, [2 1 3]);
                        [~, obj.tgtPhiCoeffs] = getPhi(set(obj,'phase',tgtphase));
                        obj.tgtPhiCoeffs = obj.tgtPhiCoeffs(:,:,:,1);
                    end


                    % this if prevents infinite loops
                    if ~isa(obj, 'StartDataset') % isa(obj, 'Dataset') says yes for subclasses
                        obj.start = StartDataset(data);
                    end
                    obj = obj.getPhi;
                    obj.pattern = {data.finalPattern};
                    obj.hash = cellfun(@(x) DataHash(x, 'SHA-1'), obj.pattern, 'UniformOutput', false);
                else
                    switch class(data)
                        case 'struct'
                            obj = repmat(obj, size(data));
                            names = fieldnames(data);
                            for n = 1:length(names)
                                try
                                    [obj.(names{n})] = data.(names{n});
                                catch
                                    disp(['There is no ', names{n}, ' property. Skipping...'])
                                end
                            end
                            if length([obj.hash]) ~= length([obj.pattern])
                                for n = 1:numel(obj)
                                    if ~iscell(obj(n).pattern)
                                        obj(n).pattern = {obj(n).pattern};
                                    end
                                    obj(n).hash = arrayfun(@(x) DataHash(x{1}, 'SHA-1'), obj(n).pattern, 'UniformOutput', false);
                                end
                            end
                        case 'StartDataset'
                            % just copying all the fields that are preset in
                            % Dataset (subclasses can have additional ones)
                            obj = repmat(obj, size(data));
                            names = properties(Dataset);
                            for n = 1:length(names)
                                [obj.(names{n})] = data.(names{n});
                            end
                        case 'Dataset'
                            obj = data;
                        otherwise
                            error('Object of class %s is not convertible to Dataset', class(data))
                    end
                end
            end
        end
        
        %------------------------------------------
        %%% End constructor
        %------------------------------------------
        
        %------------------------------------------
        %%% Basic operations methods
        %------------------------------------------
        
        function self = append(self, data, force)
        % case 1: data - object, array, or cell arry of objects of class Optimizer
        % case 2: data - struct with (at least some of) the same fields as Dataset
        % case 3: data - object, array, or cell arry of objects of class Dataset
        % force == 1,-1,0 -> if wavelength array have different step, but
        % the same boundaries, interpolate (get missing values for the 
        % Dataset with less sampled points), downsample (remove values from 
        % Dataset with more sampled points), or ask, respecively.
        %       -2, 2 -> prepend/append to the array
        % In the end, number of wavelength points has to be the same, 
        % because otherwise, some properties would be matrices of different 
        % size in the first dimension and could be contactentated. 
        % If imag(force) ~= 0, i.e. force = -1i, 1i, we ignore the object 
        % comparisoin step. It's useful when we want to visualize devices of
        % various geometries (like periods, heights), but in the same band
        
            if nargin == 2
                force = 0;
            end
            if isa(data, 'cell') 
                data = [data{:}];
            end
            % convert convertable classes to Dataset
            if ~isa(data, 'Dataset')
                data = Dataset(data);
            end
            if numel(self) > 1
                for ii = 1:length(data)
                    % adding data one by one
                    self = append(data(ii), self, force); 
                end
                return
            end
            
            %%% At this stage self is considered a scalar object
            
            % In the case where data is an array, we choose the first
            % Dataset with the same object, append to it and return the
            % whole array, or look for potential interpolate candidates.
            % One could have this code work for verctor data, but it's
            % simpler like that and it is done so that we append self
            % only to one of the data element, instead of merging them
            if numel(data) > 1
                I = find(cellfun(@(x) eq(self.object, x, false), {data.object}), 1);
                if isempty(I) 
                    % check for potential interpolation candidates
                    I = cellfun(@(x) isequal(self.wavelength([1 end]), x([1 end])), {data.wavelength});
                    % make sure that other than the wavelength objects are equal
                    II = cellfun(@(x) eq(set(self.object, 'wl', []), set(x, 'wl', []), false), {data.object});
                    I = find(I & II);
                    if I
                        % act according to force on the first match
                        % (attention note below applies)
                        data(I(1)) = append(self, data(I(1)), force);
                    else
                        switch force
                            case -2
                                data = [self, data];
                            case 2
                                data = [data, self];
                            otherwise
                                error('Cannot expand dataset, geometries appear to be different.');
                        end
                    end
                else
                    % attention: output has to be the same size as data(I)
                    % In principle, should be the case, as, even if passed
                    % with abs(force)==2, objects in data and self have
                    % already been compared.
                    data(I) = append(self, data(I), force);
                end
                self = data;
                return
            end
            
            %%% At this stage data, in addition to self, is considered a scalar object
            
            % this is rather for empty instances that cause errors otherwise
            if isempty(self)
                self = data;
                return
            elseif all(isempty(data))
                return
            end
            
            if imag(force)
                force = imag(force);
            else
                % wavelengths will be taken care of later
                assert(set(self.object,'wavelength',[]) == set(data.object,'wavelength',[]), ...
                    'Dataset objects appear to be different. Use imaginary force parameter to igone this.')
            end
            % if the wavelength arrays are not equal, check if we can
            % interpolate or downsample the values, or if we should simply
            % append to the Dataset array
            if ~isequal(self.wavelength, data.wavelength)
                if force == -2
                    self = [self, data];
                    return
                elseif force == 2
                    self = [data, self];
                    return
                end
                
                if length(self.wavelength) > 1 && length(data.wavelength) > 1 && isequal(self.wavelength([1 end]), data.wavelength([1 end]))
                    if ~force
                        while true
                            reply = input('Wavelength properties have different step sizes, interpolate [i], downsample [d], or quit [q]? ([i] by default): ', 's');
                            if isempty(reply) || strcmpi(reply, 'i')
                                force = 1;
                                break
                            elseif strcmpi(reply, 'd')
                                force = -1;
                                break
                            elseif strcmpi(reply, 'q')
                                error('Cannot expand dataset, wavelengths appear to be different')
                            else
                                disp("Didn't catch it, try again.");
                            end
                        end
                    end
                    
                    if force == 1 % force interpolation
                        % disp('Interpolating....')
                        if length(self.wavelength) < length(data.wavelength)
                            self = interpData(self, data.wavelength);
                        else
                            data = interpData(data, self.wavelength);
                        end
                    elseif force == -1 % force downsampling
                        % disp('Downsampling....')
                        if length(self.wavelength) < length(data.wavelength)
                            data = interpData(data, self.wavelength);
                        else
                            self = interpData(self, data.wavelength);
                        end
                    end 
                        
                elseif ~isempty(intersect(self.wavelength, data.wavelength))
                    if force == -1
                        self = interpData(self, intersect(self.wavelength, data.wavelength));
                        data = interpData(data, intersect(data.wavelength, data.wavelength));
                    else
                        error(['Cannot expand dataset, wavelengths appear to be different. Cannot extrapolate.\n',...
                            'Use option force == -1 or -1i to downsample at the intersection wavelengths.'])
                    end
                else
                    error('Cannot expand dataset, wavelengths appear to be different')
                end
            end
            
            % not necessarily "better" than having a case for each
            % property, but at least we do not have to hard-code their
            % names
            props = self.property;
            for ii = 1:length(props)
                if ~isempty(self.(props{ii})) % to prevent annoying 0x0xN matrices and infinite loops on 'start'
                    switch class(self.(props{ii}))
                        case 'double'
                            self.(props{ii}) = cat(3, self.(props{ii}), data.(props{ii}));
                        case 'cell'
                            self.(props{ii}) = [self.(props{ii}), data.(props{ii})];
                        case 'StartDataset'
                            self.(props{ii}) = append(self.(props{ii}), data.(props{ii}), force*1i);
                            % explanation for the force argument: if we reached this line,
                            % objects are either compared or ignored already, if force was
                            % 0, we don't ask the second time for the start Dataset
                        case {'string', 'char'}
                            % add to string array strings that are not equal
                            self.(props{ii}) = [self.(props{ii}), data(~strcmpi(self.(props{ii}), data.(props{ii}))).(props{ii})];
                    end
                end
            end
        end % function append
        
        function self = plus(self, data, varargin)
        % basically redundant, but it's tied to '+' and I'm used to append
            self = self.append(data, varargin{:});
        end
        
        function self = minus(self, data, ignoreObject)
        % when substracting with ignoreObject option set to true,
        % we compare the patterns (through hashes) only, nothing else, 
        % thus the data could be for different bands, but as
        % long as the patterns are identical, they're considered equal.
        % ignoreObject = true by default, since usually, it would mean the patterns are
        % the result of the same calculation, we just have modified the
        % band, which we can always do again
            if nargin < 3
                ignoreObject = true;
            end
            if ~isa(data, 'Dataset')
                data = Dataset(data);
            end
            if self.object ~= data.object && ~ignoreObject
                error("The two datasets seem to have different ojects, use ignoreObject=true to proceed anyway.")
            end
            self = self.remove(data.hash);
        end
        
        function self = merge(self, data, force, opt)
        % self and data can each be Dataset arrays which are merged into a
        % single Dataset object. See append for force options.
        % opt can be specified as 'nounique' to skip the use of the unique
        % function
            if nargin < 3
                force = 0;
            end
            if nargin < 2
                data = [];
            end
            if length(self) > 1
                [self, data] = deal(self(1), [self(2:end), data]);
            end
            for ii = 1:length(data)
                self = self.append(data(ii), force);
            end
            if nargin == 4 && strcmpi(opt,'nounique')
            elseif nargin == 4 && ~strcmpi(opt,'nounique')
                error("Argument #4 can only be 'nounique'")
            else
                self = unique(self);
            end
        end
        
        function self = sum(self, varargin)
        % merges array of Dataset instances into a single one. 
        % See merge for varargin options.
            self = merge(self(1), self(2:end), varargin{:});
        end
        
        function [tf, equals] = eq(obj1, obj2)
        % checks equality of two dataset objects, or arrays thereof,
        % tf is the answer whether they are equal or not, 
        % equals is a logical array of elements that are. 
        % Note: it's not very good for much, since it compare two arrays in
        % sequence, checking if they are basically identical copies of each
        % other. For sorting-independent inquiries, use functions like
        % setdiff and such
            if length([obj1.hash]) == length([obj2.hash])
                equals = strcmp(obj1.hash, obj2.hash);
                tf = all(equals);
            else
                error(['Two sets seem to have a different number of elements. ',...
                    'Use other functions, like setdiff and such, for more info.'])
            end
        end
        
        function [d, I] = reduce(d, tol_pz, tol_dp, tol_d2p)
        % reduces the number of entries in the dataset by considering all
        % the devices in the phase plane (volume) of size tol_pz x tol_dp
        % (x tol_d2p) and choosing only one with the highest absEff among
        % them. Units = radians
        % if tol_d2p is not given (or is empty), first order fit is used. 
        % If only tol_pz is given, tol_dp = tol_pz.
        % I - logical vector of length d.numdev where 1 stand for entries
        % to keep.
            if nargin < 4 || isempty(tol_d2p)
                order = 1;
            else
                order = 2;
            end
            if nargin < 2; tol_pz = 0.025*pi; end
            if nargin < 3; tol_dp = tol_pz; end
            if numel(d) > 1
                if nargin == 4
                    [d, I] = arrayfun(@(d)reduce(d, tol_pz, tol_dp, tol_d2p), d, 'UniformOutput', false);
                else
                    [d, I] = arrayfun(@(d)reduce(d, tol_pz, tol_dp), d, 'UniformOutput', false);
                end
                d = reshape([d{:}], size(d));
                return
            end
            
            for ipol = 1:size(d.phiCoeffs,2)
                if ~isempty(d.phiCoeffs)
                    coeffs = squeeze(d.phiCoeffs(:,ipol,:,order)); 
                else
                    error('phiCoeffs property empty. Run getPhi function on the Dataset to calculate them.')
                end
                avreff = squeeze(mean(d.absEff, 1));
                dpGrid = min(coeffs(2,:)) : tol_dp : max(coeffs(2,:)) + tol_dp;
                pzGrid = min(coeffs(3,:))  : tol_pz : max(coeffs(3,:)) + tol_pz;
                if order ==  1
                    d2pGrid = [0 1.e-6];
                else
                    d2pGrid = min(coeffs(1,:))  : tol_d2p : max(coeffs(1,:)) + tol_d2p;
                end
                I = false(1,d.numdev);
                for ii = 1:length(dpGrid)-1
                    for jj = 1:length(pzGrid)-1
                        for kk = 1:length(d2pGrid)-1
                            J = coeffs(2,:) >= dpGrid(ii) & coeffs(2,:) < dpGrid(ii+1) & ...
                                coeffs(3,:) >= pzGrid(jj) & coeffs(3,:) < pzGrid(jj+1) & ...
                                coeffs(1,:) >= d2pGrid(kk) & coeffs(1,:) < d2pGrid(kk+1);
                            [~, ia] = max(avreff(J));
                            f = find(J);
                            I(f(ia)) = true;
                        end
                    end
                end
                if ipol == 1
                    Itmp = I;
                else
                    I = I & Itmp;
                end
            end
            d = d.get(I);
        end
        
        %------------------------------------------
        %%% End operations methods
        %------------------------------------------
        
        %------------------------------------------
        %%% Lookup and get methods
        %------------------------------------------
        
        function tf = isempty(self)
            if numel(self) > 1
                tf = arrayfun(@isempty, self);
                return
            end
            if numel(self) == 0
                tf = true;
                return
            end
            tf = isempty(self.phase) && isempty(self.absEff) && ...
                isempty(self.relEff) && isempty(self.pattern);
        end
        
        function [tf, I] = contains(self, device)
        % device can be either a pattern, hash, or an Optimizer object, or
        % a matrix of those, or even cell array of objects of mixed kind
        % tf - True/False answer
        % I - index where match was found
        
            if numel(device) > 1
                [tf, I] = arrayfun(@(d)contains(self, d), device);
                return
            end
            if iscell(device)
                device = device{1};
            end
            if isa(device, 'Optimizer')
                device = device.finalPattern;
            elseif isnumeric(device)
                device = DataHash(device, 'SHA-1');
            elseif ~ischar(device)
                error('Wrong device format');
            end
            I = find(strcmp(self.hash, device));
            if I
                tf = true;
            else
                tf = false;
            end
        end
        
        function self = remove(self, device)
        % device can be specified as either a:
        %   - scalar index
        %   - array of indices
        %   - hash/cell arry of hashes
        %   - pattern/cell array of patterns
        %   - Optimizer object or their array/cell array
            if numel(self) > 1
                self = arrayfun(@(x) remove(x, device), self);
                return
            end
            if isa(device, 'Optimizer')
                device = device.finalPattern;
            end
            if isnumeric(device) || islogical(device)
                if all(size(device) > 1) % single pattern image
                    device = {device}; 
                else % a single or an arry of indices
                    I = device;
                end
            elseif ischar(device) || ... % single hash
                   isa(device, 'function_handle') % single generator function
                device = {device}; 
            end
            
            if iscell(device)
                if ~ischar(device{1}) % not already a hash
                    device = cellfun(@(x) DataHash(x, 'SHA-1'), device, 'UniformOutput', false);
                end
                % 'UniformOutput',false prevents error if one of the
                % devices is not found in the hash array
                I = cellfun(@(x) find(strcmp(self.hash, x)), device, 'UniformOutput', false);
                I = cell2mat(I);
            end
            % at this stage, I stores indices to be removed
            props = self.property;
            for ii = 1:length(props)
                if ~isempty(self.(props{ii})) % to prevent annoying 0x0xN matrices and infinite loops on 'start'
                    switch class(self.(props{ii}))
                        case 'double'
                            self.(props{ii})(:,:,I,:) = [];
                        case 'cell'
                            self.(props{ii})(I) = [];
                        case 'StartDataset'
                            self.(props{ii}) = remove(self.(props{ii}), I);
                    end
                end
            end
        end
        
        function [I, devs] = find(self, device)
        % device is specified by a scalar or vector of either a pattern, 
        % a hash, or an Optimizer
            if isa(device, 'Optimizer')
                device = device.finalPattern;
            end
            if ischar(device) || isnumeric(device) || ... % single hash or single pattern
                    isa(device, 'function_handle') % single generator function
                device = {device};
            end
            if ~ischar(device{1}) 
                device = cellfun(@(x) DataHash(x, 'SHA-1'), device, 'UniformOutput', false);
            end
            % 'UniformOutput',false, again, to deal with empty vectors,
            % when match not found
            I = arrayfun(@(x) find(strcmp(self.hash, x)), device, 'UniformOutput', false);
            I = cell2mat(I);
            if nargout == 2
                devs = self.get(I);
            end
        end
        
        function self = get(self, I)
        % returns subset of devices given by indices in I
            if numel(self) > 1
                if iscell(I) && numel(I) == numel(self)
                    self = arrayfun(@(x,y)get(x,y{1}),self,I);
                elseif ~iscell(I)
                    self = arrayfun(@(x)get(x,I),self);
                else
                    error('Wrong input')
                end
                return
            end
            I = I(:);
            props = self.property;
            for ii = 1:length(props)
                if ~isempty(self.(props{ii})) % to prevent annoying 0x0xN matrices and infinite loops on 'start'
                    switch class(self.(props{ii}))
                        case 'double'
                            self.(props{ii}) = self.(props{ii})(:,:,I,:);
                        case 'cell'
                            self.(props{ii}) = self.(props{ii})(I);
                        case 'StartDataset'
                            self.(props{ii}) = get(self.(props{ii}),I);
                    end
                end
            end
        end
        
        function [self, I] = rmnan(self)
        % remove patterns that are all nans
            I = find(cellfun(@(x) any(isnan(x(:))), self.pattern));
            self = remove(self, I);
        end
        
        function N = numdev(self)
        % returns the number of devices in the Dataset or their array
            if numel(self) > 1
                N = arrayfun(@numdev, self);
                N = sum(N(:));
                return
            elseif numel(self) == 0
                N = 0;
                return
            end
            if ~isempty(self.pattern)
                N = length(self.pattern);
            elseif ~isempty(self.phase)
                N = size(self.phase, 3);
            elseif ~isempty(self.absEff)
                N = size(self.phase, 3);
            else
                N = 0;
            end
        end
        
        function N = numsdev(self)
            N = numdev([self.start]);
        end
        
        function N = numsunique(self)
            s = [self.start];
            s = s.recalcHash;
            N = length(unique([s.hash]));
        end
        
        function phase = targetPhase(self)
            if numel(self) > 1
                phase = arrayfun(@targetPhase, self, 'UniformOutput', false);
                return
            end
            if size(self.tgtPhiCoeffs,1) == 3 % coefficients given
                phase = self.get_phase(self.tgtPhiCoeffs);
            else
                phase = self.tgtPhiCoeffs;
            end
        end
        
        function wl = get.wavelength(self)
            wl = self.object.wavelength;
        end
        
        function self = set.wavelength(self, wl)
            self.object.wavelength = wl;
        end
        
        function self = set.pattern(self, p)
            if ~iscell(p) && ~isempty(p), p = {p}; end
            self.pattern = p;
        end
        
        function props = get.property(self)
            m = metaclass(self);
            props = {m.PropertyList.Name};
            props = props(~[m.PropertyList.Dependent]);
        end
        
        %------------------------------------------
        %%% End lookup and get methods
        %------------------------------------------
        
        %------------------------------------------
        %%% Set methods
        %------------------------------------------
        
        function [C, ia, ib] = intersect(A, B)
        % Similar to built in intersect, overloaded for Dataset. See help intersect.
            [~, ia, ib] = intersect(A.hash, B.hash, 'stable');
            C = A.get(ia);
        end
        
        function [Lia, Locb] = ismember(A, B)
        % Similar to built in ismember, overloaded for Dataset. See help ismember.
            [Lia, Locb] = ismember(A.hash, B.hash, 'stable');
        end
        
        function [C, ia] = setdiff(A, B)
        % Similar to built in setdiff, overloaded for Dataset. See help setdiff.
            [~, ia] = setdiff(A.hash, B.hash, 'stable');
            C = A.get(ia);
        end
        
        function [C, ia, ib] = setxor(A, B)
        %  returns the data of A and B that are not in their intersection, 
        % with no repetitions, in other words, the data that occurs in A or B, but not both
            [~, ia, ib] = setxor(A.hash, B.hash, 'stable');
            C = A.get(ia) + B.get(ib);
        end
        
        function [C, ia, ib] = union(A, B)
        % Similar to built in union, overloaded for Dataset. See help union.
        % Kind of repeats basic functionality of merge, but merge is still
        % useful on dataset arrays.
            [~, ia, ib] = union(A.hash, B.hash, 'stable');
            C = A.get(ia) + B.get(ib);
        end
        
        function [C, ia, ic] = unique(A)
        % Similar to built in unique, overloaded for Dataset. See help unique.
            if numel(A) > 1
                sz = size(A);
                [C, ia, ic] = arrayfun(@unique, A, 'UniformOutput', false);
                C = [C{:}];
                C = reshape(C, sz);
                return
            end
            if numel(A.hash) ~= numel(A.pattern)
                A = recalcHash(A);
            end
            [~, ia, ic] = unique(A.hash, 'stable');
            C = A.get(ia);
        end
        
        %------------------------------------------
        %%% End set methods
        %------------------------------------------
        
        %----------------------------------------------------
        %%% Filtering and sorting methods 
        %----------------------------------------------------
        
        function self = cutwl(self, varargin)
        % cutwl(d, arr) -- arr=array of discreet wavelengths of length >= 1
        % cutwl(d, wl_min, wl_max) -- closest min and max
        % leaves data in the dataset that is related to a specific single
        % wavelength, bounded by some min and max values, or for specific 
        % wavelengths in an array arr
            if numel(self) > 1
                self = arrayfun(@(s) cutwl(s,varargin{:}), self);
                return
            end
            [wl_max, wl_min, I] = deal([]);
            switch length(varargin)
                case 1
                    [~, I] = intersect(self.wavelength, varargin{1});
                case 2
                    [wl_min, wl_max] = varargin{:};
                otherwise
                    error('Too many input arguments');
            end
            if isempty(I)
                [~, imin] = min(abs(self.wavelength - wl_min));
                [~, imax] = min(abs(self.wavelength - wl_max));
                I = imin:imax;
            end
            self.wavelength = self.wavelength(I);
            
            props = self.property;
            for ii = 1:length(props)
                if ~isempty(self.(props{ii}))
                    switch class(self.(props{ii}))
                        case 'double'
                            if ~contains(props{ii}, 'hiCoeffs') 
                                self.(props{ii}) = self.(props{ii})(I,:,:);
                            end
                        case 'StartDataset'
                            self.(props{ii}) = cutwl(self.(props{ii}), varargin{:});
                    end
                end
            end
        end
        
        function [self, I] = cuteff(self, eff, opt)
            if nargin < 3
                opt = 'mean';
            end
            switch opt
                case 'mean'
                    fun = @(x) reshape(mean(x, 1), [], 1);
                case 'min'
                    fun = @(x) reshape(min(x, [], 1), [], 1);
            end
            I = fun(self.absEff) > eff;
            self = self.get(I);
        end
        
        function [self, I] = cutphi(self, opt, varargin)
        % [self, I] = cutphi(self, opt, phi_min, phi_max)
        % [self, I] = cutphi(self, opt, [phi_min, phi_max])
        % [self, I] = cutphi(self, opt, phi_val)
        % [self, I] = cutphi(_____, polopt)
        % leaves only entries whose phase fit coefficient, specified by opt,
        % is in the given interval, or equal to a specific value phi_val.
        % opt can be 'pz' or 'dp' for phi_zero and delta_phi in linear fit,
        % or 'd2p' for delta^2_phi in the quadratic one
        % I is a logical vector such that self_new = self_old.get(I);
        % polopt is a char array which is either 'both' or 'either'. It is
        % appropriate in the case of both polarizations, when we have two
        % different sets of phase coefficient and we want to leave the
        % intersection or union of their subsets, respectively. Default:
        % 'either'.
        
            if numel(self) > 1
                [self, I] = arrayfun(@(x)cutphi(x,opt,varargin{:}), self, 'UniformOutput', false);
                self = reshape([self{:}], size(self));
                return
            end
            polopt = 'either';
            switch length(varargin)
                case 3
                    [phi_min, phi_max, polopt] = varargin{:};
                case 2
                    [phi_min, phi_max] = varargin{:};
                    if length(phi_min) == 2
                        polopt = phi_max;
                        phi_max = phi_min(2);
                        phi_min = phi_min(1);
                    end
                case 1
                    phi_min = varargin{1}(1);
                    phi_max = varargin{1}(end);
            end
            for ipol = 1:size(self.phiCoeffs,2)
                if ~isempty(self.phiCoeffs)
                    switch opt
                        case 'dp'
                            coeffs = squeeze(self.phiCoeffs(2,ipol,:,1));
                        case 'pz'
                            coeffs = squeeze(self.phiCoeffs(3,ipol,:,1));
                        case 'd2p'
                            coeffs = squeeze(self.phiCoeffs(1,ipol,:,2));
                    end
                else
                    error('phiCoeffs property empty. Run getPhi function on the Dataset to calculate them.')
                end
                I = coeffs >= phi_min & coeffs <= phi_max;
                if ipol == 1
                    Itmp = I;
                else
                    switch polopt
                        case 'both'
                            I = Itmp & I;
                        case 'either'
                            I = Itmp | I;
                        otherwise
                            error("polopt argument can be either 'both' or 'either'")
                    end
                end
            end
            self = self.get(I);
        end
        
        function [self, I] = cutpz(self, varargin)
        % alias for cutphi with the pz option
            [self, I] = cutphi(self, 'pz', varargin{:});
        end
        
        function [self, I] = cutdp(self, varargin)
        % alias for cutphi with the dp option
            [self, I] = cutphi(self, 'dp', varargin{:});
        end
        
        function [self, I] = cutoutliers(self)
        % removes outliers in terms of the delta phi value
            if numel(self) > 1
                [self, I] = arrayfun(@cutoutliers, self, 'UniformOutput', false);
                self = reshape([self{:}], size(self));
            else
                med = median(squeeze(self.phiCoeffs(2,:,:,1)));
                sig = sqrt(std(squeeze(self.phiCoeffs(2,:,:,1))));
                [self, I] = cutphi(self, 'dp', med-sig, med+sig);
            end
        end
        
        function [self, I] = cutff(self, ffmin, ffmax)
        % filters the dataset, leaving the entries whose filling fraction
        % of the material is between the ffmax and ffmin values.
        % I is a logical vector such that self_new = self_old.get(I);
            if nargin < 3 || isempty(ffmax)
                ffmax = 1;
            end
            if isempty(ffmin)
                ffmin = 0;
            end
            I = false(1,self.numdev);
            for ii = 1:numel(I)
                pat = self.gp(ii);
                ff = sum(pat(:))/numel(pat);
                I(ii) = ff >= ffmin & ff <= ffmax;
            end
            self = self.get(I);
        end
        
        function [self, I] = filtMinSize(self, varargin)
        % [self, I] = filtMinSize(self, mfs) - checks both min distance 
        % (equal to mfs) and area (equal to pi*(mfs/2)^2)
        % [self, I] = filtMinSize(self, mfs, mfa) - checks min distance mfs
        % and min area mfa
        % [self, I] = filtMinSize(self, mfs, 'd') - checks only min
        % distance mfs
        % [self, I] = filtMinSize(self, mfa, 'a') - checks only min area
        % mfa
        % ______ = filtMinSize(_______, 'inverted', true/false) - flag that
        % says whether to look at the image where pixel values are flipped
        % Units are those of self.object.unit.
        % Outputs: self - filtered Dataset, I - array of length self.numdev,
        % where true's are in the position of entries that passed one or
        % both checks, according to the input, and false's for those that
        % did not (false by default).
        
            % varargin pre-processing. At the moment there's only one real
            % name/value pair. If more are added, one should switch to
            % using updatestruct, as usual.
            if length(varargin) > 2 % inverted option given
                assert(strcmp(varargin{end-1}, 'inverted'), 'Wrong input arguments.')
                try_inverted = varargin{end};
                varargin = varargin(1:end-2);
            else
                try_inverted = false; 
            end
            if length(varargin) == 1
                mfs = varargin{1};
                mfa = pi*(mfs./2).^2;
            elseif length(varargin) == 2
                switch varargin{2}
                    case 'd'
                        mfs = varargin{1};
                        mfa = 0;
                    case 'a'
                        mfs = 0;
                        mfa = varargin{1};
                    otherwise
                        [mfs, mfa] = varargin{:};
                end
            else
                error('Wrong input arguments.')
            end
            I = true(size(self.pattern));
            gridSize = @(x) mean(size(x)./self.object.period);
            % try again only those not already rejected by the previous test
            % Narrow feature test is stricter, so it goes first
            ispost = polarity(self);
            patterns = self.gp;
            patterns(~ispost) = cellfun(@(x)~x, patterns(~ispost), 'UniformOutput', false);
            if mfa
                I = cellfun(@(x) checkNarrowFeat(x, 2*sqrt(mfa/pi)*gridSize(x), try_inverted, true), patterns);
                I(I) = cellfun(@(x) checkMinArea(x, mfa*gridSize(x)^2, try_inverted), patterns(I));
            end
            if mfs 
                I(I) = cellfun(@(x) checkNarrowFeat(~x, mfs*gridSize(x), try_inverted, true), patterns(I));
                I(I) = cellfun(@(x) checkMinDist(x, mfs*gridSize(x), try_inverted), patterns(I));
            end
            self = self.get(I);
        end
        
        function [d, dd, ispost] = polsort(d)
        % splits Dataset d into two: d - with patterns being hole-like 
        % (continuous region of 1s), and dd - post-like (continuous region of 0s).
        % ispost is a logical vector of length being the input number of 
        % patterns, true stands for posts and false for holes
            ispost = polarity(d.gp);
            if nargout > 1
                dd = d.get(ispost);
            end
            d = d.get(~ispost);
        end
        
        function [d, ia, ib] = filtpershift(d) 
        % finds equivalent patterns in d produced by periodic shifts along
        % x and/or y directions. Several options are possible:
        % - no symmetry: any periodic shift is possible; can be slow, since
        % no hash is used, instead, two patterns are considered equivalent
        % if sum of all elements of their difference is equal to zero. This
        % algorithm produces false-positives for some symmetric patterns
        % with equal filling fractions;
        % - symmetry in either x or y: shift patterns by half of the period 
        % in the corresponding direction and compare hashes;
        % - both symmetries + diagonal: shift in both directions at the
        % same time only (attention to have diagonal symmetry flag property 
        % set for square periods);
        % - both symmetries, no diagonal: try all three combinations.
        % The reason we don't have to consider arbitraty shifts for the
        % symmetric case (in addition to false-possitive issue) as well as
        % x or y only shifts for the diagonal case is that, even though such
        % patterns would really be equivalent, they would not be produced
        % by something like Optimizer.enforceSymmetry and we couldn't have
        % used symmetry planes in simulations.
        % Outputs:
        % - d - original Dataset with duplicates removed
        % - ia - indices of entries removed
        % - ib - indices of equivalent patterns left in d that correspond
        % to those in ia.
        % One can visualize all the equivalent patterns found by running:
        % [~, ia, ib] = filtpershift(d);
        % ub = unique(ib);
        % arrayfun(@(n)plotrandgrid(d.get([n ia(ib == ub)]),'ordered'),ub),
        % where image 1 will be the one that will be kept, all other -
        % removed.
        
            shift = d.object.sym(1:2)/2;
            shift = shift(:).'; % just in case
            if all(shift == 0)
                % the case of no symmetry. For symmetric case it's a bit of
                % overkill, plus, produces false-positives for some
                % patterns with equal fill fraction
                check = ones(d.numdev);
                for ii = 1:numdev(d)
                    for jj = ii+1:numdev(d)
                        if all(size(d.pattern{ii}) == size(d.pattern{jj}))
                            check(jj,ii) = sum(sum(d.gp(ii) - d.gp(jj)));
                        end
                    end
                end
                check = ~check;
            else
                if ~any(shift == 0) && d.object.sym(3) % both symmetries w/ diagonal
                    shift = [diag(shift); shift];
                end
                check = false(d.numdev);
                for ii = 1:size(shift, 1)
                    d_shifted = patternshift(d, shift(ii,:));
                    for jj = 1:numdev(d)
                        check(jj+1:end,jj) = check(jj+1:end,jj) | ...
                            strcmp(d.hash{jj}, d_shifted.hash(jj+1:end)).';
                    end
                end
            end
            f = find(check);
            [ia, ib] = ind2sub(size(check), f);
            d = d.remove(unique(ia));
        end
        
        function [d, I] = filtnoncont(d)
        % removed non-continuous patterns from Dataset D, I is a logical
        % array labeling the patterns that were kept/removed.
            brdunq = @(p) unique([p([1 end],:), p(:,[1 end]).']);
            I = false(1,d.numdev);
            for ii = 1:d.numdev
                p = d.gp(ii);
                I(ii) = numel(brdunq(p)) == 1;
                if ~I(ii)
                    I(ii) = numel(brdunq(patternshift(p, 0.5))) == 1;
                end
            end
            d = d.get(I);
        end
        
        function [d, I] = filtnoncontsuper(d)
            brdunq = @(p) unique([p([1 end],:), p(:,[1 end]).']);
            I = false(1,d.numdev);
            for ii = 1:d.numdev
                p = d.gp(ii);
                I(ii) = numel(brdunq(supercell(p, 'rot', 45, 1))) == 1;
                if ~I(ii)
                    I(ii) = numel(brdunq(supercell(patternshift(p, 0.5), 'rot', 45, 1))) == 1;
                end
            end
            d = d.get(I);
        end
        
        function d = splitpol(d)
        % Splits a single Dataset object that contains data for two
        % polarizations into an array of two Dataset objects with the data
        % for the TE and TM polarizations, respectively.
            if numel(d) > 1
                d = arrayfun(@splitpol, d, 'UniformOutput', false);
            else
                if d.object.numPol == 2
                    d = repelem(d,2);
                    props = {'phase','absEff','relEff','tgtPhiCoeffs','phiCoeffs'};
                    for i = 1:length(d)
                        for j = 1:length(props)
                            if ~isempty(d(i).(props{j}))
                                tmp = d(i).(props{j});
                                d(i).(props{j}) = tmp(:,i,:,:);
                            end
                            d(i).object.pol = i;
                        end
                    end
                    if ~isa(d,'StartDataset')
                        ds = splitpol(d(1).start);
                        for i = 1:length(d)
                            d(i).start = ds(i);
                        end
                    end
                end
            end
        end
        
        %----------------------------------------------------
        %%% End filtering and sorting methods 
        %----------------------------------------------------
        
        %----------------------------------------------------
        %%% Generator methods (calculate some new quantity)
        %----------------------------------------------------
        
        function [self, phiCoeffs] = getPhi(self, wl0)
        % calculates phiCoeffs property
        % phiCoeffs: tensor of size 3 x numPol x numDevices x order, where 
        % in the first dimension we have \Delta^2\phi, \Delta\phi, \phi_0, 
        % in accordance to matlab's polyfit convention, and in the optional
        % 4th dimension we have a results for different fits with different
        % orders
        % phase difference is taken with respect to the phase at wl0 whicih
        % is taken to be a central wavelength by default
            if nargin < 2
                wl0 = []; % don't assign it now in case self is an array
            end
            if numel(self) > 1
                if nargout == 1
                    self = arrayfun(@(x) getPhi(x, wl0), self);
                else
                    [self, phiCoeffs] = arrayfun(@(x) getPhi(x, wl0), self, 'UniformOutput', false);
                end
                return
            end
            if isempty(self) || isempty(self.phase)
                phiCoeffs = [];
                return
            end
            
            % so far, for simplicity, assume numbers of wavelengths and
            % polarizations are the same
            phiCoeffs = zeros(3, size(self.phase, 2), size(self.phase, 3), 2);
            if numel(self.wavelength) < 3
                phiCoeffs(3,:,:,:) = repmat(self.phase(end,:,:),1,1,1,2);
                phiCoeffs(2,:,:,:) = repmat(self.phase(1,:,:) - self.phase(end,:,:),1,1,1,2);
                self.phiCoeffs = phiCoeffs;
                self.start = getPhi(self.start);
                return
            end
            if isempty(wl0)
                wl0 = sum(self.wavelength([end 1]))/2;
            end
            self = centerPhase(self, wl0);
            dwl = diff(self.wavelength([end 1]));
            for order = 1:2
                for polIter = 1:size(self.phase, 2) % over polarization
                    for devIter = 1:size(self.phase, 3) % over devices
                        p = polyfit(self.wavelength' - wl0, ...
                            self.phase(:, polIter, devIter), order);
                        % converting from polynomial coeffs to derivatives
                        if order == 1
                            phiCoeffs(2, polIter, devIter, order) = p(1)*dwl;
                            phiCoeffs(3, polIter, devIter, order) = p(2);
                        else
                            phiCoeffs(1, polIter, devIter, order) = p(1)*2*(dwl/2)^2;
                            phiCoeffs(2, polIter, devIter, order) = p(2)*(dwl/2);
                            phiCoeffs(3, polIter, devIter, order) = p(3);
                        end
                    end
                end
            end
            self.phiCoeffs = phiCoeffs;
            self.start = getPhi(self.start, wl0); % may comment it and call directly separately if needed
        end % getPhi
        
        function phase = get_phase(self, coeffs)
            if numel(self) > 1
                if nargin == 1
                    phase = arrayfun(@get_phase, self, 'UniformOutput', false);
                else
                    phase = arrayfun(@(x) get_phase(x,coeffs), self, 'UniformOutput', false);
                end
                return
            end
            if nargin == 1
                coeffs = self.phiCoeffs(:,:,:,1);
            elseif size(coeffs, 4) > 1 % by default choose linear fit
                coeffs = coeffs(:,:,:,1);
            end
            if isempty(coeffs)
                phase = [];
                return
            end
            dwl = diff(self.wavelength([end 1]));
            if all(coeffs(1, :) == 0)
                coeffs(2, :, :) = coeffs(2, :, :)/dwl;
            else
                coeffs(1, :, :) = coeffs(1, :, :)/2/(dwl/2)^2;
                coeffs(2, :, :) = coeffs(2, :, :)/(dwl/2);
            end
            phase = nan(length(self.wavelength), size(coeffs, 2), size(coeffs, 3));
            if length(self.wavelength) == 2
                wl_center = self.wavelength(end);
            else
                wl_center = sum(self.wavelength([end 1]))/2; % has to match wl0 of getPhi
            end
            for polIter = 1:size(phase, 2) % over polarization
                for devIter = 1:size(phase, 3) % over devices
                    phase(:, polIter, devIter) = polyval(coeffs(:, polIter, devIter),...
                        self.wavelength' - wl_center);
                end
            end
        end
        
        function self = centerPhase(self, wl)
            if nargin == 1
                wl = sum(self.wavelength([end 1]))/2;
            end
            for ii = 1:numel(self.phase(1,:))
                ph_wl = interp1(self.wavelength, self.phase(:,ii), wl, 'pchip');
                self.phase(:,ii) = self.phase(:,ii) - floor(floor(ph_wl/pi+1e-4)/2)*2*pi;
            end
        end
        
        function self = recalcHash(self)
            if numel(self) > 1
                self = arrayfun(@recalcHash, self);
                return
            end
            self.hash = cellfun(@(x) DataHash(x, 'SHA-1'), self.pattern, 'UniformOutput', false);
        end
        
        function [self, phase] = prettifyPhase(self)
            if numel(self) > 1
                sz = size(self);
                [self, phase] = arrayfun(@(x) prettifyPhase(x), self, 'UniformOutput', false);
                self = [self{:}]; % transform back to array of Datasets
                self = reshape(self, sz);
                return
            end
            if ~isempty(self)
                phase = unwrap(self.phase, [], 1);
                N = ceil(size(phase, 1)/2);
                phi0 = mod(self.phase(N, :, :), 2*pi);
                phi0(phi0 < 0) = phi0(phi0 < 0) + 2*pi;
                phase = phase - phase(N, :, :) + phi0;
                if ~isempty(self.tgtPhiCoeffs)
                    tgtphase = self.targetPhase;
                    phase = phase + round((tgtphase(N,:,:) - phase(N,:,:))/2/pi)*2*pi;
                end
                self.phase = phase;
            else
                phase = [];
            end
            if ~isempty(self.start); self.start = prettifyPhase(self.start); end
        end
        
        function out = interpData(self, wl)
            if numel(self) > 1
                out = arrayfun(@(x) interpData(x, wl), self);
                return
            end
            out = self;
            if length(out.wavelength) ~= length(wl) || ~all(out.wavelength == wl)
                out.wavelength = wl;
            else
                % wavelength arrays are idential, nothing to do
                return
            end
            % interpolate all devices independently, interp3 won't work (I guess)
            [~, numPol, numDev] = size(self.phase); % assume we always have non-empty phase
            out.phase = zeros(length(wl), numPol, numDev);
            for polIter = 1:numPol
                for devIter = 1:numDev
                    out.phase(:, polIter, devIter) = interp1(self.wavelength, ...
                        self.phase(:, polIter, devIter), wl, 'spline');
                end
            end
            if ~isempty(self.absEff)
                out.absEff = zeros(length(wl), numPol, numDev);
                for polIter = 1:numPol
                    for devIter = 1:numDev
                        out.absEff(:, polIter, devIter) = interp1(self.wavelength, ...
                            self.absEff(:, polIter, devIter), wl, 'spline');
                    end
                end
            end
            if ~isempty(self.relEff)
                out.relEff = zeros(length(wl), numPol, numDev);
                for polIter = 1:numPol
                    for devIter = 1:numDev
                        out.relEff(:, polIter, devIter) = interp1(self.wavelength, ...
                            self.relEff(:, polIter, devIter), wl, 'spline');
                    end
                end
            end
            if ~isempty(self.start); out.start = interpData(self.start, wl); end
        end
        
        function self = dataset2struct(self)
            %dataset2struct converts the Dataset object into Matlab's struct. 
            % Not a 'clean' way of converting: using struct function we
            % copy all the fields: public, private, hidden, dependent. For
            % a 'clean' copy of public non-hidden properties, use mystruct 
            % method. Rename mystruct to struct to overload the latter and
            % make this method 'clean' as well.
            warning('off', 'MATLAB:structOnObject')
            self = arrayfun(@struct, self);
            if all(isempty([self.start]))
                [self.start] = deal([]);
            else
                tmp = arrayfun(@dataset2struct, [self.start], 'UniformOutput', false);
                [self.start] = tmp{:};
            end
            tmp = arrayfun(@optimizer2struct, [self.object], 'UniformOutput', false);
            [self.object] = tmp{:};
        end
        
        function structself = mystruct(self)
            % Converts Dataset object to struct, copying only public
            % non-hidden properties.
            if numel(self) > 1
                structself = arrayfun(@mystruct, self);
                return
            end
            structself = struct;
            names = properties(self);
            for n = 1:length(names)
                structself.(names{n}) = self.(names{n});
            end
            tmp = arrayfun(@dataset2struct, [self.start], 'UniformOutput', false);
            [structself.start] = tmp{:};
            tmp = arrayfun(@optimizer2struct, [self.object], 'UniformOutput', false);
            [structself.object] = tmp{:};
        end
        
        function self = simDataset(self, varargin) 
        % self = simDataset(self)
        % self = simDataset(self, wl)
        % self = simDataset(self, robustopt)
        % self = simDataset(self, wl, robustopt)
        % run forwardRun on all patterns in self (optionally) with new 
        % wavelength array wl, which will be set to self.wavelength
        % robustopt - struct or name/value pairs with the following options:
        %  > meanEffCutoff - min value for the average efficiency of 
        %  eroded, original, and dilated patterns [default = 0.8]
        %  > effDecCutoff - (positive number) max decrease in average
        %  efficiecny [no default, meanEffCutoff used if not given]. If
        %  given, meanEffCutoff is not checked
        %  > phaseError - (positive number in radians) maximum deviation
        %  from the phase for any wavelength [default = pi/16]
        %  > [NOT SUPPORTED] robmag - magnitude of robustness in nm (how much
        %  eroded/dilated) [default = 10]
        % if none of the robustness parameters are given, simulation is run 
        % with no robustness regardless the option in self.object
        % Note: robustness check is not very well implemented, nor useful
        % here.
%             warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')
            if numel(self) > 1
                self = arrayfun(@(x)simDataset(x,varargin{:}), self, 'UniformOutput', false);
                self = [self{:}]; % for the robustness case, where two Datasets are returned - robust and not
                return
            end
            if mod(length(varargin), 2)
                if isnumeric(varargin{1})
                    wl = varargin{1};
                    self.wavelength = wl;
                    varargin = varargin(2:end);
                end
            else
                wl = self.wavelength;
            end
            if isempty(varargin)
                robust = [];
            else
                defaultopts = struct('meanEffCutoff', 0.8, 'effDecCutoff', [], 'phaseError', pi/16, 'robmag', 10);
                robust = updatestruct(defaultopts, varargin);
            end
            o = self.object;
            if ~isempty(robust)
                if o.numRobust == 1
                    o = set(o, 'robustStartDeviation', [-2 0 2], 'robustEndDeviation', [-5 0 5], 'robustRamp', 30, 'robustWeights', [.5 1 .5]);
                end
            else
                o = set(o, 'robustStartDeviation', 0, 'robustEndDeviation', 0, 'robustWeights', 1);
            end
            I = true(self.numdev, 1); % I == true means robust device, leave in the dataset, otherwise, copy to another one
            robustInd = o.robustStartDeviation == 0;
            [outPhase, outAbsEff, outRelEff] = deal(zeros(length(wl), o.numPol, self.numdev));
            patgen = @self.gp;
            pool_autocreate_flag = Optimizer.enablePool(o.useParfor);
            o.useParfor = 1i; % saves some time in each forwardRun call
            parfor ii = 1:self.numdev
                out = o; % prevent uninitialized variable warning
                ptrn = feval(patgen, ii); %#ok actually required for parfor
                if ~isempty(ptrn) && ~any(isnan(ptrn(:)))
                    out = forwardRun(out, ptrn);
                else
                    continue
                end
                outPhase(:,:,ii) = out.results.phase(:,:,robustInd);
                outAbsEff(:,:,ii) = out.results.absEff(:,:,robustInd); 
                outRelEff(:,:,ii) = out.results.relEff(:,:,robustInd); 
                if ~isempty(robust) % check robustness
                    phaseDiff = unwrap(out.results.phase)-unwrap(out.results.phase(:,:,robustInd));
                    I(ii) = all(abs(phaseDiff) < robust.phaseError);
                    if ~isempty(robust.effDecCutoff)
                        effDiff = mean(abs(unwrap(out.results.absEff(:,:,robustInd))-unwrap(out.results.absEff)), 1);
                        I(ii) = I(ii) & all(effDiff <= robust.effDecCutoff);
                    else
                        meanEff = squeeze(mean(mean(out.results.absEff, 1), 2));
                        I(ii) = I(ii) & all(meanEff >= robust.meanEffCutoff);
                    end
                end
            end
            [self.phase, self.absEff, self.relEff] = deal(outPhase, outAbsEff, outRelEff);
            self = self.prettifyPhase;
            if length(self.wavelength) > 1; self = self.getPhi; end
            Optimizer.enablePool(o.useParfor, pool_autocreate_flag);
            if ~isempty(robust) && ~all(I)
                self = [self.get(I), self.get(~I)];
            end
        end
        
        function self = simGaps(self, wl)
            if nargin == 1; wl = unique([self.wavelength]); end
            for ii = 1:length(self)
                wldiff = setdiff(wl, self(ii).wavelength);
                d = simDataset(self(ii), wldiff);
                [wl, I] = sort([self(ii).wavelength, wldiff]);
                self(ii).phase = cat(1, self(ii).phase, d.phase);
                self(ii).phase = self(ii).phase(I,:,:);
                self(ii).absEff = cat(1, self(ii).absEff, d.absEff);
                self(ii).absEff = self(ii).absEff(I,:,:);
                self(ii).relEff = cat(1, self(ii).relEff, d.relEff);
                self(ii).relEff = self(ii).relEff(I,:,:);
                self(ii).wavelength = wl;
            end
        end
        
        function varargout = calrob(d, upscale_fac, er_rad)
        % Generates and runs forward simulation on three new datasets:
        % upscaled (images with better resolution), eroded and dilated.
        % Later to be compared versus the original Dataset d.
            if nargin < 3; er_rad = 5; end
            if nargin < 2; upscale_fac = 2; end
            [d_upscaled, d_eroded, d_dilated] = deal(d);
            gridScale = @(x) mean(size(x)./d.object.period);
            
            d_upscaled.pattern = cellfun(@(x) upscale(x,upscale_fac), d.gp,'UniformOutput',0);
            d_eroded.pattern = cellfun(@(x) imerode(x, strel('disk',er_rad*gridScale(x),4)), d_upscaled.pattern,'UniformOutput',0);
            d_dilated.pattern = cellfun(@(x) imdilate(x, strel('disk',er_rad*gridScale(x),4)), d_upscaled.pattern,'UniformOutput',0);
            
            d_upscaled = d_upscaled.simDataset;
            d_eroded = d_eroded.simDataset;
            d_dilated = d_dilated.simDataset;
            switch nargout
                case 1
                    varargout = {[d_upscaled, d_eroded, d_dilated]};
                otherwise
                    varargout = {d_upscaled, d_eroded, d_dilated};
            end
        end
        
        function p = gp(self, I)
        % Returns a cell array (or a numerical matrix for a single) of
        % patterns, specified by I (gp = "get pattern"), as 'pattern'
        % property can store data of different types.
        % I - scalar, array of indices, or array of logicals of length
        % equal to the number of devices in the dataset
        % self - has to be a single dataset, for simplicity
        % This is easier and less error-prone than creating a get method
        % for pattern (no subscripting), or overloading subsref for the
        % whole class
            if nargin < 2
                p = self.pattern;
            else
                p = self.pattern(I);
            end
            for ii = 1:length(p)
                switch class(p{ii})
                    case 'double'
                        p{ii} = gen_patt(self.object, p{ii});
                    case 'logical'
                        p{ii} = double(p{ii});
                    case 'function_handle'
                        if nargin(p{ii}) == 0
                            p{ii} = p{ii}();
                        else
                            p{ii} = p{ii}(self.object.imsize);
                        end
                    case 'uint32' 
                        p{ii} = bwunpack(p{ii}, self.object.imsize(1));
                end
            end
            if ii == 1, p = p{1}; end
        end
        
        %------------------------------------------
        %%% End generator methods
        %------------------------------------------
        
        %------------------------------------------
        %%% Plotting methods
        %------------------------------------------
        
        function plot_engine(self, opt, fig, color)
            if nargin < 4; color = []; end
            if nargin < 3; fig = []; end
            if numel(self) > 1
                if isempty(fig)
                    o = [self.object];
                    np = max([o.numPol]);
                    figure
                    for polIter = 1:np
                        ax(polIter) = subplot(1, np, polIter); %#ok
                    end
                    fig = ax;
                end
                % make sure hold is set to 'on'
                switch class(fig)
                    case 'matlab.graphics.axis.Axes'
                        hold(fig, 'on')
                    case {'matlab.ui.Figure', 'double'}
%                         figure(fig)
%                         ax = gca;
%                         hold(ax, 'on')
                end
                arrayfun(@(x) plot(x, fig, color), self);
                return
            end
            assert(~isempty(self.phiCoeffs) || numel(self.wavelength) == 2, ...
                'phiCoeffs property empty. Run getPhi function on the Dataset to calculate them.')
            numPol = self.object.numPol;
            if nargin < 3
                color = [];
            end
            if isempty(fig)
                fig = figure;
            end
            if isa(fig, 'matlab.ui.Figure') || isa(fig, 'double')
                figure(fig)
                for polIter = 1:numPol
                    ax(polIter) = subplot(1, numPol, polIter); %#ok
                    hold(ax(polIter), 'on')
                end
            elseif isa(fig, 'matlab.graphics.axis.Axes')
                ax = fig;
                if numPol > 1
                    ax = repelem(ax, numPol);
                end
            end
            if isempty(color)
                if ~isempty(self.absEff)
                    color = reshape(mean(self.absEff, 1), size(self.absEff,2), size(self.absEff,3));
                    ttl = 'Colorbar: Average Efficiency';
                else
                    color = ['b';'r'];
                    ttl = '';
                end
            elseif strcmp(color,'ff')
                    ff = @(pat)sum(pat(:))/numel(pat);
                    color = repmat(cellfun(ff,self.gp),self.object.numPol, 1);
                    ttl = 'Colorbar: Filling Fraction';
            end
            if size(color,1) < numPol
                color = repmat(color,numPol,1);
            end
            pol_str = {'TE','TM'};
            npol = self.object.npol;
            for polIter = 1:numPol
                switch opt
                    case '2D'
                        if self.object.nwl == 2
                            scatter(ax(polIter),squeeze(mod(self.phase(1,polIter,:)/pi,2)), ...
                                squeeze(mod(self.phase(2,polIter,:)/pi,2)), 56, color(polIter,:), 'LineWidth', 2);
                            xlabel(ax(polIter),['\phi at \lambda=',num2str(self.wavelength(1)),'nm, [\pi]']);
                            ylabel(ax(polIter),['\phi at \lambda=',num2str(self.wavelength(2)),'nm, [\pi]']);
                        else
                            scatter(ax(polIter),squeeze(mod(self.phiCoeffs(3,polIter,:,1)/pi,2)), ...
                                squeeze(self.phiCoeffs(2,polIter,:,1)/pi), 56, color(polIter,:), 'LineWidth', 2);
                            xlabel(ax(polIter),'\phi_0, [\pi]');
                            ylabel(ax(polIter),'\Delta\phi, [\pi]');
                        end
                    case '3D'
                        scatter3(ax(polIter),squeeze(mod(self.phiCoeffs(3,polIter,:,2)/pi,2)),...
                            squeeze(self.phiCoeffs(2,polIter,:,2)/pi), ...
                            squeeze(self.phiCoeffs(1,polIter,:,2)/pi), 56, color(polIter, :), 'LineWidth', 2);
                        xlabel('\phi_0, [pi]');
                        ylabel('\Delta\phi, [pi]');
                        zlabel('\Delta^2\phi, [pi]');
                end
                title(ax(polIter),pol_str{npol(polIter)})
                colormap(ax(polIter),flipud(hot)); colorbar;
                set(ax,'CLim',[0,1]);
                set(ax,'fontsize',16)
                if exist('ttl','var'), sgtitle(ttl); end
            end
        end
        
        function plot(self, varargin)
            plot_engine(self, '2D', varargin{:});
        end
        
        function plot3(self, varargin)
            plot_engine(self, '3D', varargin{:});
        end
        
        function plotErrorPhase(self)
            self = merge(self, [], 1i);
            err = self.get_phase(self.phiCoeffs) - self.phase;
            figure
            subplot(1,2,1)
            histogram(mean(abs(err))/2/pi);
            title('mean error, [2\pi]')
            subplot(1,2,2)
            histogram(max(abs(err))/2/pi);
            title('max error, [2\pi]')
        end
        
        function plotdev(self, I, ax)
            if nargin < 3
                ax = [];
            end
            if nargin == 1
                I = 1;
            end
            for n = 1:numel(self)
                for ii = 1:length(I)
                    patt = self(n).gp(I(ii));
                    plot(self(n).object, patt, ax);
                end
            end
        end
        
        function plotrandgrid(self, varargin)
            if numel(self) > 1
                arrayfun(@(x) plotrandgrid(x, varargin{:}), self);
                return
            end
            I = cellfun(@ischar, varargin);
            if any(I)
                opt = varargin{I};
                varargin = varargin(~I);
            else
                opt = '';
            end
            switch length(varargin)
                case 2
                    [m, n] = varargin{:};
                    total_num = min(m*n, self.numdev);
                case {1, 0}
                    if length(varargin) == 1
                        total_num = min(varargin{1}, self.numdev);
                    else
                        total_num = min(self.numdev, 128);
                    end
                    ratio = 9/16;
                    m = ceil(sqrt(total_num)*ratio)+1;
                    n = ceil(total_num/m);
                otherwise
                    error('Wrong number of arguments');
            end
            if opt % so far only one option, no explicit check required
                I = 1:numdev(self);
            else
                I = randperm(numdev(self));
            end
            figure
            for ii = 1:total_num
                ax = subplot(m,n,ii);
                self.plotdev(I(ii), ax);
                if opt; title(num2str(I(ii))); end
            end
            set(get(gcf, 'Children'), 'FontSize', 18);
        end
        
        function plotpc(self, varargin)
        % plotpc(self, phase, Ny) - "plot phase coverage"
        % generates a pcolor plot where on x axis we have wavelength, on y
        % - phase in units of 2pi from 0 to 1. Blank spaces means that no
        % such phase in the dataset, color - maximum absolute efficiency
        % among all devices that fall in this phase bracket for
        % a particular wavelength.
        % Two possible optional arguments can be provided: phase to be
        % plotted (self.phase by default, but you might want to see the
        % linear fit generated by self.get_phase from self.phiCoeffs) and
        % Ny which is the number of steps in phase between 0 and 2pi. The
        % arguments can be given in any order. Phase can be given either in
        % radians or in the interval [0 1] signifying units of 2*pi 
        % (size = [number of wavelengths, number of entries]
            pol_str = {'TE','TM'};
            npol = self.object.npol;
            for i = 1:length(npol)
                ph = mod(squeeze(self.phase(:,i,:))/2/pi, 1);
                Ny = 320; % 320 shows a nice picture. Setting it to size(ph,2)
                % would give spot to ech device in the dataset
                for ii = 1:length(varargin)
                    if iscalar(varargin{ii})
                        Ny = varargin{ii};
                    else
                        ph = varargin{ii};
                        if max(ph(:)) > 1
                            ph = mod(squeeze(ph)/2/pi, 1);
                        end
                    end
                end
                M = nan(length(self.wavelength), Ny);
                vals = linspace(0,1,Ny+1);
                vals = vals(1:Ny);
                for ii = 1:length(self.wavelength)
                    dif = ph(ii,:).' - vals;
                    [~,I]=min(mod(abs(dif),1),[],2);
                    for jj = 1:Ny
                        if any(I==jj)
                            M(ii,jj) = max(squeeze(self.absEff(ii,1,I==jj)));
                        end
                    end
                end
                figure
                pcolor(self.wavelength, vals, M.');
                axis([self.wavelength([1 end]) 0 1])
                shading flat
                set(gca,'CLim',[0 1]);
                colormap(flipud(hot));
                set(gca,'fontSize',16)
                xlabel '\lambda, nm'
                ylabel 'phase, [2\pi]'
                title(pol_str{npol(i)})
            end
        end
        
        function plotpc2(d, delta_phi)
        % plotpc2(d) - "plot phase coverage 2"
        % See plotpc for general idea. The difference here is that the
        % dataset d here is supposed to be smaller (e.g. what we've
        % selected to use for a deflector or other device). Plots the
        % calculated phase and the ideal phase, corresponding to the linear
        % approximation of these elements.
        % delta_phi is like in the output of getopt dev: if it's length 3,
        % stores coefficients delta^2 phi, delta phi, phi zero shift for
        % quadratic fit and only the last two for linear. Phi zero shift is
        % a number between 0 and 1 that specifies the fraction of the phi
        % zero step (2*pi/number of discrete elements) by which our
        % sampling position is shifted with respect to the lower boundary
        % of the interval.
        % It is optional for linear fit, in which case 0.5 is used (middle
        % of the interval).
            pol_str = {'TE','TM'};
            npol = d.object.npol;
            for i = 1:length(npol)
                if nargin < 2
                    delta_phi = [0, mean(d.phiCoeffs(2,i,:,1)), 0.5];
                end
                switch numel(delta_phi)
                    case 1
                        delta_phi = [0 delta_phi 0.5]; %#ok
                    case 2
                        delta_phi = [0 delta_phi(:).'];
                end
                step_dp = numdev(d);
                coeffs = [repmat(delta_phi(1),1,step_dp); ...
                          repmat(delta_phi(2),1,step_dp); ...
                          ((1:step_dp)-delta_phi(3))*2*pi/step_dp];
                coeffs = reshape(coeffs, 3, 1, step_dp);
                ph_ideal = d.get_phase(coeffs);
                figure, hold on
                h1 = plot(d.wavelength,mod(squeeze(ph_ideal(:,1,:))/2/pi,1),'ro','MarkerSize',5,'MarkerFaceColor','r');
                h2 = plot(d.wavelength,mod(squeeze(d.phase(:,i,:))/2/pi,1),'bo','MarkerSize',5,'MarkerFaceColor','b');
                legend([h1(1), h2(1)], {'fit', 'calc'})
                set(gca,'fontSize',16)
                xlabel '\lambda, nm'
                ylabel 'phase, [2\pi]'
                title(pol_str{npol(i)})
            end
        end

        %------------------------------------------
        %%% End plotting methods
        %------------------------------------------
        
        %------------------------------------------
        %%% Additinoal methods
        %------------------------------------------
        
        function self = set(self, varargin)
        % A standard-looking set method when you pass properties and their
        % values to set as a Name/Value pair. Could be useful if a property's
        % SetAccess is set to private or protected.
            if numel(self) > 1
                self = arrayfun(@(x) set(x, varargin{:}), self);
                return
            end
            fields_given = varargin(1:2:end);
            vals_given = varargin(2:2:end);
            for ii = 1:length(fields_given)
                self.(fields_given{ii}) = vals_given{ii};
            end
        end
        
        function self = convertClass(self, classname)
        % might be useful for dynamic class convertion when multiple
        % subclasses, but probably not the best way to do it
            eval(['self = ',classname,'(self);']);
        end
        
        function write(self, file)
        % Inputs: 
        %   self - Dataset array which will be written (or appended)
        %   file - filename, or matfile object to write into
        % The file will have two variables: d (the Dataset object) and
        % isWriting - a flag that shows that the file is being written to.
        % The function attempts to write 10 times with a random wait of
        % beetween 1 and 10 seconds beteween attemps. It might fail if some
        % other process takes too long to write to the file, or (more
        % probably) something went wrong and the file is either corruped or
        % isWriting is stuck at 1, when it should've been 0.
        
            if ischar(file) && ~isempty(file)
                if length(file) < 4 || ~strcmp(file(end-3:end), '.mat')
                    file = [file, '.mat'];
                end
                if isfile(file)
                    file = matfile(file, 'Writable', true);
                else
                    isWriting = 0;
                    d = self;
                    save(file, 'd', 'isWriting', '-v7')
                    return
                end
            elseif ~isa(file, 'matlab.io.MatFile')
                error('Wrong file parameter')
            end
            for tries = 1:10
                if file.isWriting
                    pause(rand*10);
                else
                    try
                        file.isWriting = 1;
                        file.d = file.d + self;
                    catch ME
                        file.isWriting = 0;
                        rethrow(ME);
                    end
                    file.isWriting = 0;
                    break
                end
            end
        end
        
        function self = circshift(self,K)
            self.phase = circshift(self.phase,K,3);
            self.absEff = circshift(self.absEff,K,3);
            self.relEff = circshift(self.relEff,K,3);
            self.phiCoeffs = circshift(self.phiCoeffs,K,3);
            self.pattern = circshift(self.pattern,K);
            self.hash = circshift(self.hash,K);
        end
        
        function self = bwpack(self) 
            for ii = 1:numel(self)
                self(ii).pattern = cellfun(@(x) bwpack(round(x)), self(ii).gp, 'UniformOutput', false);
            end
        end
        
    end % methods
    
    methods (Static = true)
        function d = load(file) 
        % Loads Dataset object(s) into array d form file(s) provided as an
        % argument string.
            if ischar(file)
                if length(file) < 4 || ~strcmp(file(end-3:end), '.mat')
                    file = [file, '.mat'];
                end
                files = cellfun(@(x)split(ls(x)), file, 'UniformOutput', false);
                files = [files{:}];
                files = files(~cellfun(@isempty, files));
                d = []; % some files can have one dataset, some their arrays,
                % we don't know the length in advance
                for ii = 1:length(files)
                    f = load(files{ii});
                    % find all Dataset instances in f struct
                    isDataset = cellfun(@(x) isa(x, 'Dataset'), struct2cell(f));
                    fnames = fieldnames(f);
                    fnames = fnames(isDataset);
                    for jj = 1:length(fnames)
                        d = [d f.(fnames{jj})]; %#ok - having a cell array doesn't appear to be much quicker
                    end
                end
            elseif isa(file, 'matlab.io.MatFile')
                d = file.d;
            else
                error('Wrong file parameter')
            end
        end
    end
    
end