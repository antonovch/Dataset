classdef StartDataset < Dataset
% Dataset's subclass, whose constructurr is tuned to extract starting geometry
% from Optimizer specifically. It usually resides in Dataset's start
% property.
    methods
        function obj = StartDataset(data)
            if nargin == 0  || builtin('isempty', data)
                data = struct;
            end
            % if it's an Optimizer object, it means we want to record
            % the start pattern only
            if isa(data, 'Optimizer')
                fields = {'patternHistory', 'absEff', 'relEff', 'phase'};
                for n = 1:numel(data)
                    for f = 1:length(fields)
                        val = data(n).(fields{f});
                        if ~isempty(val)
                            if iscell(val)
                                data(n) = set(data(n), fields{f}, val(1,1));
                            else
                                data(n) = set(data(n), fields{f}, val(1,:,:,:)); 
                            end
                        end
                    end
                end
            elseif isa(data, 'Dataset')
                data = dataset2struct(data);
            end
            obj@Dataset(data)
            % removing unneccesary data and prevent nested empty datasets
            obj = set(obj, 'hash', [], 'tgtPhiCoeffs', []);
        end % constructor
        
    end

end