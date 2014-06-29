classdef ClA_operator
% ClA_OPERATOR Class of linear operator
    properties
        fhA
        fhAt
    end
    methods
        function A_oper = ClA_operator(A, At)
            switch nargin
                case 1
                    if isa(A, 'A_operator')
                        A_oper = A;
                    else
                        error('Input must be a class_A_operator');
                    end
                case 2
                    if isa(A, 'function_handle')
                        A_oper.fhA = A;
                        A_oper.fhAt = At;
                    else
                        error('Input must be a function_handel');
                    end
                otherwise
                    error('Syntax: A = A_operator( A, A_transpose )');
            end
        end
        function At = ctranspose(A)
            At = ClA_operator(A.fhAt, A.fhA);
        end
        function At = transpose(A)
            At = ClA_operator(A.fhAt, A.fhA);
        end
        function r = mtimes(A,x)
            r = A.fhA(x);    
        end
        function s = saveobj(objA)
            s = [];
            s.sA = func2str(objA.fhA);
            s.sAt = func2str(objA.fhAt);
        end
    end
%     methods (Static = true)
%         function obj = loadobj(A)
%             if isstruct(A)
%                 tmpfhA = struct2fh(A.sA);
%                 tmpfhAt = struct2fh(A.sAt);
%                 obj = ClA_operator(tmpfhA, tmpfhAt);
%             else
%                 obj = ClA_operator(A);
%             end
%         end
%     end    
end
% 
% function fun___=struct2fh(sfun___)
% 
%     evalvar(sfun___);
% 
%     fun___ = sfun___.function;
%     clear sfun___;
%     fun___=eval(fun___);
% 
% end
% 
% % nested function
% function evalvar(sfun)
% 
% for n=1:length(sfun.workspace)
%     wn = sfun.workspace{n};
%     vlist=fieldnames(wn);
%     for k=1:length(vlist)
%         val = wn.(vlist{k});
%         assignin('caller',vlist{k}, val);
%     end
% end
% 
% end