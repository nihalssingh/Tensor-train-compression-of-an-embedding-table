load('tt_embedding.mat', 'A');
A = double(A);

% A = rand(5,8,25,25,4,8); % d-dimensional tensor A

% e = 2; % prescribed accuracy epsilon
% G = TTSVD(A,e);
% Arc = Rec(G);
% Fn = norm(A-Arc, "fro");
% S = CompRat(A,G);
% scatter(e, S,20,'filled',"black");
% set(gca,'xscale','log');
% hold on
% xlabel('Error Bound');
% ylabel('Compression Ratio');

for e = [0.1 0.2 0.3 0.4 0.5]
    G = TTSVD(A,e);
    Arc = Rec(G);
    Fn = norm(A-Arc, "fro");
%     disp(Fn);
    S = CompRat(A,G);
    grid on
    scatter(e, S,20,'filled',"black");
%     set(gca,'xscale','log');
    hold on
%     xlim([0.001 100])
%     axis 'auto y'
    xlabel('Error Bound');
    ylabel('Compression Ratio');
end

function Arc = Rec(G)
    d = size(G,2);
    for i=1:d-1
        if i == 1
            Arc = cell2mat(G(1));
        end
        Arc = tensorprod(Arc, cell2mat(G(i+1)), i+1, 1);
    end
end


function G = TTSVD(A,e)
    d = size(size(A),2);
    n = size(A);
    
    % Initialization
    del = e/sqrt(d-1)*norm(A, "fro");
    
    % Temp tensor
    C = A;
    r = [];
    r(1) = 1;
    G = {};
    for k=1:d-1
        C = reshape(C,r(k)*n(k), []);
        t = 2;
        % (Ef > del) &&
        while t <= rank(C)
            [U,S,V] = svds(C,t);
            Cest = U*S*transpose(V);
            E = C - Cest;
%             disp(e);
            Ef = norm(E, "fro");
%             P1 = sprintf('Ef = %d, k = %d, t=%d',Ef,k,t);
%             disp(P1);
%             P2 = sprintf('del = %d, k = %d, t=%d',del,k,t);
%             disp(P2);
            r(k+1) = rank(C);
            if Ef <= del 
                r(k+1) = rank(Cest);
%                 disp(r(k+1));
%                 disp(k+1);
                break;
            end
            t = t+1;
        end
        G{k} = reshape(U,r(k),n(k),[]);
        if k == 1
            G{k} = reshape(cell2mat(G(k)), n(k),[]);
        end
        C = S*transpose(V);
    end
    G{d} = C;
%     disp(G);
end

function S = CompRat(A,G)
    d = size(G,2);
    Gele = 0;
    for ii=1:d
        Gele = Gele + numel(cell2mat(G(ii)));
    end
    S = Gele/numel(A);
    disp(Gele);
    disp(numel(A));
end