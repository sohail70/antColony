%%
clear;
close all;
clc;
G = [1	2   0.2;
    1	4	0.1;
    2	5   0.1;
    4  7   0.2;
    5  8   0.8;
    6  11  0.1
    7  2   0.1;
    9  3   0.7;
    10 5   0.8;
    2  6   0.3;
    2  10  0.1;
    2  9   0.2
    11 8   0.3;
    6  8   0.3;
    1  12  0.1;
    12 13  0.1;
    1  14  0.1;
    14 15  0.1;
    1  16  0.1;
    16 17  0.1;
    ];
g = digraph(G(:,1),G(:,2),G(:,3));
p=plot(g,'EdgeLabel',g.Edges.Weight);
% [path1,d,edgepath] = shortestpath(g,1,8,'Method','positive');
% highlight(p,'Edges',edgepath); % 'EdgeColor','r'


%% ACO params
% maxIt = 100;nAnt = 5;alpha = 1;beta = 1;rho = 0.2;Q = .1;%------> answer [1 2 5 8] 
% maxIt = 100;nAnt = 10;alpha = 1;beta = 0;rho = 0.2;Q = .1;%------> answer [1 2 6 8] or [1 2 6 11 8]

maxIt = 100;
nAnt = 10;

alpha = 1;
beta = 0;

rho = 0.2;
Q = .1;
tau0 =10*Q/(sum(g.Edges.Weight)/size(g.Nodes,1)); 
%% Init
bestCost = zeros(maxIt,1);
emptyAnt.Path = [];
emptyAnt.Cost = [];
ant = repmat(emptyAnt,nAnt,1);

g.Edges.eta(:,1) = 1./g.Edges.Weight(:,1);  
g.Edges.tau(:,1) = tau0; %pheromone
startNode = 1;
goalNode = 8;

bestAnt.Cost = inf;
store = [];
%% Main Loop
for it = 1:maxIt
    %move ants
    for k=1:nAnt
        curNode = startNode;
        ant(k).Path = curNode;
        cost = 0;
        while curNode~=goalNode
            suc = successors(g,curNode);
            index=findedge(g,curNode,suc); 
            P = g.Edges.tau(index).^alpha .* g.Edges.eta(index).^beta;
            P = P/sum(P) ; 
            
            % roulettewheelSelection
            r=rand;
            C=cumsum(P);
            j=find(r<=C,1,'first');
            
            ant(k).Path =[ant(k).Path , suc(j)];
            cost = cost + g.Edges.Weight(index(j));
            curNode = suc(j);
            if isempty(curNode) %Deadends!
                cost = inf;
                break;
            end
        end
        ant(k).Cost = cost;
        if ant(k).Cost < bestAnt.Cost
            bestAnt = ant(k);
        end
    end
    
    % Phromone Update
    for k = 1:nAnt
        for i = 1:size(ant(k).Path,2)-1
            index=findedge(g,ant(k).Path(i),ant(k).Path(i+1));
            g.Edges.tau(index) = g.Edges.tau(index) + Q/ant(k).Cost;
        end
    end
    %Evaporation
    g.Edges.tau(:) = (1-rho).*g.Edges.tau(:);
    %storing the best cost
    bestCost(it) = bestAnt.Cost;
    %Info
    disp(['Iter' num2str(it) ': Best Cost = ' num2str(bestCost(it))])
end

for i =1:size(bestAnt.Path,2)-1
    in = findedge(g,bestAnt.Path(i),bestAnt.Path(i+1));
    store = [store;in];
    highlight(p,'Edges',in,'EdgeColor','r'); %
end
