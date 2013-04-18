function Result = symDIRECT(fun,x_L,x_U,GLOBAL,PriLev)
format long;
if nargin < 5;
    PriLev = [];
    if nargin < 4
        GLOBAL = [];
        if nargin < 3
            x_U = [];
            if nargin < 2
                x_L = [];
                if nargin < 1
                    fun = [];
                end
            end
        end
    end
end
if isempty(PriLev), PriLev=1; end
if isempty(fun) || isempty(x_L) || isempty(x_U)
    disp(' gblSolve requires at least three nonempty input arguments')
    return;
end
if isempty(GLOBAL)
    T = 50; % Number of iterations
    epsilon = 1e-4; % global/local weight parameter.
    tol = 0.01; % Error tolerance parameter.
else
    if isfield(GLOBAL,'iterations') % Number of iterations
        T = GLOBAL.iterations;
    else
        T = 50;
    end
    if isfield(GLOBAL,'epsilon') % global/local weight parameter
        epsilon = GLOBAL.epsilon;
    else
        epsilon = 1E-4;
    end
    if isfield(GLOBAL,'tolerance') % Convergence tolerance
        tol = GLOBAL.tolerance;
    else
        tol = 0.01;
    end
    if isfield(GLOBAL,'optimalvalue') % Known optimal function value
        optimalvalue = GLOBAL.optimalvalue;
    else
        optimalvalue = 0.;
    end
    if isfield(GLOBAL,'optimalvector') % Known optimal vector
        optimalvector = 0.01*GLOBAL.optimalvector;
    end
end
nFunc=0;
x_L = x_L(:);
x_U = x_U(:);
n = length(x_L); % Problem dimensioResult = gblSolve_(fun,x_L,x_U,GLOBAL,PriLev)n
tolle2 = 1E-12;
%
% STEP 1, Initialization
%
if isfield(GLOBAL,'C') && ~isempty(GLOBAL.C)
    % Restart with values from previous run.
    F = GLOBAL.F;
    m = length(F);
    if PriLev > 0
        fprintf('\n Restarting with %d sampled points from previous run\n',m);
    end
    D = GLOBAL.D;
    L = GLOBAL.L;
    d = GLOBAL.d;
    d_min = GLOBAL.d_min;
    f_min = min(F);
    E = max(epsilon*abs(f_min),1E-8);
    [dummy i_min] = min( (F - f_min + E)./D );
    % Must transform Prob.GLOBAL.C back to unit hypercube
    for i = 1:m
        C(:,i) = ( GLOBAL.C(:,i) - x_L )./(x_U - x_L);
    end
else
    % No restart, set first point to center of the unit hypercube.
    m = 1; % Current number of rectangles
    C = ones(n,1)./2; % Matrix with all rectangle centerpoints
    % All C_coordinates refers to the n-dimensional hypercube.
    x_m = x_L + C.*(x_U - x_L); % Transform C to original search space
    if (nargout(fun)>1)
        [logf_min, f_min] = fun(x_m);  % Logaritmic function log(1+f(x)) and function f(x) values at x_m
    else
        f_min = feval(fun, x_m); % Function value at x_m
    end
    nFunc=nFunc+1;
    i_min = 1; % The rectangle which minimizes (F - f_min + E)./D where
    % E = max(epsilon*abs(f_min),1E-8)
    L = ones(n,1)./2; % Matrix with all rectangle side lengths in each dimension
    D = sqrt(sum(L.^2)); % Vector with distances from centerpoint to the vertices
    d = D; % Row vector of all different distances, sorted
    if (nargout(fun)>1)
        F = [logf_min]; % Vector with logaritmic function - log(1 + f(x)) values
        FF = [f_min]; % Vector with function f(x) values
        d_min = logf_min; % Row vector of minimum logaritmic function value for each distance
    else
        F = [f_min]; % Vector with function values
        d_min = f_min; % Row vector of minimum function value for each distance
    end
end
% ITERATION LOOP
t = 1; % t is the iteration counter
if (optimalvalue == 0.0)
    PE=100*f_min; % Percent error
else
    PE=100*(f_min - optimalvalue)/abs(optimalvalue); % Percent error
end
stopRule=0; % Logical stop rule
% Only for Clustering problems
%{
if isfield(GLOBAL,'optimalvector') % Known optimal vector
    PE2 = norm(C(:,1) - optimalvector,Inf);
end
%}

%while PE > tol 
while t <= 4
    %
    % STEP 2 Identify the set S of all potentially optimal rectangles
    %
    S = []; % Set of all potentially optimal rectangles
    S_1 = [];
    idx = find(d==D(i_min));
    %idx = find(abs( d-D(i_min) ) <= tolle );
    if isempty(idx)
        if PriLev >= 0
            fprintf('\n WARNING: Numerical trouble when determining S_1\n');
        end
        return;
    end
    
    for i = idx : length(d)
        idx2 = find( (abs( F-d_min(i) ) <= tolle2  ) & ( D==d(i) ) );
        %idx2 = find( ( F==d_min(i) ) & ( D==d(i) ) );
        S_1 = [S_1 idx2];
        %S_1 = [S_1 idx2(1)]; % Without repeats equals values
    end
   
    %{
    for i = idx : length(d) % PLOR Algorithm
        if (i==idx)||(i==length(d))
            idx2 = find( (abs( F-d_min(i) ) <= tolle2  ) & ( D==d(i) ) );
            S_1 = [S_1 idx2];1.0
            %S_1 = [S_1 idx2(1)]; % Without repeats equal values
        end
    end
    %}
    
    % S_1 now includes all rectangles i, with D(i) >= D(i_min)
    % and F(i) is the minimum function value for the current distance.
    % Pick out all rectangles in S_1 which lies below the line passing through
    % the points: ( D(i_min), F(i_min) ) and the lower rightmost point.
    S_2 = [];
    if length(d)-idx > 1
        a1 = D(i_min);
        b1 = F(i_min);
        a2 = d(length(d));
        b2 = d_min(length(d));
        % The line is defined by: y = slope*x + const
        slope = (b2-b1)/(a2-a1);
        const = b1 - slope*a1;
        for i = 1 : length(S_1)
            j = S_1(i);
            if F(j) <= slope*D(j) + const + 1e-08% + tolle2
                S_2 = [S_2 j];
            end
        end
        % S_2 now contains all points in S_1 which lies on or below the line
        % Find the points on the convex hull defined by the points in S_2
        xx = D(S_2);
        yy = F(S_2);
        h = conhull(xx,yy); % conhull is an internal subfunction
        S_3 = S_2(h);
    else
        S_3 = S_1;
    end
    S = S_3;

    % % % % % % % % % % % % % % (VERSION 1: a part) % % % % % % % % % % % % % % %%
    % Graphical visualization of DIRECT algorithm on separate subplot
    % DEFAULT: t = 4 iterations
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %{
    if PriLev > 2
        subplot(2,2,t,'FontSize',18.0);
        axis([0.0 1.0 0.0 1.0]);
        hold on;
      
        % Draw potentially optimal rectangles
        for jj = 1:length(S)
            j = S(jj);
            rectangle('Position',[C(1,j)-L(1,j),C(2,j)-L(2,j),2*L(1,j),2*L(2,j)],'LineWidth',1.0,'EdgeColor','b','FaceColor','y');
        end
        
        % Draw all rectangles
        for kk = 1:m
            rectangle('Position',[C(1,kk)-L(1,kk),C(2,kk)-L(2,kk),2*L(1,kk),2*L(2,kk)],'LineWidth',1.0,'EdgeColor','b');
            plot(C(1,kk),C(2,kk),'.b','MarkerSize',18.0);
            text(C(1,kk),C(2,kk)+0.05,...
                ['',num2str(F(kk),'%10.2f')],'HorizontalAlignment','center','FontSize',16);
            
            text(C(1,kk),C(2,kk)-0.05,...
                        ['\color{blue}(',num2str(C(1,kk),'%6.2f'),',',num2str(C(2,kk),'%6.2f'),')'],'HorizontalAlignment','center','FontSize',14);
        end
        
        set(gca,'XTick',[0 0.33 0.66 1.0]);
        set(gca,'YTick',[0 0.33 0.66 1.0]);
        % Add line: x_1 = x_2 = ... = x_n
        %subplot(2,2,1);
        %line([0,1],[0,1],'LineWidth',1,'LineStyle','-.')
        %set(gcf,'PaperPositionMode','auto');
        %saveas(gcf,'symDIRECT-Potent-Opt-Rect.eps','psc2');
    end
    %}
    

    % % % % % % % % % % % % % % (VERSION 2) % % % % % % % % % % % % % % % %
    % Simplified graphical visualization of symDIRECT algorithm: 
    % without FaceColor, funct. values. Suitable to show final partitioning
    % Draw only specified iteration (t) together will below part
    % On left side - partitioning; On right side - potent. opt. rectangles
    % In current situation used for comparison DIRECT vs symDIRECT
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %{
    if PriLev > 2 && t == 4
        subplot(1,2,2,'FontSize',18.0);
        hold on;
        axis([-0.0001 1.0001 -0.0001 1.0001]);
      
        % Draw all rectangles 
        for kk = 1:m
            rectangle('Position',[C(1,kk)-L(1,kk),C(2,kk)-L(2,kk),2*L(1,kk),2*L(2,kk)],'LineWidth',1.0,'EdgeColor','b');
            plot(C(1,kk),C(2,kk),'.b','MarkerSize',6.0);
        end
        
        set(gca,'XTick',[0 0.33 0.66 1.0]);
        set(gca,'YTick',[0 0.33 0.66 1.0]);
        %subplot(1,2,2);
        %line([0,1],[0,1],'LineWidth',1,'LineStyle','-.')
        %set(gcf,'PaperPositionMode','auto');
        %saveas(gcf,'DIRECT-vs-symDIRECT.eps','psc2');
    end
    %}
      
    % STEP 3, 5 Select any rectangle j in S
    for jj = 1:length(S) % For each potentially optimal rectangle
        j = S(jj);
        %
        % STEP 4 Determine where to sample within rectangle j and how to
        % divide the rectangle into subrectangles. Update f_min
        % and set m=m+delta_m, where delta_m is the number of new
        % points sampled.
        % 4:1 Identify the set I of dimensions with the maximum side length.
        % Let delta equal one-third of this maximum side length.
        max_L = max(L(:,j));
        I = find( L(:,j)==max_L );
        % I = find( abs( L(:,j) - max_L ) < tolle);
        delta = 2*max_L/3;
        % 4:2 Sample the function at the points c +- delta*e_i for all
        % i in I.
        w=[];
        c = [];  % Centerpoint for new rectangle
        fv = []; % Function value at centerpoints
        log_fv = []; % Logaritmic function value at centerpoints
        for ii = 1:length(I) % for each dimension with maximum side length
            i = I(ii);
            e_i = [zeros(i-1,1);1;zeros(n-i,1)];

            c(:,2*ii-1) = C(:,j) + delta*e_i; % Centerpoint for new rectangle
            x_m1 = x_L + c(:,2*ii-1).*(x_U - x_L); % Transform c to original search space

            if (nargout(fun)>1)
                [log_fv(2*ii-1), fv(2*ii-1)] = fun(x_m1); % Logaritmic function and function values at x_m1 
            else
                fv(2*ii-1) = feval(fun, x_m1); % Function value at x_m1
            end
            nFunc=nFunc+1;
            
            PE_t=100*(fv(2*ii-1) - optimalvalue)/abs(optimalvalue);
            if (PE_t < tol)&&(stopRule == 0)
               nFuncEval=nFunc;
               stopRule = 1;
            end

            c(:,2*ii) = C(:,j) - delta*e_i; % Centerpoint for new rectangle
            x_m2 = x_L + c(:,2*ii).*(x_U - x_L);

            if (nargout(fun)>1)
                [log_fv(2*ii), fv(2*ii)] = fun(x_m2); % Logaritmic function and function values at x_m1 
            else
                fv(2*ii) = feval(fun, x_m2); % Function value at x_m1
            end
            nFunc=nFunc+1;

            PE_t=100*(fv(2*ii) - optimalvalue)/abs(optimalvalue);
            if (PE_t < tol)&&(stopRule == 0)
               nFuncEval=nFunc;
               stopRule = 1;
            end
            w(ii) = min(fv(2*ii-1),fv(2*ii));        
        end
        % 4:3 Divide the rectangle containing C(:,j) into thirds along the
        % dimension in I, starting with the dimension with the lowest
        % value of w(ii)
        [a b] = sort(w);
        for ii = 1:length(I)
            i = I(b(ii));
            L(i,j) = delta/2; 
            D(j) = sqrt(sum(L(:,j).^2));
               
            % Generate all rectangle vertices in terms of O and 1
            ok = 1; % true
            b_ = zeros(n,1);       
            k_ = 1;
            V(:,k_) = b_;
            k_ = k_ + 1;
            while (ok == 1)
                j_ = n;
                while ((j_>0)&&(b_(j_) == 1))
                    b_(j_) = 0;
                    j_ = j_ - 1;
                end
                if (j_ > 0)
                    b_(j_) = 1;
                    V(:,k_) = b_';
                    k_=k_+1;
                else
                    ok = 0; 
                end
            end
            % Transform to real coordinates
            for i_=1:2^n
                for j_=1:n
                    if V(j_,i_) == 0
                        V(j_,i_) = c(j_,2*b(ii)-1) - L(j_,j);
                    else
                        V(j_,i_) = c(j_,2*b(ii)-1) + L(j_,j);
                    end
                end
            end
            % Checking the symmetry for first rectangle
            delta_ = 0.0; % parameter for avoiding distances between two cluster centers smaller than delta_
            
            discard_1 = true; 
            for i_=1:2^n
                check = 0;
                for j_=1:n-1
                    %if V(j_,i_) - delta_ >= V(j_+1,i_)  % Without x=y border
                    if V(j_,i_) >= V(j_+1,i_)           % With x=y border
                        check = check + 1;
                    end
                end
                if check == n-1
                    discard_1 = false; % Cannot discard rectangle
                    %V(:,:)
                    break;
                end
            end
            %{
            if discard_1 == true
                V(:,:) 
            end
            %}
                                      
            % % % % % % % % % % % % % % (VERSION 1: b part) % % % % % % % % % % % 
            % Graphical visualization of symDIRECT algorithm
            % Using after upper part (VERSION 1a) on the same figures
            % Draw (add) rejected rectangles to the graph for all iteration
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %{
            if PriLev > 2
                if discard_1 == true
                    for tt = t+1:4
                        subplot(2,2,tt,'FontSize',18.0);
                        hold on;

                        rectangle('Position',[V(1,1),V(2,1),2*L(1,j),2*L(2,j)],'LineWidth',1.0,'EdgeColor','r');
                        plot(c(1,2*b(ii)-1), c(2,2*b(ii)-1),'xr','MarkerSize',14);
                    
                        text(c(1,2*b(ii)-1),c(2,2*b(ii)-1)-0.05,...
                            ['\color{red}(',num2str(c(1,2*b(ii)-1),'%6.2f'),',',num2str(c(2,2*b(ii)-1),'%6.2f'),')'],'HorizontalAlignment','center','FontSize',14);
                    end
                end
            end
            %}
            
            % % % % % % % % % % % % % % (VERSION 2) % % % % % % % % % % % % % % % %
            %
            if PriLev > 2
                if discard_1 == true
                    subplot(1,2,2,'FontSize',18.0);
                    hold on;
                    % Draw the first rectangle (works for n = 2 test functions)
                    rectangle('Position',[V(1,1),V(2,1),2*L(1,j),2*L(2,j)],'LineWidth',1.0,'EdgeColor','r');
                    plot(c(1,2*b(ii)-1), c(2,2*b(ii)-1),'xr','MarkerSize',10);
                end
            end
            %  
            
            if discard_1 == false
                m = m + 1;
                L(:,m) = L(:,j);
                D(m) = D(j);
                C(:,m) = c(:,2*b(ii)-1);
                if (nargout(fun)>1)
                    F(m) = log_fv(:,2*b(ii)-1); % Vector with logaritmic function values           
                    FF(m) = fv(:,2*b(ii)-1); % Vector with function values
                else
                    F(m) = fv(:,2*b(ii)-1); % Vector with function values
                end
            end

            
            % Generate all rectangle vertices
            ok = 1; % true
            b_ = zeros(n,1);       
            k_ = 1;
            V(:,k_) = b_;
            k_ = k_ + 1;
            while (ok == 1)
                j_ = n;
                while ((j_>0)&&(b_(j_) == 1))
                    b_(j_) = 0;
                    j_ = j_ - 1;
                end
                if (j_ > 0)
                    b_(j_) = 1;
                    V(:,k_) = b_';
                    k_=k_+1;
                else
                    ok = 0; 
                end
            end
            % Transform to real coordinates
            for i_=1:2^n
                for j_=1:n
                    if V(j_,i_) == 0
                        V(j_,i_) = c(j_,2*b(ii)) - L(j_,j);
                    else
                        V(j_,i_) = c(j_,2*b(ii)) + L(j_,j);
                    end
                end
            end
            % Checking the symmetry for second rectangle
            discard_2 = true; 
            for i_=1:2^n
                check = 0;
                for j_=1:n-1
                    %if V(j_,i_) - delta_ >= V(j_+1,i_) % Without x=y border
                    if V(j_,i_) >= V(j_+1,i_)           % With x=y border
                        check = check + 1;
                    end
                end
                if check == n-1
                    discard_2 = false; % Cannot discard rectangle
                    %V(:,:)
                    break;
                end
            end    
            %{
            if discard_2 == true
                V(:,:) 
            end
            %}
            
            % % % % % % % % % % % % % % (VERSION 1: b part) % % % % % % % % % % % 
            % Graphical visualization of symDIRECT algorithm
            % Using after upper part (VERSION 1a) on the same figures 
            % Draw (add) rejected rectangles to the graph for all iteration
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %{
            if PriLev > 2
                if discard_2 == true
                    for tt = t+1:4
                        subplot(2,2,tt,'FontSize',18.0);
                        hold on;

                        rectangle('Position',[V(1,1),V(2,1),2*L(1,j),2*L(2,j)],'LineWidth',1.0,'EdgeColor','r');
                        plot(c(1,2*b(ii)), c(2,2*b(ii)),'xr','MarkerSize',14);
                    
                        text(c(1,2*b(ii)),c(2,2*b(ii))-0.05,...
                            ['\color{red}(',num2str(c(1,2*b(ii)),'%6.2f'),',',num2str(c(2,2*b(ii)),'%6.2f'),')'],'HorizontalAlignment','center','FontSize',14);
                    end
                end
            end
            %}
            
            % % % % % % % % % % % % % % (VERSION 2) % % % % % % % % % % % % % % % %
            %
            if PriLev > 2
                if discard_2 == true
                    subplot(1,2,2,'FontSize',18.0);
                    hold on;

                    rectangle('Position',[V(1,1),V(2,1),2*L(1,j),2*L(2,j)],'LineWidth',1.0,'EdgeColor','r');
                    plot(c(1,2*b(ii)), c(2,2*b(ii)),'xr','MarkerSize',10);
                end
            end
            %
            
            if discard_2 == false
                m = m + 1;
                L(:,m) = L(:,j);
                D(m) = D(j);
                C(:,m) = c(:,2*b(ii));
                if (nargout(fun)>1)
                    F(m) = log_fv(:,2*b(ii)); % Vector with logaritmic function values           
                    FF(m) = fv(:,2*b(ii));    % Vector with function values
                else
                    F(m) = fv(:,2*b(ii));     % Vector with logaritmic function values
                end
            end
            
            if ii == length(I)   
                % Generate all rectangle vertices
                ok = 1; % true
                b_ = zeros(n,1);       
                k_ = 1;
                V(:,k_) = b_;
                k_ = k_ + 1;
                while (ok == 1)
                    j_ = n;
                    while ((j_>0)&&(b_(j_) == 1))
                        b_(j_) = 0;
                        j_ = j_ - 1;
                    end
                    if (j_ > 0)
                        b_(j_) = 1;
                        V(:,k_) = b_';
                        k_=k_+1;
                    else
                        ok = 0; 
                    end
                end
                % Transform to real coordinates
                for i_=1:2^n
                    for j_=1:n
                        if V(j_,i_) == 0
                            V(j_,i_) = C(j_,j) - L(j_,j);
                        else
                            V(j_,i_) = C(j_,j) + L(j_,j);
                        end
                    end
                end
                % Checking the symmetry for last rectangle
                discard_3 = true; 
                for i_=1:2^n
                    check = 0;
                    for j_=1:n-1
                        % if V(j_,i_) - delta_ >= V(j_+1,i_)
                        if V(j_,i_) >= V(j_+1,i_)
                            check = check + 1;
                        end
                    end
                    if check == n-1
                        discard_3 = false; % Cannot discard rectangle
                        %V(:,:)
                        break;
                    end
                end
                %{
                if discard_3 == true
                    V(:,:) 
                end
                %}
                                
            % % % % % % % % % % % % % % (VERSION 1: b part) % % % % % % % % % % % 
            % Graphical visualization of symDIRECT algorithm
            % Using after upper part (VERSION 1a) on the same figures 
            % Draw (add) rejected rectangles to the graph for all iteration
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %{
            if PriLev > 2
                if discard_3 == true
                    for tt = t+1:4
                        subplot(2,2,tt,'FontSize',18.0);
                        hold on;

                        rectangle('Position',[V(1,1),V(2,1),2*L(1,j),2*L(2,j)],'LineWidth',1.0,'EdgeColor','r');
                        plot(C(1,j),C(2,j),'xr','MarkerSize',14);
                    
                        text(C(1,j),C(2,j)-0.05,...
                        ['\color{red}(',num2str(C(1,j),'%6.2f'),',',num2str(C(2,j),'%6.2f'),')'],'HorizontalAlignment','center','FontSize',14); 
                    end
                end
            end
            %}
                
            % % % % % % % % % % % % % % (VERSION 2: b part) % % % % % % % % % % % % % % % %
            %
            if PriLev > 2
                if discard_3 == true
                    subplot(1,2,2,'FontSize',18.0);
                    hold on;
                    
                    rectangle('Position',[V(1,1),V(2,1),2*L(1,j),2*L(2,j)],'LineWidth',1.0,'EdgeColor','r');
                    plot(C(1,j),C(2,j),'xr','MarkerSize',10);
                end
            end
                %
                
                if discard_3 == true
                    %disp('rado');
                    D(j) = 0;
                    if (nargout(fun)>1)
                        F(j) = 1E6; % Vector with logaritmic function values           
                        FF(j) = 1E6;    % Vector with function values
                    else
                        F(j) = 1E6;     % Vector with logaritmic function values
                    end
                end
                
            end
                     
        end
    end
    % UPDATE:
    if (nargout(fun)>1)
        f_min = min(FF); % Minimum function f(x) value
        logf_min = min(F); % Minimum function log(1 + f(x)) value
        E = max(epsilon*abs(logf_min),1E-8);
        [dummy i_min] = min( (F - logf_min + E)./D );
    else
        [f_min, f_min_ind] = min(F);
        E = max(epsilon*abs(f_min),1E-8);
        [dummy i_min] = min( (F - f_min + E)./D );
    end
    i = 1;
    if (optimalvalue == 0.0)
        PE=100*f_min; % Percent error
    else
        PE=100*(f_min - optimalvalue)/abs(optimalvalue); % Percent error
    end
    % For clustering problems
    %{
    if isfield(GLOBAL,'optimalvector') % Known optimal vector
        PE2 = norm(sort(C(:,f_min_ind),'descend') - optimalvector,Inf);
    end
    %}
    D = roundn(D,-12);
    d = unique(D);
    d_min = [];
    for i = 1:length(d);
        idx1 = find(D==d(i));
        %idx1 = find( abs( D-d(i) ) <= tolle );
        d_min(i) = min(F(idx1));
    end
    if (nargout(fun)>1)
        if PriLev > 1
            fprintf('\n Iteration: %d logf_min: %15.10f f_min: %15.10f Sampled points: %d Percent error: %8.4f',...
                t,logf_min, f_min, nFunc, PE);
        end
    else
        if PriLev > 1
            fprintf('\n Iteration: %d f_min: %15.10f Sampled points: %d Percent error: %8.4f',...
                t,f_min,nFunc,PE);
        end
    end
    t = t + 1;
end % ITERATION LOOP
% SAVE RESULTS
Result.f_k = f_min; % Best function value
Result.Iter = t-1;  % Number of iterations
Result.FuncEv=nFunc; % Number of functions evaluations needed for convergence
CC = [];
for i = 1:m % Transform to original coordinates
    CC = [CC x_L+C(:,i).*(x_U-x_L)];
end
Result.GLOBAL.C = CC; % All sampled points in original coordinates
Result.GLOBAL.F = F; % All function values computed
Result.GLOBAL.D = D; % All distances
Result.GLOBAL.L = L; % All lengths
Result.GLOBAL.d = d;
Result.GLOBAL.d_min = d_min;
% Find all points i with F(i)=f_min
if (nargout(fun)>1)
    idx = find(FF==f_min);
    Result.x_k = CC(:,idx); % All points i with FF(i)=f_min
else
    idx = find(F==f_min);
    Result.x_k = CC(:,idx); % All points i with F(i)=f_min
end

function h = conhull(x,y);
% conhull returns all points on the convex hull, even redundant ones.
%
% conhull is based on the algorithm GRAHAMSHULL pages 108-109
% in "Computational Geometry" by Franco P. Preparata and
% Michael Ian Shamos.
%
% Input vector x must be sorted i.e. x(1) <= x(2) <= ... <= x(length(x)).
%
format long;
x = x(:);
y = y(:);
xyAllRound = roundn([x y], -12);
xyUnique = unique(xyAllRound, 'rows');
x = xyUnique(:,1);
y = xyUnique(:,2);
m = length(x);
if length(x) ~= length(y)
    disp('Input dimension must agree, error in conhull-gblSolve');
    return;
end
if m == 2
    h = [];
    xy = [x y];
    for i = 1:size(xy,1)
        h = [h; find(xy(i,1) == xyAllRound(:,1) & xy(i,2) == xyAllRound(:,2))];
    end
    %h = [1 2];
    return;
end
if m == 1
    h = [];
    xy = [x y];
    for i = 1:size(xy,1)
        h = [h; find(xy(i,1) == xyAllRound(:,1) & xy(i,2) == xyAllRound(:,2))];
    end
    %h = [1];
    return;
end
START = 1;
v = START;
w = length(x);
flag = 0;
h = [1:length(x)]'; % Index vector for points in convex hull
while ( next(v,m)~=START ) | ( flag==0 )
    if next(v,m) == w
        flag = 1;
    end
    a = v;
    b = next(v,m);
    c = next(next(v,m),m);
    %[ x(a) y(a); x(b) y(b); x(c) y(c) ]
    %Det = det([ x(a) y(a) 1 ; x(b) y(b) 1 ; x(c) y(c) 1 ])
    if det([ x(a) y(a) 1 ; x(b) y(b) 1 ; x(c) y(c) 1 ]) >= 0
        leftturn = 1;
    else
        leftturn = 0;
    end
    if leftturn
        v = next(v,m);
    else
        j = next(v,m);
        x = [x(1:j-1);x(j+1:m)];
        y = [y(1:j-1);y(j+1:m)];
        h = [h(1:j-1);h(j+1:m)];
        m=m-1;
        w=w-1;
        v = pred(v,m);
    end
end
h = [];
xy = [x y];
for i = 1:size(xy,1)
    h = [h; find(xy(i,1) == xyAllRound(:,1) & xy(i,2) == xyAllRound(:,2))];
end


function h = conhull2(x,y);
% conhull returns all points on the convex hull, even redundant ones.
%
% conhull is based on the algorithm GRAHAMSHULL pages 108-109
% in "Computational Geometry" by Franco P. Preparata and
% Michael Ian Shamos.
%
% Input vector x must be sorted i.e. x(1) <= x(2) <= ... <= x(length(x)).
%
x = x(:);
y = y(:);
m = length(x);
if length(x) ~= length(y)
    disp('Input dimension must agree, error in conhull-gblSolve');
    return;
end
if m == 2
    h = [1 2];
    return;
end
if m == 1
    h = [1];
    return;
end
START = 1;
v = START;
w = length(x);
flag = 0;
h = [1:length(x)]'; % Index vector for points in convex hull
while ( next(v,m)~=START ) | ( flag==0 )
    if next(v,m) == w
        flag = 1;
    end
    a = v;
    b = next(v,m);
    c = next(next(v,m),m);
    if det([ x(a) y(a) 1 ; x(b) y(b) 1 ; x(c) y(c) 1 ]) >= 0
        leftturn = 1;
    else
        leftturn = 0;
    end
    if leftturn
        v = next(v,m);
    else
        j = next(v,m);
        x = [x(1:j-1);x(j+1:m)];
        y = [y(1:j-1);y(j+1:m)];
        h = [h(1:j-1);h(j+1:m)];
        m=m-1;
        w=w-1;
        v = pred(v,m);
    end
end
function i = next(v,m);
    if v==m
        i = 1;
    else
        i = v + 1;
    end
function i = pred(v,m);
    if v==1
        i = m;
    else
        i = v - 1;
    end


