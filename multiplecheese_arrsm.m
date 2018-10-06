% AGGRESSIVE-RESPONSE RESIDUAL SPACE MAPPING ALGORITHM
% MULTIPLE-CHEESE CUTTER ILLUSTRATION
%% House keeping
clc; close all;
clearvars;
co = [4 0]; fo =[4 4];
for ic = 1:length(co)
    %% Inits
    % desired fine model volume response
    Raim = 10;
    % candidate designable parameters
    % initial guess; length and width
    % n.b: if c ~= f then there is a misalignment
    % to get the root solution, first calculate without misalignment, c == f
    % then calculate with misalignment.
    l=1; c = co(ic); w_c = 1;
    l_c = moptm_coarse(Raim,[l, c, w_c]);
    %
    l_f = l_c; f1=fo(ic); f2=fo(ic); w_f=1; w_f1=1; w_f2=1;
    
    %% Initial Model Guess
    %coarse model
    R_c = Rcoarse([l_c, c, w_c]);
    %fine model
    R_f = Rfine([l_f, f1, f2, w_f, w_f1, w_f2]);
    Raim = [Raim Raim Raim];
    % display
    fprintf('l:%g\n', l_c)
    fprintf('w:%g\n', w_c)
    fprintf('R_c: %g\n',R_c)
    fprintf('R_f: %g\n',R_f)
    fprintf('Fine aim: %g\n',Raim)
    
    id = 1;
    
    % Parameter Extraction -> new coarse length mapping
    rng default % For reproducibility
    fun_x = @(x)norm(R_f-Rcoarse([l_c+x, c, w_c])); % cost function
    options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton');
    x = fminunc(fun_x,1,options);
    %
    l = l_c + x;
    f = l-l_c; % error vector
    
    % Mapping Jacobian
    h = Inf;
    B = eye(1,1); % unit mapping
    % init. surrogate and residual response
    R_s = [0 0 0]; dR = 0;

    chck = norm(Raim - R_f); % [] store error
    Rf = R_f; % [] store fine model response
    itTab = []; % store parameters per iteration
    %%
    while norm(h) > 1e-4
        %% Prediction, Evaluate fine length
        % Inverse mapping of coarse length;
        h = -(f) / B; % quasi-newton step in fine space
        l_f = l_f + h; % update
        %% Verification
        % coarse model
        R_c = Rcoarse([l_c,c,w_c]);
        % fine model
        R_f = Rfine([l_f, f1, f2, w_f, w_f1, w_f2]);
        Rf = [Rf; R_f];
        % save tab
        itTab= [itTab; [id f B h l_c R_c R_s l_f R_f (chck(id)*100/(sum(Raim)/3))]];
        %% Next Iterate Prediction
        fun_x = @(x)norm(R_f-Rcoarse([l_c+x, c w_c])); % cost function
        options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton');
        x = fminunc(fun_x,1,options); % alternative: fminunc, but slower convergence
        %
        l = l_c + x;
        f = l-l_c; % update error vector
        
        if norm(f) < 1e-4
            % update response residual
            dR = R_f - R_c;
            % optimise coarse length-> minimise surrogate
            l_c = optim_sug(Raim,[l_c, c, w_c],dR);
            f = l - l_c; % update error vector
            % Inverse mapping of coarse length;
            h = -(f) ./ B; % quasi-newton step in fine space
            l_f = l_f + h;
            %
            R_s = Rsurrogate([l_c, c, w_c],dR);
            B = eye(1,1); % constrain B to an identity mapping
        else
            % broyden rank-one update
            B = B + ((f.*h')/(h'.*h));
        end
        %
        % display
        fprintf('\nIteration.%g\n', id)
        fprintf('f: %g\n',f)
        fprintf('l:%g\n', l_c)
        fprintf('R_c: %g\n',R_c)
        fprintf('R_s: %g\n',R_s)
        fprintf('lf:%g\n', l_f)
        fprintf('R_f: %g\n',R_f)
        fprintf('Fine aim: %g\n',Raim)
        %
        chck = [chck; norm(Raim - R_f)]; %#ok<*AGROW>
        itTab= [itTab; [id f B h l_c R_c R_s l_f R_f (chck(id)*100/(sum(Raim)/3))]];
        % stop if limit attractor reached
        % not converging
        if (id >= 2)
            if abs(chck(id) - chck(id-1)) <= 1e-6
                break;
            end
        end
        id = id + 1;
    end
    fprintf('Fine Model Parameters are: l=%g, w=%g, h=%g, volume=%g\n',...
        l_f,w_f,1,R_f(1))
    fprintf('Fine Model Parameters are: l=%g, w=%g, h=%g, volume=%g\n',...
        l_f,w_f1,1,R_f(2))
    fprintf('Fine Model Parameters are: l=%g, w=%g, h=%g, volume=%g\n',...
        l_f,w_f2,1,R_f(3))
    if ic == 1
        chck_opt = chck;
        l_opt = l_f;
    end
end

%% Visualization
figure(1);
% subplot 1
subplot(211)
plot(1:numel(chck),chck,'-.sr','LineWidth',1.25)
grid on;
xlabel('Iteration','Interpreter','latex')
ylabel('Error, $$||R_{f}-{R_{f}}^{\ast}||$$',...
    'FontSize',12,'Interpreter','latex')
title('Multiple Cheese Cutter:Response Residual Aggressive Space Mapping Optimization',...
    'FontSize',10,'Interpreter','latex')
axis([1,inf,min(chck)-0.01,inf])
% subplot 2
subplot(212)
% figure(2);
for ix = 1:id
    %     l_err(ix) = norm(l_opt-itTab(ix,3));
    l_err(ix) = norm(chck_opt(end)- chck(ix));%#ok<SAGROW>
end
plot(l_err,'Marker','.','MarkerSize',10,'LineWidth',1.25)
% plot(Rf,'-.o','LineWidth',1.25)
grid on;
xlabel('Iteration','Interpreter','latex')
ylabel('Error, $$||l_{f}-{l_{f}}^{\ast}||$$',...
    'FontSize',12,'Interpreter','latex')
% ylabel('$$R_{f}$$',...
%     'FontSize',12,'Interpreter','latex')
title('Multiple Cheese Cutter: Fine-Model Response',...
    'FontSize',10,'Interpreter','latex')
axis([1,inf,min(l_err)-0.01,inf])
% legend({'Cheese1','Cheese2','Cheese3'})


