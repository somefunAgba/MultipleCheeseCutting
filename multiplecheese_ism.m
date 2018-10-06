% IMPLICIT SPACE MAPPING ALGORITHM
% SINGLE-CHEESE CUTTER ILLUSTRATION
%% House keeping
clc; close all;
clearvars;
co = [4 2]; fo =[4 4];
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
    chck = norm(Raim - R_f); % [] store error
    Rf = R_f; % [] store fine model response
    itTab = []; % store parameters per iteration
    %%
    while norm(Raim - R_f) > 1
        %% Parameter Extraction -> width
        rng default % For reproducibility
        fun_x = @(x)norm(R_f-Rcoarse([l_c, c, w_c+x])); % cost function
        % options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton')
        % options = optimoptions(@simulannealbnd, 'HybridFcn', 'fminunc')
        % options = optimoptions(@particleswarm, 'SwarmSize',64,'HybridFcn', @fminsearch)
        % fminsearch - derivative free- simplex method.
        [x,val] = fminunc(fun_x,w_c); % alternative:fminsearch, or fminunc, but slower convergence
        %
        w_c = w_c+x;
        %% Verification 1
        % coarse model
        R_c = Rcoarse([l_c,c,w_c]);
        % fine model
        R_f = Rfine([l_f, f1, f2, w_f, w_f1, w_f2]);
        % display
        fprintf('\nIteration.%g\n', id)
        fprintf('q1: %g\n',x)
        fprintf('l:%g\n', l_c)
        fprintf('w:%g\n', w_c)
        fprintf('R_c: %g\n',R_c)
        fprintf('R_f: %g\n',R_f)
        fprintf('Fine aim: %g\n',Raim)
        itTab= [itTab; [id x l_c w_c R_c R_f (chck(id)*100/(sum(Raim)/3))]];
        %% Prediction -> PE length
        [l_c, xef] = moptm_coarse(Raim,[l_c, c, w_c]);
        l_f = l_c;
        %% Verification 2
        % coarse model
        R_c = Rcoarse([l_c, c, w_c]);
        % fine model
        R_f = Rfine([l_f, f1, f2, w_f, w_f1, w_f2]);
        Rf = [Rf; R_f];
        % display
        fprintf('l:%g\n', l_c)
        fprintf('w:%g\n', w_c)
        fprintf('R_c: %g\n',R_c)
        fprintf('R_f: %g\n',R_f)
        fprintf('Fine aim: %g\n',Raim)
        
        %% Iterate
        chck = [chck; norm(Raim - R_f)]; %#ok<*AGROW>
        itTab = [itTab; [id x l_c w_c R_c R_f (chck(id+1)*100/(sum(Raim)/3))]];
        %% Stop if limit attractor reached
        % not converging
        if (id >= 2)
            %         if abs(itTab(id,2:end) - itTab(id-1,2:end)) <= 1e-2
            if abs(chck(id) - chck(id-1)) <= 1e-9
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
title('Multiple Cheese Cutter: Implicit Space Mapping Optimization',...
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

