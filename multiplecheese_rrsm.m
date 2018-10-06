% RESPONSE RESIDUAL SPACE MAPPING ALGORITHM
% MULTIPLE-CHEESE CUTTER ILLUSTRATION
%% House keeping
clc; close all;
clearvars;
co = [0 2]; fo =[0 0];
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
    
    %% Initial 
    % coarse model
    R_c = Rcoarse([l_c, c, w_c]);
    % fine model
    R_f = Rfine([l_f, f1, f2, w_f, w_f1, w_f2]);   
    % residual
    dR = R_f - R_c;
    R_s = Rsurrogate([l_c,c,w_c],dR);
    % aim
    Raim = [Raim Raim Raim];
    % display
    fprintf('\nIteration.%g\n', 0)
    fprintf('l:%g\n', l_c)
    fprintf('w:%g\n', w_c)
    fprintf('R_c: %g\n',R_c)
    fprintf('R_f: %g\n',R_f)
    fprintf('Fine aim: %g\n',Raim)
    
    id = 1;
    itTab = []; % store parameters per iteration
    chck = norm(Raim - R_f); % [] store error
    itTab = [itTab; [0 dR l_c w_c R_c R_s R_f (chck(id)*100/(norm(Raim)))]];
    Rf = R_f; % [] store fine model response

    %% SM - Iteration
    while norm(Raim - R_f) > 0.001
        %% Optimise Surrogate, Prediction -> length
        l_c = optim_sug(Raim,[l_c, c, w_c],dR);
        l_f = l_c;
        %% Fine Model Verification / Evaluation
        % coarse model
        R_c = Rcoarse([l_c, c, w_c]);
        % fine model
        R_f = Rfine([l_f, f1, f2, w_f, w_f1, w_f2]);
        Rf = [Rf; R_f];
        % residual
        dR = R_f - R_c;
        R_s = Rsurrogate([l_c,c,w_c],dR);
        % display   
        fprintf('\nIteration.%g\n', id)
        fprintf('l:%g\n', l_c)
        fprintf('w:%g\n', w_c)
        fprintf('R_s: %g\n',R_s)
        fprintf('R_c: %g\n',R_c)
        fprintf('R_f: %g\n',R_f)
        fprintf('Fine aim: %g\n',Raim)
        
       %% Termination Condition
        chck = [chck; norm(Raim - R_f)]; %#ok<*AGROW>
        itTab = [itTab; [id dR l_c w_c R_s R_c R_f (chck(id)*100/(norm(Raim)))]];
        % stop if limit attractor reached
        if (id >= 2)
            if abs(chck(id) - chck(id-1)) <= 1e-6
                break;
            end
        end
        id = id + 1;
    end
    % store the root fine model solution
    if ic == 1
        chck_opt = chck;
        l_opt = l_f;
    else
        % display
        formatSpec = 'Fine Model Parameters are: %g, %g, %g (units)\n';
        fprintf(formatSpec,[l_f; w_f; 1])
        formatSpec = 'Fine Model Cheese Volume are: %g, %g, %g (cubic units)\n';
        fprintf(formatSpec,R_f')
    end
    
end

%% Visualization
figure(1);
% subplot 1
subplot(211)
plot(1:numel(chck),chck,'-.sr','MarkerSize',5,'LineWidth',0.75)
grid on;
xlabel('Iteration','Interpreter','latex')
ylabel('$$Norm-2 Error$$',...
    'FontSize',10,'Interpreter','latex')
title('Multiple Cheese Cutter: Response Residual Space Mapping Optimization-2',...
    'FontSize',10,'Interpreter','latex')
axis([1,inf,min(chck)-0.01,inf])
%
subplot(212)
for ix = 1:id
    l_err(ix) = norm(l_opt - itTab(ix,6));%#ok<SAGROW>
end
plot(l_err,'-.','Marker','.','MarkerSize',20,'LineWidth',0.25)
grid on;
xlabel('Iteration','Interpreter','latex')
ylabel('Error, $$||l_{f}-{l_{f}}^{\ast}||$$',...
    'FontSize',10,'Interpreter','latex')
title('Multiple Cheese Cutter: Fine-Model Response',...
    'FontSize',10,'Interpreter','latex')
axis([1,inf,min(l_err)-0.01,inf])
