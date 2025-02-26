function par_data = profile2par(bc_base, vals, vars)
%PROFILE2PAR convert a given combination of vel components into parameters

% search for the velocity components and convert them
[vals, vars] = convert2angl(vals, vars);

% preallocate an array of data objects
par_data = cell(1,numel(vars));

% loop over the variables
for i=1:numel(vars)
    % loop over the different cases
    par_data{i} = projProf(bc_base, vals{i})';
end
end

function par_data = projProf(bc_base, vals)
% preallocate the par_data array
par_data = zeros(size(bc_base,2), size(vals,2));

% define a default s vector
s_vec = linspace(0,1,size(vals,1));

% loop over the columns of vals
for i=1:size(vals,2)
    % initialize the lambda value
    lam = 0.1;
    
    % start a while loop to project the profiles
    while true
        % plot the current profile
         clf; hold on; plot(vals(:,i), s_vec);
        
        % get the projection via L1 stabilization
        [sv, st] = lasso(bc_base, vals(:,i), 'Lambda', lam);
        
        % get the vector of the projection
        plot(bc_base*[st.Intercept; sv(2:end)], s_vec);
        
        % show the current lambda and dof
        text(0.1, 0.9, {['current lambda: ' num2str(lam)], ...
            ['current dof: ' num2str(st.DF)]}, 'Units', 'normalized');
        
        % plot the projected profile
        br = input('do you want to choose another lambda?', "s");
        
        % see what answer we got
        if isempty(br), break;
        else, lam = str2double(br); end
    end
    
    % save the identified projection
    par_data(:,i) = [st.Intercept; sv(2:end)];
end
end

function [vals, vars] = convert2angl(vals, vars)
% check if there are the columns vel, yaw and pitch
[~,ids] = ismember({'c_ax', 'c_circ', 'c_rad'}, vars);

%  if they are not available return
if any(ids==0), return; end

% get the components vel yaw and pitch
c_ax = vals{ids(1)}; c_circ = vals{ids(2)}; c_rad = vals{ids(3)};

% calculate the velocity and angle profiles
vals{end+1} = sqrt(c_ax.^2 + c_circ.^2 + c_rad.^2); 
vals{end+1} = atand(c_circ./c_ax); vals{end+1} = atand(c_rad./c_ax);
vars = horzcat(vars, {'vel', 'yaw', 'pitch'});
end
% 
% function [vel, yaw, pitch] = single_profile2par(c_ax, c_circ, c_rad, lam)
% % calculate the velocity profile
% vel = sqrt(c_ax.^2 + c_circ.^2 + c_rad.^2); vsel = ~isnan(vel);
% yaw = atand(c_circ./c_ax); ysel = ~isnan(yaw);
% pitch = atand(c_rad./c_ax); psel = ~isnan(pitch);
% 
% % scale the lambda values for the different coords
% mval = max([vel, yaw, pitch],[],'all');
% lam1 = lam*(max(vel)/mval); lam2 = lam*(max(yaw)/mval); lam3 = lam*(max(pitch)/mval);
% 
% % define the base onto which the profiles are projected
% n = length(c_ax); [M, Bl, Bu, G, ~] = bc_functions.inletbc_base(n); C_base = [M,Bl,Bu,G];
% 
% % project the velocity profile onto the defined base
% [vp_vec, v_stat] = lasso(C_base(vsel,:), vel(vsel), 'lambda', lam1);    % use L1 regularized minimization
% vel = [v_stat.Intercept; vp_vec(2:end)];                                     
% 
% % project the yaw profile onto the defined base
% [yp_vec, y_stat] = lasso(C_base(ysel,:), yaw(ysel), 'lambda', lam2);	% use L1 regularized minimization
% yaw = [y_stat.Intercept; yp_vec(2:end)];
% 
% % project the pitch profile onto the defined base
% [pp_vec, p_stat] = lasso(C_base(psel,:), pitch(psel), 'lambda', lam3);  % use L1 regularized minimization
% pitch = [p_stat.Intercept; pp_vec(2:end)];
% end
% 
% 
% 
% 
% 
% 
% 
% 
% function v_par = dec_velprofile(vel, M, B, ~)
% % initialize a counter
% c  = 1;
% 
% % loop over all the boundarylayer functions
% for i=1:size(B,2)-1
%     for j=1:size(B,2)-1
%         % get the current building blocks
%         c_base = [M, B(:,i:i+1), flipud(B(:,j:j+1))];
%         
%         % project the velocity onto the current base
%         c_svec = c_base\vel;
%         
%         % save the idetified parameters
%         v_par(c).mean = c_svec(1); v_par(c).lean = c_svec(2);
%         sid = size(M,2) + 1; eid = sid + size(B,2) - 1;
%         v_par(c).up_blayer = [i i+1]/size(B,2); 
%         v_par(c).lo_blayer = [j j+1]/size(B,2);
%         
%         % get the current error and count one up
%         res(c) = norm(vel-c_base*c_svec); c = c + 1;
%     end
% end
% 
% % get the best approximation
% [~,mid] = min(res); v_par = v_par(mid);
% end
% 
% function pos = dec_transl(sc_vec, s_vec)
% 
% % sum the piecewise functions
% intv = sum([reshape(sc_vec(1:end-1),[],1), reshape(sc_vec(2:end)',[],1)],2);
% 
% % get the interval with the highest sum
% [~,mid] = max(abs(intv)); s_new = linspace(s_vec(1), s_vec(end), length(sc_vec));
% 
% % get the position
% pos = dot(sc_vec(mid:mid+1)/sum(sc_vec(mid:mid+1)),s_new(mid:mid+1));
% end
% 
% function sc_vec = enc_transl(pos, s_vec)
% 
% % create hat functions
% hf = eye(length(s_vec));
% 
% % loop over the hat functions
% for i=1:length(s_vec), sc_vec(i) = interp1(s_vec, hf(i,:), pos); end
% 
% end