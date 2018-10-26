function NetworkGUI_ThinFilm_V1(~,~)
% A GUI for drawing freestyle fluidics networks and finding the resulting 
% drainage timescale behaviour. The pressure in the drops and the conduit 
% are modelled with the thin-film approximation  



%% To do
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Add transparency to drops and channels to make them look nicer,
% this will require that they don't overlap!
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Allow for intersecting channels. If a new channel is defined that
% overlaps an existing one:
% a) Find the point at which the mid-lines intersect
% b) Save the connections  
% c) Delete the channels
% d) Put a node at the intersection
% e) Define 4 new channels connected to it 
% f) Update everything
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
%% Functions names 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are a lot of functions in this code, here are all their names to
% make them easier to find:
% -------------------------------------------------------------------------
% ScrollScale
% PlaceDropsAndNodes
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Warning management
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Matlab function polyshape in OverlapTest produces warnings when it
% makes the drops, nodes and channels into polygons. I think it can be
% ignored
id = 'MATLAB:polyshape:repairedBySimplify';
warning('off',id)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Positioning axes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a figure to hold the GUI. Remember that position is 
% [left bottom width height] (either in pixels or normaized units)
main = figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6], 'name', 'Network GUI', ...
              'NumberTitle', 'off', 'Color', [.5 .5 .5], 'MenuBar', 'none', 'ToolBar', 'none');
        
% Subplot for the circuit diagram
circuit_plot = Subplot2D([0.02 0.1 0.7 0.8],'','',[],[],[],0.05); hold on; legend off; box on; grid on

% Use the fontsize and colour of this figure
fontSize = get(circuit_plot,'FontSize'); colourT = get(circuit_plot,'XColor');
        
% Add axis in mm with lots of grid lines but no labels or ticks
axis equal
x_limit = 80; y_limit = 50;
xlim([0 x_limit]); ylim([0 y_limit]);

% Set ticks for grid
xticks(0:1:x_limit); yticks(0:1:y_limit);

% Add empty labels, i.e. don't label them but still show the lines
x_labels = cell(1,x_limit); for k = 1:x_limit, x_labels{k} = ''; end; xticklabels(x_labels);
y_labels = cell(1,y_limit); for k = 1:y_limit, y_labels{k} = ''; end; yticklabels(y_labels);

% Make the ticks have zero length
set(circuit_plot,'TickLength',[0, 0])

% Add a line indicating distance between tick marks
x_indicator = [1 11]; y_indicator = [47 47];
line(x_indicator,y_indicator,'LineWidth',1,'color','k');
line([x_indicator(1) x_indicator(1)],[y_indicator(1)-0.4 y_indicator(2)+0.4],'LineWidth',1,'color','k');
line([x_indicator(2) x_indicator(2)],[y_indicator(1)-0.4 y_indicator(2)+0.4],'LineWidth',1,'color','k');

% Label the line
text(mean(x_indicator),mean(y_indicator)+1,'10 mm','HorizontalAlignment','center','Interpreter','Latex','FontSize',fontSize);

% Add functions to the plots which are run when e.g. the mouse is clicked 
set(circuit_plot, 'ButtonDownFcn', @PlaceDropsAndNodes);
set(main, 'WindowButtonMotionFcn', @GridPointIndicator);
set(main, 'WindowscrollWheelFcn', @ScrollScale);
set(main, 'WindowKeyPressFcn', @KeyPress)
set(main, 'WindowKeyReleaseFcn', @DeactivateModifier)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Default physical parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 2.6e-2; rho1 = 1e3; rho2 = 1.85e3; g = 9.81; mu1 = 1e-3;
B = (rho2-rho1)*g/gamma; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% The buttons     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first one makes Add_drop_on true so when the circuit_plot is cliked 
% on a drop is drawn there. The next two do the same thing for nodes and
% drops
uicontrol('Style', 'pushbutton', 'String', 'Add drop', 'units', 'normalized', ...
          'Position', [0.73 0.6 0.12 0.1], 'FontName', 'Times New Roman',...
          'FontSize', 15, 'FontWeight', 'bold','Callback', @Add_drop);
uicontrol('Style', 'pushbutton', 'String', 'Add node', 'units', 'normalized', ...
          'Position', [0.86 0.6 0.12 0.1], 'FontName', 'Times New Roman',...
          'FontSize', 15, 'FontWeight', 'bold','Callback', @Add_node);
uicontrol('Style', 'pushbutton', 'String', 'Add conduit', 'units', 'normalized', ...
          'Position', [0.73 0.5 0.12 0.1], 'FontName', 'Times New Roman',...
          'FontSize', 15, 'FontWeight', 'bold', 'Callback', @Add_channel);
% Clears the current circuit    
uicontrol('Style', 'pushbutton', 'String', 'Clear', 'units', 'normalized', ...
          'Position', [0.86 0.5 0.12 0.1], 'FontName', 'Times New Roman',...
          'FontSize', 15, 'FontWeight', 'bold', 'Callback', @ClearCircuit);
% Runs the ODE solver     
uicontrol('Style', 'pushbutton', 'String', 'Run', 'units', 'normalized', ...
          'Position', [0.73 0.1 0.25 0.1], 'FontName', 'Times New Roman',...
          'FontSize', 15, 'FontWeight', 'bold', 'Callback', @ODEsolve);
% Set the physical parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% uicontrol('Style','text', 'String', 'Valley Position','Position', [400 80 100 20], 'BackgroundColor',[.70 .70 .85]);
% vSpot = uicontrol('Style','edit', 'String', '3','Position', [510 80 50 20],'Callback',@fitnesslandscape);



     


%% Secret options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select which model to use for the pressure and flux
model = 1;

% Show lots of labels in the circuit figure
extra_labels = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







% set the behaviour for when the mouse is over a graphics object (patch)
% pointerBehavior.enterFcn    = @HightlightPatch;
% pointerBehavior.exitFcn     = @NotOver;
% pointerBehavior.traverseFcn = [];
 

 

%% Boolean variables and various other global parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Booleans for which button has been pressed
Add_drop_on = false;
Add_channel_on = false;
Add_node_on = false;
select_on = false;
          
% Variables for counting the number of drops, nodes and channels and cell 
% arrays for storing the data in 
drop_count = 0;
drop_data = cell(0);
node_count = 0;
node_data = cell(0);
channel_count = 0;
channel_data = cell(0);
             
% Cell arrays for holding the pressure functions and flux functions for
% each drop, there size will change 
P_fcn = {};
ODE_RHS = {};

% Maximum and minimum volumes 
min_v = 2.5; max_v = 30;

% Cell arrays for labels on drops and nodes
Vol_label = cell(0);
if extra_labels
    pressure_label = cell(0);
    radius_label = cell(0);
    node_label = cell(0);
    width_label = cell(0);
    length_label = cell(0);
end

% Colour of the circuit and edges
fill_colour = [0.8 0.8 0.9]; edge_colour = [0.1 0.1 0.5];



              channel_clicks = 1;
              channel_ends = [];
%               end_pressures = [];
              channel_arrow = cell(0);
              drop_radius = 0;
              drop_volume = 0;
              
              channel_width = 0;
              channel_length = 0;
              channel_connections = [];
              

             
             vol_scroll = false;
             
             current_object = [];
             
pressure_plots = {};

highlight_plot = [];

colours = [];

theta = linspace(0,2*pi,50);

ODE_RHS = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
      
%% The functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 




%% Detect mouse scroll
% Modify either the parameters of the current object by scrolling the mouse
% wheel and update the image
function  ScrollScale(~, event)       
    
    % If scroll is used while a channel is being defined reset the click
    % count. This is because a non-channel thing may be getting altered
    channel_clicks = 1;
    
    % Only do something if the current object is a drop or channel 
    if ~isempty(current_object)
        
        % Update drops
        if current_object.Tag(1) == 'd'
            
            % Detect scroll and direction and modify either volume or radius
            % depending on whether shift key is held down 
            if event.VerticalScrollCount > 0    
                if vol_scroll == true
                    drop_volume = min(drop_volume+0.5,max_v); 
                else
                    drop_radius = min(drop_radius+0.1,4.2);            
                end
            else 
                if vol_scroll == true
                    drop_volume = max(drop_volume-0.5,min_v); 
                else
                    drop_radius = max(drop_radius-0.1,0.5);            
                end
            end

            % The drops may now have different pressures and so too may the
            % nodes; both of these need to be updated
            UpdateDrop(str2double(current_object.Tag(2:end)),drop_radius,drop_volume);
            NodePressures()
            
            % The function DrawArrow is very slow when using either the 
            % arrow or quiver functions (but it can draw a line quickly). 
            % It seems that drawing arrows is just something that Matlab 
            % can't do quickly, so I'll update the channels using a timer
            % that resets whenever a scroll happens
            t = timer('TimerFcn', 'stat=false; disp(''Timer!'')',... 
                 'StartDelay',10);
            
            % quiver
            % The direction of the initial fluxes may be altered as too may
            % be the length of the channels, so they are all updated. ITS 
            % QUITE LIKELY THAT EVERY SINGLE CHANNEL WONT NEED TO BE
            % UPDATED, WE ONLY NEED TO CHECK THE ONES CONNECTED TO THIS
            % DROP OR CONNECTED TO THIS DROP VIA CHANNELS CONNECTED TO 
            % NODES
            for i0_0 = 1:channel_count  
                UpdateChannel(channel_data{i0_0}{1},channel_data{i0_0}{2}(1)) 
            end
             
        % Update channels
        elseif current_object.Tag(1) == 'c'

            % Detect scroll and direction and modify channel width
            if event.VerticalScrollCount > 0    
                channel_width = min(channel_width+0.1,2); 
            else 
                channel_width = max(channel_width-0.1,0.2); 
            end
        
            % This channel now needs to be updated, but this won't alter
            % the pressures or anything
            UpdateChannel(current_object,channel_width)  

        end
            
        % The pressures in the drop
        % Make separate functions for updating geometry and for updating
        % initial flux directions
        

%         
        % Update the information in the drop data array and change the 
        % figure 
        
                
        % Then recalculate the initial pressures in the nodes if there are
        % any
%         if node_count > 0
%             NodePressures()
%         end
        

        
        % Update the information in the drop data array and change the 
        % figure 
            
            
            
    
                % If some small amount of time has elapsed since the scroll wheel
        % was used we update the solution. Otherwise updating it on each
        % scroll takes too much computation
%         if Update_now == true
%             UpdateDrop()
%             UpdateChannel(current_object)  
%         end

    end
            
%             
%         type = Object.Tag(1);
%         
%             
%             
%     % Remove the true part get object data when it is clicked on or created
%     if plot_drop == true
%         
% 
%         
%         
%         
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
      
        


function KeyPress(~, eventdata)
    
    % Use shift to alternate between changing the volume and radius
    if strcmp(eventdata.Key,'shift')
           vol_scroll = true;
    end
    
    % Use backspace key to delete current object (if it exists), all 
    % channels connected to it need to be deleted also and this data 
    % removed from the arrays 
    if ~isempty(current_object) && strcmp(eventdata.Key,'backspace')         

        if  current_object.Tag(1) == 'c'
            DeleteChannel(current_object)            
        elseif current_object.Tag(1) == 'd' || 'n'
            DeleteHubs(current_object)
        end
        
        % Return the current object variable to its initial state 
        current_object = [];
        
        
        
%         % Find what type the object is
%         type = current_object.Tag(1);
%         pos = str2double(current_object.Tag(2:end));
%         
%         % If it is a node or drop remove connecting channels also
%         if type == 'd' && ~isempty(drop_data{pos}{3})
%             
%             % Indices for connected channels 
%             connections = drop_data{pos}{3};
%             
%             for i = connections
%                 % Delete the channels and remove the data from the array
%                 delete(channel_data{i}{1});
%                 channel_data{i} = [];
%             end
%             
%             % Now make the entry in the data array empty
%             drop_data{pos} = [];
%             
%         elseif type == 'n' && ~isempty(node_data{pos}{3})
%                         
%             % Indices for connected channels 
%             connections = node_data{pos}{3};
%             
%             for i = connections
%                 % Delete the channels and remove the data from the array
%                 delete(channel_data{i}{1});
%                 channel_data{i} = [];
%             end
%             
%             % Now make the entry in the data array empty
%             node_data{pos} = [];
%         end
%         
%         % Delete the graphics object
%         delete(current_object)
        
    end
end


        
        
        
        

function DeactivateModifier(~, eventdata)
    if strcmp(eventdata.Key,'shift')
           vol_scroll = false;
    end
end





function ODEsolve(~, ~)
   
    % Empty vector for initial volumes
    V_initial = zeros(1,drop_count);    
 
    % Empty cell array for the flux functions into each drop and for the
    % pressure functions
    Q_fcn = cell(1,drop_count);
    P_fcn = cell(1,drop_count);
 
    %% Pressures
    % First go through all the drops and find the pressure in each as a 
    % function of volume, loading the data for the thick model if neccesary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if model == 1
        
        for i0 = 1:drop_count
            
            % We need the drop base radius
            base_radius = 1e-3*drop_data{i0}{2}(1);
            
            % Pressure is linear function of volume for the thin-film models 
            beta = B*gamma*besselj(0,sqrt(B)*base_radius) / (pi*base_radius^2*besselj(2,sqrt(B)*base_radius));
            P_fcn{i0} = @(V) beta*V;
        
        end
        
    elseif model == 2
        
        % Load the data, this was scaled using the drop base radius
        thickPressure = load('Pressures');
        
        for i0 = 1:drop_count
            
            % We need the drop base radius
            base_radius = 1e-3*drop_data{i0}{2}(1);
            
            % Contour3 also plots the contour (the function contourc gets just the 
            % data but the input needs to be vectors and doesn't work with the
            % multivalued stuff I have here. So I'll create a figure and delete it,
            % this is STUPID!
     
            stupid = figure; hold on
        
            % We need to extract to curves for the pressure in each drop, we do
            % this by finding the abropriate contour in each drop
            M1 = contour3(thickPressure.V,thickPressure.p,thickPressure.Bo,[B*base_radius^2 B*base_radius^2]);
            
            delete(stupid)
    
            % Extract the data from the contours, the point at the origin seems to
            % be missing so we add this on
            V_full = [0 M1(1,2:end)]; P_full = [0 M1(2,2:end)];

            % Also the output may contain non-unique values, which the interpolator
            % doesn't like. These are removed from the volume vector and the
            % corressponding pressure vector, while we're at it we also
            % dimensionalise the data        
            [~, ind] = unique(V_full);
            V_full = base_radius^3*V_full(ind);
            P_full = gamma/base_radius*P_full(ind);
               
            % Interpolation is then used to find the pressure for any volume (in
            % range) 
            P_fcn{i0} = @(V) interp1(V_full,P_full,V,'spline');
    
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Fluxes into drops
    % For passive pumping the conduits are assumed to always be thin, so
    % the flux is the same function for either model. For each drop we find
    % the total flux into it  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i1 = 1:drop_count
        
        % The initial values for the ODE are the drop volumes 
        V_initial(i1) = 1e-9*drop_data{i1}{2}(4);
            
        % The function is built iteratively so we need to start with an
        % empty function
        Q_fcn{i1} = @(t,V) 0;

        % Find all the channels connected to it and find the flux in them
        for i2 = drop_data{i1}{3}
            
            % Length and half width of the channel  
            a = 1e-3*channel_data{i2}{2}(1)/2;
            L = 1e-3*channel_data{i2}{2}(2);
            
            % The flux is multiplied by a constant in the thin-film models
            sigma = a^7/(105*gamma^3*mu1*L); 
        
            % Find if the other end is a drop or node, so first find which
            % end we are dealing with.
            % If the first entry in the channel connections is this
            % drop we look at the second entry, otherwise we want the
            % first
            str = ['d' num2str(i1)];
            if strcmp(channel_data{i2}{3}{1},str)                
                type = channel_data{i2}{3}{2}(1);
                pos = str2double(channel_data{i2}{3}{2}(2:end));                             
            else                
                type = channel_data{i2}{3}{1}(1);
                pos = str2double(channel_data{i2}{3}{1}(2:end));
            end             

            if type == 'd'  
                % If it's connected to a drop we use the pressure function
                % to find the flux as afunction of volume 
                Q_fcn{i1} = @(t,V) Q_fcn{i1}(t,V) - sigma*(P_fcn{i1}(V(i1))^4 - P_fcn{pos}(V(pos))^4);    

            else
                % Otherwise it's connected to a node and the pressure is an
                % unknown 
                Q_fcn{i1} = @(t,V) Q_fcn{i1}(t,V) - sigma*(P_fcn{i1}(V(i1))^4 - V(drop_count+pos)^4); 
            end                                          
            
        end
    end               
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Fluxes through nodes
    % Next go through all the nodes and find the flux into each of them
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sys_RHS = zeros(1,node_count); % Vector for node pressures
    node_connections_matrix = zeros(node_count); % Matrix for known pressures

    for i3 = 1:node_count
        
        % The function is built iteratively so we need to start with an
        % empty function
        Q_fcn{drop_count+i3} = @(t,V) 0;
        
        % The initial values for the pressure in the node are given by an
        % algebraic equation. The node may be connected to another node
        % (with unknown pressure) so we will to solve a linear system 
        % (I - M)PJ^4 = [sum(Ai*Pi^4)/sum(Ai)]
        % for i in connections, where Ai is -a_i^7/(105*gamma^3*mu1*Li) and
        % M(i,k) is Ai/sum(Ai) if node i is connected to node k and zero
        % otherwise. So here are some vectors for storing 
        % sum(Ai*Pi^4) and sum(Ai) for each node 
        pressure_sum = 0;
        flux_const_sum = 0;
 
        % Find all the channels connected to it and find the flux in them
        for i4 = node_data{i3}{3}
            
            % Length and half width of the channel 
            a = 1e-3*channel_data{i4}{2}(1)/2;
            L = 1e-3*channel_data{i4}{2}(2);
            
            % The flux is multiplied by a constant in the thin-film models
            sigma = a^7/(105*gamma^3*mu1*L); 
            
            % Find if the other end is a drop or node, so first find which
            % end we are dealing with
            % If the first entry in the channel connections is this
            % node we look at the second entry, otherwise we want the
            % first
            str = ['n' num2str(i3)];
            if strcmp(channel_data{i4}{3}{1},str)                 
                type = channel_data{i4}{3}{2}(1);
                pos = str2double(channel_data{i4}{3}{2}(2:end));
            else                
                type = channel_data{i4}{3}{1}(1);
                pos = str2double(channel_data{i4}{3}{1}(2:end));               
            end
            
            if type == 'd'
                % If it's connected to a drop we use the pressure function
                % to find the flux as a function of volume
                Q_fcn{drop_count+i3} = @(t,V) Q_fcn{drop_count+i3}(t,V) - sigma*(V(drop_count+i3)^4 - P_fcn{pos}(V(pos))^4);    
                  
                flux_const_sum = flux_const_sum - sigma; 
                pressure_sum = pressure_sum - sigma * P_fcn{pos}(1e-9*drop_data{pos}{2}(4))^4;
            else
                % Otherwise it's connected to a node and the pressure is an
                % unknown and we need to add an entry to the matrix
                Q_fcn{drop_count+i3} = @(t,V) Q_fcn{drop_count+i3}(t,V) - sigma*(V(drop_count+i3)^4 - V(drop_count+pos)^4);%+ flux_fcn(a,L,R_i,R_e,[],V(drop_count+i4),[],V(drop_count+pos),model);
 
                flux_const_sum = flux_const_sum - sigma;%flux_const(a,L,model);
                node_connections_matrix(i3,pos) = -sigma;%flux_const(a,L,model);                   
            end                  
                
        end
                    
        node_connections_matrix(i3,:) = node_connections_matrix(i3,:)/flux_const_sum;
        sys_RHS(i3) = pressure_sum/flux_const_sum;          
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     
    
     % Initial values for the node pressures
    if node_count > 0
        
        sys_LHS = diag(ones(node_count,1)) - node_connections_matrix;
        node_initial = (sys_LHS\sys_RHS').^(1/4);
        V_initial = [ V_initial node_initial' ];
        
    end
    
    % Definition of the ODE, the function ODE_RHS is a global variable so
    % it can be used in the ZeroFlux events function
    ODE_RHS = @(t,V) Q_fcn{1}(t,V); 
    for i5 = 2:length(Q_fcn)
        
        ODE_RHS = @(t,V)[ ODE_RHS(t,V); Q_fcn{i5}(t,V) ];
        
    end
     
%     V_initial
%     [ P_fcn{1}(V_initial(1)) P_fcn{2}(V_initial(2)) V_initial(3) V_initial(4) ]
%     ODE_RHS(1,V_initial)
%     [Q_fcn{1}(1,V_initial) Q_fcn{2}(1,V_initial) Q_fcn{3}(1,V_initial) Q_fcn{4}(1,V_initial)]

%     efhbefhf
%     V_initial(1)
%     Q_fcn{1}(1,V_initial)
%     sys_LHS
%     V_initial
%     node_initial
%     node_connections_matrix
    
%     V_initial
%      ODE_RHS(1,V_initial)
     
% V_initial
% [Q_fcn{1}(1,V_initial) Q_fcn{2}(1,V_initial) Q_fcn{3}(1,V_initial) Q_fcn{4}(1,V_initial)]


% efefefef
% etteryh
% V_initial = [4.1762    0.0483    4.1762    0.0483 1]
%     VV = linspace(0,4e-8,100);
%     figure
%     plot(VV,    P_fcn{2}(VV))
% 
%     P_fcn{2}
%     
% P_fcn{1}(VV(1))^4 - P_fcn{2}(VV(2))^4
%             
% 
% dgjhf
%     V_initial
% %     ODE_RHS(1,V_initial)
%     Q_fcn{1}(1,V_initial)
%     Q_fcn{2}(1,V_initial)

%     ODE_RHS
  
    % We have to solve algebraic equations for the nodes, so we need to define
    % a mass matrix
    M = diag([ones(1,drop_count) zeros(1,node_count)]);

    options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Mass',M);%,'Events', @ZeroFlux);

    % For a dumbbell set-up the time-scale is 1/epsilon^7. So we find the
    % smallest epsilon (1/2*channel_width/channel_length) and solve on this
    % time-scale 
    epsilon_min = 1; 
    for i7 = 1:channel_count
        % Length and half width of the channel 
        epsilon = channel_data{i7}{2}(1)/2/channel_data{i7}{2}(2);
        epsilon_min = min(epsilon,epsilon_min);
    end
%     channel_data{i7}{2}(1)/2
%     channel_data{i7}{2}(2)
%     epsilon_min
max_time = 8000000;%1/epsilon_min^7;
    [ts,vs] = ode15s(ODE_RHS,[0 max_time],V_initial,options);
    
    
    
    
    
    % For the drops we find the pressure ratio so a vector of pressure can
    % be written. This might be a bit of a silly way to do it, maybe this
    % information should be added to the data array
    ps = vs;
    for i8 = 1:drop_count        
        % Drop radius
        R = 1e-3*drop_data{i8}{2}(1);
        % Ratio fo drop volume to presure
        beta = B*gamma*besselj(0,sqrt(B)*R)/(pi*R^2*besselj(2,sqrt(B)*R));                                    
        % The preesure as a function of time
        ps(:,i8) = beta*ps(:,i8);
    end
    
    % We also need to find the flux in the channels with the arrows
    % indicating the initial direction?
    fs = zeros(length(ts),channel_count);
    size(channel_connections)
    for i8 = 1:channel_count        
        % Find the connected objects for each channel. This is very similar
        % to the code used when the chanel (and arrows) were defined
        P = cell(1,2);
        for i9 = 1:2
%             channel_data{i8}{3}{i9}
% %             sfsf
            type = channel_data{i8}{3}{i9}(1);
            pos = str2double(channel_data{i8}{3}{i9}(2:end));
            if type == 'd'
                P{i9} = ps(:,pos);
            elseif type == 'n'
                P{i9} = ps(:,pos);
            end
        end
        
        % Length and half width of the channel 
        a = 1e-3*channel_data{i8}{2}(1)/2;
        L = 1e-3*channel_data{i8}{2}(2);
            
        % The direction of the flux is determined from the initial
        % condition 
        if P{1}(1) > P{2}(1)
            fs(:,i8) = a^7/(105*gamma^3*mu1*L) * (P{1}.^4 - P{2}.^4);
        else
            fs(:,i8) = a^7/(105*gamma^3*mu1*L) * (P{2}.^4 - P{1}.^4);
        end
            
        
    end
    
    
    % We give different colours to the lines each drop, channel and node 
    colours = jet(node_count+drop_count+channel_count);
    
    
    
ODEPressurePlot(ts,ps,pressure_plot)


ODEVolumePlot(ts,vs(:,1:drop_count),volume_plot)


ODEFluxPlot(ts,fs,flux_plot)

end

 
% Stop the ODE solver when the RHS of the equation is sufficiently small
% (i.e. we are close to steady state). This is done by finding the norm of
% the right hand side, when it is less than 1e-18 we get roughly the right
% behaviour
function [position,isterminal,direction] = ZeroFlux(t,v)
    dt = ODE_RHS(t,v); % Right hand side of the ODE
    position = norm(dt) - 2e-18; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end

%% Functions for plotting the ODE solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the flux 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ODEPressurePlot(ts,ps,ax)
    
    % Set the current axis and clear any existings things from it
    cla(ax);
    subplot(ax); hold on 
    
    % Remove entries from global variable pressure_plots 
    pressure_plots = {};
    
    % Plot the solution for each pressure and update the pointer behaviour
    % for the associated patch. I got this from https://uk.mathworks.com/help/images/ref/iptsetpointerbehavior.html
    % I don't know exactly what it all does! The time-scale is hours
    for i = 1:size(ps,2)
        pressure_plots{i} = plot(ts/3600,ps(:,i),'color',colours(i,:));
%         iptSetPointerBehavior(pressure_plots{i}, pointerBehavior);
%         iptPointerManager(gcf)
    end
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ODEVolumePlot(ts,vs,ax)
    
    % Set the current axis and clear any existings things from it
    cla(ax);
    subplot(ax); hold on 
    
    % Remove entries from global variable volume_plots 
    volume_plots = cell(1,size(vs,2));
    
    % Plot the solution for each pressure and update the pointer behaviour
    % for the associated patch. I got this from https://uk.mathworks.com/help/images/ref/iptsetpointerbehavior.html
    % I don't know exactly what it all does! The time-scale is hors and the
    % volume is measured in micro-litres
    for i = 1:size(vs,2)
        volume_plots{i} = plot(ts/3600,1e9*vs(:,i),'color',colours(i,:));
%         iptSetPointerBehavior(volume_plots{i}, pointerBehavior);
%         iptPointerManager(gcf)
    end
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ODEFluxPlot(ts,fs,ax)
    
    % Set the current axis and clear any existings things from it
    cla(ax);
    subplot(ax); hold on 
    
    % Remove entries from global variable flux_plots 
    flux_plots = cell(1,size(fs,2));
    
    % Plot the solution for each pressure and update the pointer behaviour
    % for the associated patch. I got this from https://uk.mathworks.com/help/images/ref/iptsetpointerbehavior.html
    % I don't know exactly what it all does! The time-scale is hours and
    % the flux is measured in nano-litres per second 
    for i = 1:size(fs,2)
        flux_plots{i} = plot(ts/3600,1e12*fs(:,i),'color',colours(node_count+drop_count+i,:));
%         iptSetPointerBehavior(flux_plots{i}, pointerBehavior);
%         iptPointerManager(gcf)
    end
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












% Function
    function HightlightPatch(fig, ~)
        % This function is called when the pointer enters a patch object in 
        % the circuit diagram.
        
        % The mouse position is used to work out what object the mouse is
        % over
        C = get(fig.CurrentAxes,'CurrentPoint');
        x = C(1,1); y = C(1,2);
        
        % For each of the graphics objects check if this point is in them
%         for i0 = 1:2
           % The pointer function gets activated as the patch is entered so
           % usually gives a coordinate just outside it, so we define a
           % slightly larger polygon and see if the point is in this
           verts = channel_data{1}{1}.Vertices; % Vertices of the patch 
           poly1 = polyshape(verts(:,1)',verts(:,2)'); % Define a new polygon
           [xc,yc] = centroid(poly1); % Find its centre and use this to rescale it 
           poly2 = scale(poly1,2,[xc,yc]); % Make a bigger version (1.1 times bigger)
           verts2 = poly2.Vertices; % Find the new vertices
           
           % Is the mouse in our near this patch
           [in, on] = inpolygon(x,y,verts2(:,1),verts2(:,2));

           subplot(circuit_plot)
           highlight_plot = plot(poly2);
%            plot(pressure_plots{lin_num}.XData,pressure_plots{lin_num}.YData,'linewidth',8,'color','w');
           uistack(highlight_plot,'bottom');

%         end
        
        
        
%         x;
%         y;
%         
%         verts(:,1);
%         verts(:,2);
%         
%         
%         
%         % Make slightly larger polygon 
%         polyin = polyshape(verts(:,1)',verts(:,2)');
%         polyin
%         scale(polyin,2)
% %  poly1 = 
%         poly1.Vertices
        
%         for i0 = 1:channel
        % We use the position of the figure to work out which axis we are
        % in (it can be rescaled). Then the current point of that axis is 
        % used to find which object it is  
% fig        % The output of pointerBehavior is in pixels, so we get the pixel
        % position of the relevant figure
%         

% circuit_pixel_pos = getpixelposition(circuit_plot)



%         
%         currentPoint
%         if  circuit_pixel_pos(1) <= currentPoint(1) && currentPoint(1) <= circuit_pixel_pos(1) + circuit_pixel_pos(3) && ...
%             circuit_pixel_pos(2) <= currentPoint(2) && currentPoint(2) <= circuit_pixel_pos(2) + circuit_pixel_pos(4)        
              
% subplot(circuit_plot)

                
%         end
        
        %         currentPoint
        
        % Current pixel position of the figure and the subplots
%         pos_main = get(main, 'Position')
%         poo = getpixelposition(flux_plot);
%         pos_sub1 = getpixelposition(flux_plot);%get(flux_plot, 'Position');
%         pos_sub2 = getpixelposition(shear_plot);%get(shear_plot, 'Position');
%         
%         if pos_sub1(1) <= currentPoint(1) && currentPoint(1) <= pos_sub1(1) + pos_sub1(3) && ...
%            pos_sub1(2) <= currentPoint(2) && currentPoint(2) <= pos_sub1(2) + pos_sub1(4)
%                     C = get(flux_plot,'CurrentPoint');
%                     subplot(flux_plot)
%         elseif pos_sub2(1) <= currentPoint(1) && currentPoint(1) <= pos_sub2(1) + pos_sub2(3) && ...
%                pos_sub2(2) <= currentPoint(2) && currentPoint(2) <= pos_sub2(2) + pos_sub2(4)
%             C = get(shear_plot,'CurrentPoint');
%             subplot(shear_plot)
%         end
%     
    
%         
        
%           x
% %           y
%           
%           % Find the line that has this point in it. The lines will all
%           % have the same x-values in common, find where this is then
%           % compare the y-values to see which is closest to the current 
%           % point  
%           
%           % Array of lines
%           lines = get(shear_plot,'children');
%           pressure_plots{1}.XData;
%           % The index of x-value nearest to click
%           [~,I] = min(abs(pressure_plots{1}.XData - x));
%           
%           % Create a vector of y-values for this x
%           y_vals = zeros(1,length(pressure_plots));
%           for i = 1:length(pressure_plots)
%               y_vals(i) = pressure_plots{i}.YData(I);% lines(i).YData(I);                             
%           end
%           
% %            
%           % The index of the plotted line 
%           [~,lin_num] = min(abs(y_vals-y));
% %           y_vals
%           lin_num;
%           % Find which line has the closest y-value to the current point
% %           lines(:).YData(1)
%           
%           highlight_plot = plot(pressure_plots{lin_num}.XData,pressure_plots{lin_num}.YData,'linewidth',8,'color','w');
%           uistack(highlight_plot,'bottom');
%           
% %         currentPoint
%         
%         % Use this to 
%         
%         % Gte the position of the figure then use the current point to work
%         % out which subplot it is then get current point in that axis??
%         


%         C = get(flux_plot,'CurrentPoint');
        
        
%         fig;
%         get(fig,'Children');
%         currentPoint
%          C = get(fig,'CurrentPoint')
% pos = get(fig, 'Position'); 

% Now find what object this line relates to in the circuit figure
    end

% 
function NotOver(fig, currentPoint)
    delete(highlight_plot)

end

   function pooopoo(Object,Event)
 C = get(Object,'CurrentPoint');
   end

% Empty object to highlight point where drops or nodes will be placed
grid_hightlight{1} = rectangle('Position',[0,0,0,0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Highlights the point on the grid where the drop or node will be placed,
% at the moment the grid is divide by mm
function GridPointIndicator(~,~)
    
    % Delete the current gird highlighter  
    delete(grid_hightlight{1});
        
    % Get the location of the mouse
    C = get(circuit_plot,'CurrentPoint');
    x = C(1,1); y = C(1,2);
    
    % Find nearest grid point, they are in units of 1 so this is easy 
    x_grid = round(x); y_grid = round(y);
    
    % Plot something on this point using a marker with pickable parts set
    % to none so that it does not respond to mouse clicks
    grid_hightlight{1} = plot(x_grid,y_grid,'d', 'LineWidth',2,...
                         'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],...
                         'MarkerSize',10,'PickableParts', 'none');    
                     
end
% For your information the possible markers are
% '+','o','*','.','x','s','d','^','v','>','<','p','h' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the length of a channel given the connecting objects and the 
% channel width and end points
function L = ChannelLength(xi,yi,connections,W)
    
    % Length of initial channel
    L_temp = sqrt( (xi(2)-xi(1))^2 + (yi(2)-yi(1))^2 );
    
    % Find the extra distance the channel goes into the drop and subtract
    % it, for this the drop radius is needed. If it is a node we do nothing
    if connections{1}(1) == 'd'
        R_i = drop_data{str2double(connections{1}(2:end))}{2}(1);
        alpha_1 = sqrt( R_i^2 + (W/2)^2 );
    else
        alpha_1 = 0;
    end
    if connections{2}(1) == 'd'
        R_i = drop_data{str2double(connections{2}(2:end))}{2}(1);
        alpha_2 = sqrt( R_i^2 + (W/2)^2 );
    else
        alpha_2 = 0;
    end

    % Subtract this to find the length along the edge of the channel
    L = L_temp - alpha_1 - alpha_2;
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

%% Corners
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the corners of a channel given the origins of the two connected
% drops
function [x, y] = corners2(xi,yi,connections,W)
    
    % Length of initial channel
    L_temp = sqrt( (xi(2)-xi(1))^2 + (yi(2)-yi(1))^2 );
    
    % Find the extra distance the channel goes into the drop and subtract
    % it, for this the drop radius is needed. If it is a node we do nothing
    if connections{1}(1) == 'd'
        R_i = drop_data{str2double(connections{1}(2:end))}{2}(1);
        alpha_1 = sqrt( R_i^2 + (W/2)^2 );
    else
        alpha_1 = 0;
    end
    if connections{2}(1) == 'd'
        R_i = drop_data{str2double(connections{2}(2:end))}{2}(1);
        alpha_2 = sqrt( R_i^2 + (W/2)^2 );
    else
        alpha_2 = 0;
    end

    % Subtract this to find the length along the edge of the channel
    L = L_temp - alpha_1 - alpha_2;
    
    % Angle channel makes with horizontal
    phi = atan2(yi(2)-yi(1),xi(2)-xi(1));
    
    x = [ xi(1) - W/2*sin(pi-phi) ...
          xi(1) + W/2*sin(pi-phi) ...
          xi(2) + W/2*sin(pi-phi) ...
          xi(2) - W/2*sin(pi-phi) ];
      
    y = [ yi(1) - W/2*cos(pi-phi) ...
          yi(1) + W/2*cos(pi-phi) ...
          yi(2) + W/2*cos(pi-phi) ...
          yi(2) - W/2*cos(pi-phi) ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  





%% Corners
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the corners of a channel given the line along the middle and the
% width
function [x, y] = Corners(xi,yi,W)
    
    % Angle channel makes with horizontal
    phi = atan2(yi(2)-yi(1),xi(2)-xi(1));
    
    x = [ xi(1) - W/2*sin(pi-phi) ...
          xi(1) + W/2*sin(pi-phi) ...
          xi(2) + W/2*sin(pi-phi) ...
          xi(2) - W/2*sin(pi-phi) ];
      
    y = [ yi(1) - W/2*cos(pi-phi) ...
          yi(1) + W/2*cos(pi-phi) ...
          yi(2) + W/2*cos(pi-phi) ...
          yi(2) - W/2*cos(pi-phi) ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  


%% Arrow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The arrow is a polygon with 5 corners 
function H = Arrow(xi,yi,width,side,length_scale)    
    
    % length of the line, this is used to scale the arrow 
    L = sqrt((xi(1)-xi(2))^2+(yi(1)-yi(2))^2);

    % Distance the arrow is plotted from the line depends on the channel
    % width. The length of the arrow is half the length of the line
    dist = 1+width; length = L/length_scale;
    
    % Angle channel makes with horizontal
    phi = atan2(yi(2)-yi(1),xi(2)-xi(1));
    
    % Coordinates of centre line of the arrow
    xa = [ xi(1) - side*dist*sin(pi-phi) - (L-length)/2*cos(pi-phi) ...
          xi(2) - side*dist*sin(pi-phi) + (L-length)/2*cos(pi-phi)];         
    ya = [ yi(1) - side*dist*cos(pi-phi) + (L-length)/2*sin(pi-phi)...
          yi(2) - side*dist*cos(pi-phi) - (L-length)/2*sin(pi-phi)];          
    
    %   
    x = [ xa(1) - W/2*sin(pi-phi) ...
          xa(1) + W/2*sin(pi-phi) ...
          xa(2) + W/2*sin(pi-phi) ...
          xa(2) - W/2*sin(pi-phi) ];
      
    y = [ ya(1) - W/2*cos(pi-phi) ...
          ya(1) + W/2*cos(pi-phi) ...
          ya(2) + W/2*cos(pi-phi) ...
          ya(2) - W/2*cos(pi-phi) ];
      
  
      
    H = fill(xd,yd,'k','edgecolor','k');

      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Functions attached to the buttons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% They mostly just switch what happens when you click the 
% corresponding button
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Add_drop(~, ~)
    Add_drop_on = true;
    Add_node_on = false;
    Add_channel_on = false;
    select_on = false;
    channel_clicks = 1;
end
function Add_node(~, ~)
    Add_drop_on = false;
    Add_node_on = true;
    Add_channel_on = false;
    select_on = false;
    channel_clicks = 1;
end
function Add_channel(~, ~)
    channel_connections = [];
    Add_drop_on = false;
    Add_node_on = false;
    Add_channel_on = true;
    select_on = false;
    channel_clicks = 1;
end
function Select(~, ~)
    Add_drop_on = false;
    Add_node_on = false;
    Add_channel_on = false; 
    select_on = true;
    channel_clicks = 1;
    % Make the current object empty
    current_object = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  


%% Functions for deleting things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The objects are labelled An, where A = {c,d,n} for channel, drop or
% node and n is a number. To delete An we need to remove An from all 
% arrays in which it occurs, then for all i>n we redefine An = A(n-1).
% Finally we decrease the relevant count for this object by 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deletes a channel, removes it from all arrays and relabels everything 
% accordingly
function DeleteChannel(Object)
    
    % Determine what it's position in the data array is 
    pos = str2double(Object.Tag(2:end));

    % This label needs to be removed from all connected nodes and drops
    % and they need to have the other labels updated
    for i = 1:drop_count

        % Find all connected channels for this drop
        connections = drop_data{i}{3};

        % Remove the entry for the deleted channel
        connections(connections == pos) = [];

        % The other channels are relabelled by subtracting 1 from their
        % number 
        connections(connections > pos) = connections(connections > pos) - 1;

        % Put the new connections into the data array
        drop_data{i}{3} = connections;

    end

    for i = 1:node_count

        % Find all connected channels for this drop
        connections = node_data{i}{3};

        % Remove the entry for the deleted channel
        connections(connections == pos) = [];

        % The other channels are relabelled by subtracting 1 from their
        % number 
        connections(connections > pos) = connections(connections > pos) - 1;

        % Put the new connections into the data array
        node_data{i}{3} = connections;

    end

    % Update the labels for each entry after the one selected
    for i = pos+1:channel_count
        channel_data{i}{1}.Tag =  ['c' num2str(i-1)];
    end

    % Delete the arrow associated with this channel    
    delete(channel_arrow{pos});
    channel_arrow{pos} = [];
    
    % Delete the graphics object
    delete(Object)

    % Remove this data from the array
    channel_data(pos) = [];

    % Decrease the channel count
    channel_count = channel_count - 1;    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deletes a drop or node
function DeleteHubs(Object)

    % Find what we are dealing with
    type = Object.Tag(1);
    pos = str2double(Object.Tag(2:end));
    
    % If it's a drop
    if type == 'd'
        
        % First find all the channels conncted to this drop and delete
        % those by calling the DeleteChannel function. The channel numbers
        % are updated every time one is deleted, so just take the first
        % entry in this array each time it is updated until it is empty
        things_to_delete = length(drop_data{pos}{3});
        for i = 1:things_to_delete
            DeleteChannel(channel_data{drop_data{pos}{3}(1)}{1});
        end
        
        % The labels for this drop only appeared in the deleted channels, 
        % so we only need renumber the labels for each entry after the one 
        % selected
        for i = pos+1:drop_count
            drop_data{i}{1}.Tag =  ['d' num2str(i-1)];
        end

        % Delete the graphics object
        delete(Object)

        % We also need to delete any associated labels 
        delete(Vol_label{pos}); 
        
        % Do the same again for the extra labels if they exist
        if extra_labels
            delete(pressure_label{pos});
            delete(radius_label{pos});
        end
            
        % Remove this data from the array
        drop_data(pos) = [];

        % Decrease the drop count
        drop_count = drop_count - 1;
        
    % If it's a node we do pretty much the same thing
    elseif type == 'n'
        
        % First find all the channels conncted to this node and delete
        % those by calling the DeleteChannel function 
        things_to_delete = length(node_data{pos}{3});
        for i = 1:things_to_delete
            DeleteChannel(channel_data{node_data{pos}{3}(1)}{1});
        end
        
        % The labels for this node only appeared in the deleted channels, 
        % so we only need renumber the labels for each entry after the one 
        % selected
        for i = pos+1:node_count
            node_data{i}{1}.Tag =  ['n' num2str(i-1)];
        end

        % Delete the graphics object
        delete(Object)

        % Remove this data from the array
        node_data(pos) = [];

        % Decrease the drop count
        node_count = node_count - 1;
        
    end
    
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% ClearCircuit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removes graphics objects from the axis, clears the data arrays and sets 
% the counters to zero 
function ClearCircuit(~, ~)

    % Delete the objects from the axis first 
    cla(circuit_plot);
    
    % Add a line indicating distance between tick marks
    x_indicator = [1 11]; y_indicator = [47 47];
    line(x_indicator,y_indicator,'LineWidth',1,'color','k');
    line([x_indicator(1) x_indicator(1)],[y_indicator(1)-0.4 y_indicator(2)+0.4],'LineWidth',1,'color','k');
    line([x_indicator(2) x_indicator(2)],[y_indicator(1)-0.4 y_indicator(2)+0.4],'LineWidth',1,'color','k');

    % Label the line
    text(mean(x_indicator),mean(y_indicator)+1,'10 mm','HorizontalAlignment','center','Interpreter','Latex','FontSize',fontSize);

    % Return the variables to their original state
    drop_count = 0;
    drop_data = cell(0);
    node_count = 0;
    node_data = cell(0);
    channel_count = 0;
    channel_data = cell(0);
    channel_arrow = cell(0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 

%% UpdateDrop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changes the volume and base radius of the drop specified 
function UpdateDrop(pos,radius,volume)
        
    % Coordinates for the new drop
    xd = drop_radius*sin(theta) + drop_data{pos}{2}(2);
    yd = drop_radius*cos(theta) + drop_data{pos}{2}(3);

    % Update the vertices
    current_object.Vertices = [xd', yd'];
            
    % The drop base radius may have changed so we need to update the
    % pressure function
    P_fcn{pos} = @(V) pressure_fcn(V,1e-3*radius,model);
    
    % Make sure this is also in the data array
    drop_data{pos}{2}(1) = radius;
    drop_data{pos}{2}(4) = volume;
    drop_data{pos}{2}(5) = P_fcn{pos}(1e-9*volume);
            
    % Delete the old labels and replace with updated ones
    delete(Vol_label{pos});
    Vol_label_string = [num2str(volume) '\fontname{Times}' char(181) 'l']; 
    Vol_label{pos} = text(drop_data{pos}{2}(2),drop_data{pos}{2}(3),Vol_label_string,'HorizontalAlignment','center',...
                     'Interpreter','tex','FontSize',fontSize,'PickableParts','none');    
    
    % Do the same again fr the extra ones if they exist
    if extra_labels
        delete(pressure_label{pos});
        delete(radius_label{pos});
        pressure_label_string = [num2str(drop_data{pos}{2}(5)) '\fontname{Times} Pa']; 
        pressure_label{pos} = text(drop_data{pos}{2}(2),drop_data{pos}{2}(3)-2.5,pressure_label_string,'HorizontalAlignment','center',...
                         'Interpreter','tex','FontSize',fontSize,'PickableParts','none');                    
        radius_label_string = [num2str(drop_radius) '\fontname{Times} mm']; 
        radius_label{pos} = text(drop_data{pos}{2}(2),drop_data{pos}{2}(3)+2.5,radius_label_string,'HorizontalAlignment','center',...
                         'Interpreter','tex','FontSize',fontSize,'PickableParts','none');                                                                                            
    end
    
    % Test if this now overlaps any existing drops, channels or nodes and
    % print a warning somewhere if so 
    if OverlapTest()
        fprintf('Print this warning somewhere \n')
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updates the width of a channel 
function UpdateChannel(object,width)
    
    % Check that this is a channel, if not we do nothing
    if object.Tag(1) == 'c'
    
        % Find position in channel_data array
        pos1 = str2double(object.Tag(2:end));

        % The end points are unchanged, but the corners need to be
        % recalculated
        ends = channel_data{pos1}{4};
        connections = channel_data{pos1}{3};

        % The corners of the channel are given by              
        [xc, yc] = Corners(ends(:,1), ends(:,2), width);

        % Find the length of the channel along the outer edge
        length = ChannelLength(ends(:,1), ends(:,2), connections, width);

        % Update the vertices
        object.Vertices = [xc', yc'];

        % Make sure this is also in the data array
        channel_data{pos1}{2}(1) = width;
        channel_data{pos1}{2}(2) = length; 

        % Delete the arrow associated with this channel    
        delete(channel_arrow{pos1});
        channel_arrow{pos1} = [];

        % Update width and length labels if they exist
        if extra_labels
            
            % Delete the old labels
            delete(width_label{pos1});
            delete(length_label{pos1});
            
            % Angle channel makes with horizontal
            phi = atan2(ends(2,2)-ends(1,2),ends(2,1)-ends(1,1));
             
            % Mid-point of channel
            x_mid = ends(1,1) + length/2*cos(phi);
            y_mid = ends(1,2) + length/2*sin(phi);
        
            % The labels
            width_label_string = [num2str(width) '\fontname{Times} mm']; 
            width_label{pos1} = text(x_mid,y_mid+1,width_label_string,'HorizontalAlignment','center',...
                                    'Interpreter','tex','FontSize',fontSize,'PickableParts','none');      
            length_label_string = [num2str(length) '\fontname{Times} mm']; 
            length_label{pos1} = text(x_mid,y_mid-1,length_label_string,'HorizontalAlignment','center',...
                                    'Interpreter','tex','FontSize',fontSize,'PickableParts','none');         
        
        end
 
        % Work out which direction the initial flow is in. It goes from
        % high pressure to low, so for the two connected things find
        % the highest pressure
        P = zeros(1,2);
        for i = 1:2
            type = connections{i}(1);
            pos2 = str2double(connections{i}(2:end));
            if type == 'd'
                P(i) = drop_data{pos2}{2}(5);
            elseif type == 'n'
                P(i) = node_data{pos2}{2}(3); 
            end        
        end

        if P(1) > P(2)  
            channel_arrow{pos1} = DrawArrow(ends(:,1), ends(:,2), width, 1, 4);
        elseif P(1) < P(2)
            channel_arrow{pos1} = DrawArrow(flipud(ends(:,1)), flipud(ends(:,2)), width, 1, 4);
        else
            channel_arrow{pos1} = [];
        end

        % Test if this now overlaps any existing drops, channels or nodes and
        % print a warning somewhere if so 
        if OverlapTest()
            fprintf('Print this warning somewhere \n')
        end    

    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Updates the width of a channel using global variables (which are
% % altered by scrolling) 
% function UpdateChannel(object)
%     
%     % Check that this is a channel, if not we do nothing
%     if object.Tag(1) == 'c'
%     
%         % Find position in channel_data array
%         pos1 = str2double(object.Tag(2:end));
% 
%         % The end points are unchanged, but the corners need to be
%         % recalculated
%         channel_ends = channel_data{pos1}{4};
%         channel_connections = channel_data{pos1}{3};
% 
%         % The corners of the channel are given by              
%         [xc, yc] = corners(channel_ends(:,1), channel_ends(:,2), channel_width);
% 
%         % Find the length of the channel along the outer edge
%         channel_length = channelLength(channel_ends(:,1), channel_ends(:,2), channel_connections, channel_width);
% 
%         % Update the vertices
%         object.Vertices = [xc', yc'];
% 
%         % Make sure this is also in the data array
%         channel_data{pos1}{2}(1) = channel_width;
%         channel_data{pos1}{2}(2) = channel_length; 
% 
%         % Delete the arrow associated with this channel    
%         delete(channel_arrow{pos1});
%         channel_arrow{pos1} = [];
% 
%         % Work out which direction the initial flow is in. It goes from
%         % high pressure to low, so for the two connected things find
%         % the highest pressure
%         P = zeros(1,2);
%         for i = 1:2
%             type = channel_connections{i}(1);
%             pos2 = str2double(channel_connections{i}(2:end));
%             if type == 'd'
%                 P(i) = drop_data{pos2}{2}(5);
%             elseif type == 'n'
%                 P(i) = node_data{pos2}{2}(3); 
%             end        
%         end
% 
%         if P(1) > P(2)  
%             channel_arrow{pos1} = DrawArrow(channel_ends(:,1), channel_ends(:,2), channel_width, 1, 4);
%         elseif P(1) < P(2)
%             channel_arrow{pos1} = DrawArrow(flipud(channel_ends(:,1)), flipud(channel_ends(:,2)), channel_width, 1, 4);
%         else
%             channel_arrow{pos1} = [];
%         end
% 
% 
%     end
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% NodePressures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the pressure in each node. The values for the pressure in a node is
% given by an algebraic equation. The node may be connected to another node
% (with unknown pressure) so we will to solve a linear system:
% (I - M)PJ^4 = [sum(Ai*betai^4Vi^4)/sum(Ai)] for i in connections, 
% where Ai is ai^7/(105*gamma^3*mu1*Li) and M(i,k) is Ai/sum(Ai) if node i
% is connected to node k and zero otherwise 
function NodePressures()
     
    % Empty vector and matrix 
    sys_RHS = zeros(1,node_count);  
    node_connections_matrix = zeros(node_count);  

    for i0_0 = 1:node_count

        % We need to sum some values 
        pressure_sum = 0;
        flux_const_sum = 0; 
 
        % Find all the channels connected to it and find the flux in them
        for i0_1 = node_data{i0_0}{3}
            
            % Length and half width of the channel 
            a = 1e-3*channel_data{i0_1}{2}(1)/2;
            L = 1e-3*channel_data{i0_1}{2}(2);
            
            % The flux is multiplied by a constant in the thin-film models
            sigma = a^7/(105*gamma^3*mu1*L);             
            
            % Find if the other end is a drop or node, so first find which
            % end we are dealing with
            % If the first entry in the channel connections is this
            % drop we look at the second entry, otherwise we want the
            % first
            str = ['n' num2str(i0_0)];
            if strcmp(channel_data{i0_1}{3}{1},str)                 
                type = channel_data{i0_1}{3}{2}(1);
                pos = str2double(channel_data{i0_1}{3}{2}(2:end));
            else                
                type = channel_data{i0_1}{3}{1}(1);
                pos = str2double(channel_data{i0_1}{3}{1}(2:end));               
            end
                         
            if type == 'd'
                % If it's connected to a drop we use the pressure function 
                flux_const_sum = flux_const_sum - sigma; 
                pressure_sum = pressure_sum - sigma * P_fcn{pos}(1e-9*drop_data{pos}{2}(4))^4;
            else
                % Otherwise it's connected to a node and the pressure is an
                % unknown  
                flux_const_sum = flux_const_sum - sigma; 
                node_connections_matrix(i0_0,pos) = -sigma;                   
            end
                        
        end
                    
        node_connections_matrix(i0_0,:) = node_connections_matrix(i0_0,:)/flux_const_sum;
        sys_RHS(i0_0) = pressure_sum/flux_const_sum;          
        
    end
                 
    % Initial values for the node pressures 
    node_pressures = ((diag(ones(node_count,1))-node_connections_matrix)\sys_RHS').^(1/4);
            
    % There may be NaNs in this solution, which correspond to nodes that
    % are not connected to anything, the pressure i these is set to zero
    node_pressures(isnan(node_pressures)) = 0;
    
    for i0_0 = 1:node_count       
       % Add the pressures to node data array
       node_data{i0_0}{2}(3) = node_pressures(i0_0);  
       
       % Update the labels for each node if they exist
       if extra_labels
           delete(node_label{i0_0});       
           node_label_string = [num2str(node_pressures(i0_0)) '\fontname{Times} Pa']; 
           node_label{i0_0} = text(node_data{i0_0}{2}(1),node_data{i0_0}{2}(2),node_label_string,'HorizontalAlignment','center',...
                                    'Interpreter','tex','FontSize',fontSize,'PickableParts','none'); 
       end
    end        

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for drawing arrows using the matlab fill function, 7 points need
% to be defined to make the arrow polygon. Input determines which side of 
% the line the arrow is plotted on, it takes values -1 or 1 
function H = DrawArrow(xi,yi,width,side,length_scale)
  
    % Width of the body of the arrow and of the head
    W_body = 0.2; W_head = 0.6;

    % length of the line, this is used to scale the arrow 
    L = sqrt((xi(1)-xi(2))^2+(yi(1)-yi(2))^2);

    % Distance the arrow is plotted from the line depends on the channel
    % width. The length of the arrow is half the length of the line
    dist = 1+width; length = L/length_scale;
    
    % Angle channel makes with horizontal
    phi = atan2(yi(2)-yi(1),xi(2)-xi(1));
    
    % Coordinates of centre line of the arrow
    xa = [ xi(1) - side*dist*sin(pi-phi) - (L-length)/2*cos(pi-phi) ...
           xi(2) - side*dist*sin(pi-phi) + (L-length)/2*cos(pi-phi)];         
    ya = [ yi(1) - side*dist*cos(pi-phi) + (L-length)/2*sin(pi-phi)...
            yi(2) - side*dist*cos(pi-phi) - (L-length)/2*sin(pi-phi)];          
    
        
    % Coordinates for arrow corners       
    x = [ xa(2) + cos(phi)...         
          xa(2) + W_head*sin(pi-phi) ...
          xa(2) + W_body/2*sin(pi-phi) ... 
          xa(1) + W_body/2*sin(pi-phi) ...
          xa(1) - W_body/2*sin(pi-phi) ...
          xa(2) - W_body/2*sin(pi-phi) ...
          xa(2) - W_head*sin(pi-phi) ];      
    y = [ ya(2) + sin(phi) ...      
          ya(2) + W_head*cos(pi-phi) ...
          ya(2) + W_body/2*cos(pi-phi) ...
          ya(1) + W_body/2*cos(pi-phi) ...
          ya(1) - W_body/2*cos(pi-phi) ...
          ya(2) - W_body/2*cos(pi-phi) ...
          ya(2) - W_head*cos(pi-phi) ];
      
    % Define the arrow patch object  
    H = fill(x,y,'k','edgecolor','k');
            
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% flux_fcn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for defining the flux constant and the pressure as a function 
% of volume for use in the ODEs. Different models can be used for each
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q = flux_fcn(channel_width,channel_length,drop_R,connected_drop_R,drop_V,node_P,connected_drop_V,connected_node_P,model)
          
    if model == 1 || 2 
        % The flux is multiplied by a constant in the thin-film models
        sigma = channel_width^7/(105*gamma^3*mu1*channel_length);  
        % Supply either the unknown pressure for a node or the volume for a
        % drop, DON'T GIVE BOTH. the there are four possibilities
        if ~isempty(drop_V) && ~isempty(connected_drop_V)
            Q = -sigma * (pressure_fcn(drop_V,drop_R,model)^4 - pressure_fcn(connected_drop_V,connected_drop_R,model)^4); 
        end
        if ~isempty(drop_V) && isempty(connected_drop_V)
            Q = -sigma * (pressure_fcn(drop_V,drop_R,model)^4 - connected_node_P^4); 
        end
        if isempty(drop_V) && ~isempty(connected_drop_V)
            Q = -sigma * (node_P^4 - pressure_fcn(connected_drop_V,connected_drop_R,model)^4); 
        end
        if  isempty(drop_V) && isempty(connected_drop_V)
            Q = -sigma * (node_P^4 - connected_node_P^4); 
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% pressure_fcn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure as a function of volume for each model
function P = pressure_fcn(volume,base_radius,model)
    % Pressure is linear function of volume for the thin-film models 
    if model == 1
        
        beta = B*gamma*besselj(0,sqrt(B)*base_radius) / (pi*base_radius^2*besselj(2,sqrt(B)*base_radius));
        P = beta*volume;
        
    elseif model == 2 
        
    % Load the data, this was scaled using the drop base radius
%     thickPressure = load('Pressures');

    % Contour3 also plots the contour (the function contourc gets just the 
    % data but the input needs to be vectors and doesn't work with the
    % multivalued stuff I have here. So I'll create a figure and delete it,
    % this is STUPID!
     
    stupid = figure; hold on
        
    % We need to extract to curves for the pressure in each drop, we do
    % this by finding the abropriate contour in each drop
    M1 = contour3(thickPressure.V,thickPressure.p,thickPressure.Bo,[B*base_radius^2 B*base_radius^2]);
    
    delete(stupid)
    
    % Extract the data from the contours, the point at the origin seems to
    % be missing so we add this on
    V_full = [0 M1(1,2:end)]; P_full = [0 M1(2,2:end)];

    % Also the output may contain non-unique values, which the interpolator
    % doesn't like. These are removed from the volume vector and the
    % corressponding pressure vector, while we're at it we also
    % dimensionalise the data
    for i0 = 1:2
        [~, ind] = unique(V_full);
        V_full = base_radius^3*V_full(ind);
        P_full = gamma/base_radius*P_full(ind);
    end   

    % Interpolation is then used to find the pressure for any volume (in
    % range) 
    P = interp1(V_full,P_full,volume,'spline');
    
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Callback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defines the clicked graphics object as the current object and maybe does
% stuff to it
function Callback(object,~)
    
    % Make clicked object the current object
    current_object = object; 
    
    % Put variables for clicked object into global variables
    if current_object.Tag(1) == 'd'
        drop_radius = drop_data{str2double(current_object.Tag(2:end))}{2}(1);
        drop_volume = drop_data{str2double(current_object.Tag(2:end))}{2}(4);    
    elseif current_object.Tag(1) == 'c'
        channel_width = channel_data{str2double(current_object.Tag(2:end))}{2}(1);
    end
    
    % If plot_channel is true we get the data from the clicked object and 
    % store it in a global variable. When two (different) things have been 
    % clicked they are connected by a channel. We also require the clicked
    % thing to be a node or a drop.
    if Add_channel_on == true && (current_object.Tag(1) == 'd' || current_object.Tag(1) == 'n')
        
        % Two points are needed to specify the channel so we count how many
        % clicks there have been
        if channel_clicks == 1 % first point
            
            % Clear the vector containing the ends of the channel 
            channel_ends = [];
            
            % Get the tag of the current object and add it to the
            % channel_connections variable
            channel_connections = object.Tag;
                        
            % Get ready for next click
            channel_clicks = 2;
            
         elseif channel_clicks == 2 % second point
        
            % Get the tag of the current object and add it to the 
            % channel_connections variable
            channel_connections = {channel_connections, object.Tag};

            % Check that there does not already exist a channel
            % connecting these things
            if channel_count > 0
                for i0_0 = 1:channel_count 
                    if (strcmp(channel_connections{1},channel_data{i0_0}{3}{1}) || ...
                        strcmp(channel_connections{1},channel_data{i0_0}{3}{2})) && ...
                       (strcmp(channel_connections{2},channel_data{i0_0}{3}{1}) || ...
                        strcmp(channel_connections{2},channel_data{i0_0}{3}{2}))
                        % There is a channel connecting these things
                        % already and we reset the click counter
                        channel_clicks = 1;   
                    end
                end
            end  
                        
            % We only do something if it's not the same as the previous 
            % thing clicked and if such a channel does not already exist
            if ~strcmp(channel_connections{1},channel_connections{2}) && channel_clicks == 2
                
                % Increase the channel counter by 1
                channel_count = channel_count + 1;

                % Get the coordinates of the centres of the connected objects
                % at each end of the channel 
                for i1_0 = 1:2
                    type = channel_connections{i1_0}(1);
                    pos = str2double(channel_connections{i1_0}(2:end));
                   if type == 'd'

                        % Get the centre of this object (where this data is depends
                        % on the object type)
                        channel_ends = [channel_ends; drop_data{pos}{2}(2:3)];

                        % Add this connection to the drop data array
                        drop_data{pos}{3} = [drop_data{pos}{3} channel_count];

                    elseif type == 'n'

                        % Get the centre of the node
                        channel_ends = [channel_ends; node_data{pos}{2}(1:2)];

                        % Add this connection to the node data array
                        node_data{pos}{3} = [node_data{pos}{3} channel_count];

                   end
                end

                % The initial channel width is 
                width = 1;
                
                % Draw the channel 
                PlaceChannels(channel_ends(:,1),channel_ends(:,2),width,channel_connections)
                
                % This channel will change the initial pressure values if
                % it is connected to a node so we update the nodes if they 
                % exist
                if node_count > 0
                    NodePressures()
                end
                
                % This also means that the fluxes will need to be updated,
                % it probably isn't neccessary to update every channel
                for i2_0 = 1:channel_count  
                    UpdateChannel(channel_data{i2_0}{1},channel_data{i2_0}{2}(1)) 
                end
                
            end
            
        end
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      


%% OverlapTest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function tests if a new or existing object overlaps an existing drop 
% channel or node. A pre-existing drop or node may be connected to a 
% channel, in which case they will overlap so any 'legitametly' connected 
% things need to be ignored
function overlap = OverlapTest()
    
    % Matrix of things that should be connected and empty vector for the
    % polyshape objects
    connections = eye(drop_count+node_count+channel_count);
    polyvec = polyshape.empty(drop_count+node_count+channel_count,0);
    
    
    % I will get the vertices of each object, input them into the 
    % polyshape function and then find if these polygons overlap, this 
    % seems a bit convoluted! 
    for i0_0 = 1:drop_count
        polyvec(i0_0) = simplify(polyshape(drop_data{i0_0}{1}.Vertices));%,'Simplify',true);        
        % Find which channels, if any, this is connected to and add a 1 in
        % the connections matrix (it will be symmetric) 
        for i0_1 = drop_data{i0_0}{3}
            connections(i0_0,drop_count+node_count+i0_1) = 1; connections(drop_count+node_count+i0_1,i0_0) = 1;
            % All channels meeting in this drop should also be connected
            for i0_2 = drop_data{i0_0}{3}
                connections(drop_count+node_count+i0_2,drop_count+node_count+i0_1) = 1; connections(drop_count+node_count+i0_1,drop_count+node_count+i0_2) = 1;
            end
        end
    end
    for i1_0 = 1:node_count
        polyvec(drop_count+i1_0) = simplify(polyshape(node_data{i1_0}{1}.Vertices));%,'Simplify',true);
        % Find which channels, if any, this is connected to and add a 1 in
        % the connections matrix (it will be symmetric) 
        for i1_1 = node_data{i1_0}{3}
            connections(drop_count+i1_0,drop_count+node_count+i1_1) = 1; connections(drop_count+node_count+i1_1,drop_count+i1_0) = 1;
            % All channels meeting in this node should also be connected
            for i1_2 = node_data{i1_0}{3}
                connections(drop_count+node_count+i1_2,drop_count+node_count+i1_1) = 1; connections(drop_count+node_count+i1_1,drop_count+node_count+i1_2) = 1;
            end
        end
    end
    for i2_0 = 1:channel_count
        polyvec(drop_count+node_count+i2_0) = simplify(polyshape(channel_data{i2_0}{1}.Vertices));%,'Simplify',true);
    end

    % The matlab function overlaps then gives a matrix where m(i,j) is
    % 1 if objects i and j overlap and 0 otherwise
    TF = overlaps(polyvec);
        
    % We subtract all the known connections from this, if there are any
    % non-zero entries there is an overlap
    overlap = sum(sum(TF - connections)) > 1;
       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% PlaceDropsAndNodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Plots drops or nodes when circuit_plot is clicked on by finding the x-y
% coordinate of mouse click and using this to position the appropriate
% object           
function PlaceDropsAndNodes(Object,~)

    % Object is the axis, current point gives extra dimensions, extract
    % only x and y and round them to the nearest grid point
    C = get(Object,'CurrentPoint');
    x = round(C(1,1)); y = round(C(1,2));
    
    % The booleans are used to determine if anything should be done when 
    % the axis is clicked    
    % Add a drop
    if Add_drop_on == true
        
        % Initial radius in millimetres and volume in microlitres
        drop_radius = 2.5;
        drop_volume = 10;
        
        % The coordinates for the drop origin are x and y, so the edge is 
        xd = drop_radius*sin(theta) + x;
        yd = drop_radius*cos(theta) + y;        
        
        % Plot the new drop and check that it does not overlab with any
        % existing objects
        new_drop = fill(xd,yd,fill_colour,'edgecolor',edge_colour,'ButtonDownFcn',@Callback);

        % Increase the drop count by one and plot it with the a function 
        % handle allowing a function to be run when the object is clicked
        drop_count = drop_count + 1;
        
        % Use the volume (in metres^3) and radius (in metres) to find the 
        % pressure function for this drop (in Pascals)
        P_fcn{drop_count} = @(V) pressure_fcn(V,1e-3*drop_radius,model);

        % Find the initial pressure
        initial_pressure = P_fcn{drop_count}(1e-9*drop_volume);
               
        % Non-clickable text boxes with the base radius, Laplace pressure
        % and volume, THIS IS MORE FOR TEST PURPOSES
        Vol_label_string = [num2str(drop_volume) '\fontname{Times}' char(181) 'l']; 
        Vol_label{drop_count} = text(x,y,Vol_label_string,'HorizontalAlignment','center',...
                                'Interpreter','tex','FontSize',fontSize,'PickableParts','none'); 
        % Only add the radius and pressure ones if wanted
        if extra_labels                            
            pressure_label_string = [num2str(initial_pressure) '\fontname{Times} Pa']; 
            pressure_label{drop_count} = text(x,y-2.5,pressure_label_string,'HorizontalAlignment','center',...
                         'Interpreter','tex','FontSize',fontSize,'PickableParts','none');
            radius_label_string = [num2str(drop_radius) '\fontname{Times} mm']; 
            radius_label{drop_count} = text(x,y+2.5,radius_label_string,'HorizontalAlignment','center',...
                         'Interpreter','tex','FontSize',fontSize,'PickableParts','none');                                                          
        end
        
        % Add a tag so it can be identified and then add the new drop to
        % the cell array with the geometric parameters and the labels of
        % the connected channels (initially there are none)
        connections = [];
        geometry = [drop_radius x y drop_volume initial_pressure];
        new_drop.Tag = ['d' num2str(drop_count)];
        drop_data{drop_count} = {new_drop geometry connections};

        % Make the new drop the current object so it's volume and radius
        % can be modified by scrolling
        current_object = new_drop;
        
        % Test if this overlaps any existing drops, channels or nodes and
        % print a warning if so            
        if OverlapTest()
            fprintf('warning\n')
        end                
        
    end
    
    % Add a node
    if Add_node_on == true
        
        % This is very similar to adding a drop, could be done differently
        xn = sin(theta) + x;
        yn = cos(theta) + y;        
                
        % Plot the new node and check that it does not overlab with any
        % existing objects
        new_node = fill(xn,yn,fill_colour,'edgecolor',edge_colour,'ButtonDownFcn',@Callback);
                        
        % When the node is not connected to anything the pressure is zero, 
        % this will be updated when channels are connected
        node_pressure = 0;

        % Increase the node count by one and plot it with the function
        % handle so stuff can be done when it is clicked, e.g. adding a
        % channel
        node_count = node_count + 1;        
                
        % Non-clickable text boxes with the Laplace pressure 
        if extra_labels
            node_label_string = [num2str(node_pressure) '\fontname{Times} Pa']; 
            node_label{node_count} = text(x,y,node_label_string,'HorizontalAlignment','center',...
                                    'Interpreter','tex','FontSize',fontSize,'PickableParts','none');      
        end
        
        % Add a tag so it can be identified and then add the new node to
        % the cell array with the geometric parameters and the labels of
        % the connected channels (initially there are none)
        connections = [];
        geometry = [x y node_pressure];
        new_node.Tag = ['n' num2str(node_count)];
        node_data{node_count} = {new_node geometry connections};
        
        % Make the new node the current object 
        current_object = new_node;

        % Test for overlaps             
        if OverlapTest()
            fprintf('warning\n')
        end
                 
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% PlaceChannels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Plots a channel between to two specified objects and defines the channel
% object       
function PlaceChannels(x,y,width,connections)                

    % The corners of the channel are given by              
    [xc, yc] = Corners(x, y, width);
                
    % The length along the outer edge is 
    length = ChannelLength(x, y, connections, width);
                                                
    % Plot the new channel with callback function attached
    new_channel = fill(xc,yc,fill_colour,'edgecolor',edge_colour,'ButtonDownFcn',@Callback);

    % Plot the channels on the bottom layer so the drops and nodes
    % can be clicked on more easily
    uistack(new_channel,'bottom');    
        
    % Add a tag so it can be identified and then add the new 
    % channel to the cell array with the geometric parameters and 
    % the labels of the things it connects 
    geometry = [width length];
    new_channel.Tag = ['c' num2str(channel_count)];
    channel_data{channel_count} = {new_channel geometry connections [x y]};

    % Make the new channel the current object so it's width 
    % can be modified by scrolling
    current_object = new_channel;

    % Non-clickable text boxes with the width and length of the channel
    if extra_labels
        
        % Angle channel makes with horizontal
        phi = atan2(y(2)-y(1),x(2)-x(1));
        
        % Mid-point of channel
        x_mid = x(1) + length/2*cos(phi);
        y_mid = y(1) + length/2*sin(phi);

        % The labels
        width_label_string = [num2str(width) '\fontname{Times} mm']; 
        width_label{channel_count} = text(x_mid,y_mid+1,width_label_string,'HorizontalAlignment','center',...
                                'Interpreter','tex','FontSize',fontSize,'PickableParts','none');      
        length_label_string = [num2str(length) '\fontname{Times} mm']; 
        length_label{channel_count} = text(x_mid,y_mid-1,length_label_string,'HorizontalAlignment','center',...
                                'Interpreter','tex','FontSize',fontSize,'PickableParts','none');         
    
    end
        
    % Reset the click counter
    channel_clicks = 1;     

    % Work out which direction the initial flow is in. It goes from
    % high pressure to low, so for the two connected things find
    % the highest pressure
    P = zeros(1,2);
    for i0_0 = 1:2
        type = connections{i0_0}(1);
        pos = str2double(connections{i0_0}(2:end));
        if type == 'd'
            P(i0_0) = drop_data{pos}{2}(5);
        elseif type == 'n'
            P(i0_0) = node_data{pos}{2}(3);
        end
    end             

    if P(1) > P(2)  
        channel_arrow{channel_count} = DrawArrow(x, y, width, 1, 4);
    elseif P(1) < P(2)
        channel_arrow{channel_count} = DrawArrow(flipud(x), flipud(y), width, 1, 4);
    else
        channel_arrow{channel_count} = [];
    end
    
    % Test if this now overlaps any existing drops, channels or nodes and
    % print a warning somewhere if so 
    if OverlapTest()
        fprintf('Print this warning somewhere \n')
    end    
        
                                
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





end