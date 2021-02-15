%This code is desgined to simulate system 4.1. Intolerance case.
%By Belen Serrano Anton
%Created 25/02/2020
%Last Modified 31/03/2020

%Variable definition
%t_cycle = 0.15;                         %Time lap between the restriction point and cell division
t_cycle = 0.05; 
%t_apo = 0.2;                           %Time lap between the deactivation of Bcl-2 and cell death
t_apo = 0.1;
%t_next = 0.3;                           %Time step in this simulation
t_next = 0.15;

%Parameters: Pathogen
%alpha = 6;                              %Pathogen proliferation rate
alpha = 6.4;
%beta = 0.04;                            %Pathogen death rate
beta = 0.22;
%Parameters: Effector T cells

lambda_pd = 0.05;                           %Change rate in membrane receptor Rd, due to Rp signals
lambda_taup = 6*10^(-5);                    %Change rate in membrane receptor Rd, due to TCR signals
lambda_pp = 0.5*10^(-4);                    %Change rate in membrane receptor Rp, due to Rp signals                           
mu_pc = 8;                                  %Change rate in inhibitor molecule Rb, due to receptor Rc
mu_da = 15;                                 %Change rate in inhibitor molecule Bcl-2, due to receptor Rc

%Parameters: Memory T cells
lambda_pd_mem = 0;                      %Change in membrane receptor Rd, due to Rp signals
lambda_taup_mem = 10^(-5);              %Change rate in membrane receptor Rd, due to TCR signals
lambda_pp_mem = 2*10^(-2);              %Change rate in membrane receptor Rp, due to Rp signals
mu_pc_mem = 13;                         %Change rate in inhibitor molecule Rb, due to receptor Rc

%Define the final time we will simulate to
T_final = 25;

%Define the initial number of particles
N_init = 25;                            %N will represent T cells                
Y_init = 5;                             %Y will represent pathogen

%Define how long the recording vectors should be
num_rec_steps = round(T_final/t_next);

%Initialise the index which will tell us where to write the current values
rec_ind=1;

%Define the maximum number of t cells
num_max_cells=10^7;

%Instantiate a vector which will hold the time varying values of T cells
%and pathogen
rec_vector_N_eff = -ones(1,num_rec_steps);      %For effector T cells
rec_vector_N_mem = -ones(1,num_rec_steps);      %For memory T cells
rec_vector_Y = -ones(1,num_rec_steps);          %For the pathogen

%Write the initial condition to these vectors
rec_vector_N_eff(rec_ind) = N_init;             %Asymmetric division of naive T cells
rec_vector_N_mem(rec_ind) = N_init;             %Asymmetric division of naive T cells
rec_vector_Y(rec_ind) = Y_init;

%Instantiate a vector which will hold the t cells
t_cell_matrix=zeros(num_max_cells,6);

%Write the initial condition to this vector
a0 = 0.3;
c0 = 0.08;
c0_mem = 0.04;

t_cell_matrix(1:2:2*N_init,1)=1;                %type 1: Effector T cell
t_cell_matrix(1:2:2*N_init,2)=a0;
t_cell_matrix(1:2:2*N_init,3)=c0;

t_cell_matrix(2:2:2*N_init,1)=2;                %type 2: Memory T cell
t_cell_matrix(2:2:2*N_init,3)=c0_mem; 


%Initialise a vector which will hold the times when reactions occur
time_vec=zeros(1,num_rec_steps);

%Initialise the number of particles for this repeat
N_eff = N_init;
N_mem = N_init;
N = N_eff + N_mem;
Y = Y_init;

%Initialise index for t_cell_matrix
rec_ind_tcell_matrix = N+1;

%Define the initial time to be zero
t=0;

%Gone is true if the pathogen is dead and false otherwise
gone = 0;

%Maximum number of cells:
max_t_cell = 0;
time_max_t_cell = 0;
time_0_t_cell = 0;
max_y_cell = 0;
time_max_y_cell = 0;
time_0_y_cell = 0;

while t < T_final
    
    %Increase the recording index
    rec_ind=rec_ind+1;
    
    if(gone==0)
        %Calculate Y
        Y = Y_init*exp(t*(alpha - N_eff*beta));
        Y = max(Y,0);
        if(Y < 10^(-6)) %condition that states when the pathogen is defeated
            Y = 0;
            gone = 1;
        end
    end
    
    %Fate decision for each T cell
    %Initialise indexes
    nCell=1;
    ind_N = 1;
    
    while nCell < rec_ind_tcell_matrix
        v_rand = rand(N,1)/N; %vector of N random numbers
       
        if(t_cell_matrix(nCell,1) == 1 || t_cell_matrix(nCell,1) == 2)
            rho = v_rand(ind_N);
            r_tau=rho*Y;
            ind_N = ind_N + 1;
        end
        
        %Killer T cell
        if(t_cell_matrix(nCell,1) == 1 || t_cell_matrix(nCell,1) == 3)
            if(t_cell_matrix(nCell,6) > 0)
                %In division phase
                t_cell_matrix(nCell,6) = max(t_cell_matrix(nCell,6)-t_next,0);
                
                %Division phase completed
                if(t_cell_matrix(nCell,6) == 0 && t_cell_matrix(nCell,1) == 3) 
                    N_eff = N_eff + 1;
                    t_cell_matrix(nCell,1) = 1;
                end
            else
                %Initial conditions
                p0_sys = t_cell_matrix(nCell,4);
                d0_sys = t_cell_matrix(nCell,5);
                c0_sys = t_cell_matrix(nCell,3);
                a0_sys = t_cell_matrix(nCell,2);
                
                %Explicit solutions for system 4.1
                [c,a,p,d] = sys_4_1_sol(t,lambda_taup,lambda_pp, r_tau, p0_sys, lambda_pd, d0_sys, mu_pc, c0_sys, mu_da, a0_sys);
                
                %Desision state
                if( a > 0 && c > 0)
                    d=max(d,0);
                    p=max(p,0);
                    t_cell_matrix(nCell,4) = p;
                    t_cell_matrix(nCell,5) = d;
                    t_cell_matrix(nCell,3) = c;
                    t_cell_matrix(nCell,2) = a;
                else
                    if(a <= 0)      %Initiate apoptosis
                        t_cell_matrix(nCell,6) = t_apo;
                        t_cell_matrix(nCell,1) = 4;
                        
                    elseif(c <= 0)  %Initiate division                       
                        %Membrane receptors are divided between 2 daughter
                        %cells
                        delta_P_child_1 = 0.4+(0.6-0.4)*rand();
                        delta_P_child_2 = 1 - delta_P_child_1;
                        delta_D_child_1 = 0.4+(0.6-0.4)*rand();
                        delta_D_child_2 = 1 - delta_D_child_1;
                        
                        r_p_child_1 = delta_P_child_1 * p;
                        r_p_child_2 = delta_P_child_2 * p;
                        
                        r_d_child_1 = delta_D_child_1 * d;
                        r_d_child_2 = delta_D_child_2 * d;
                        
                        %Actualization for daughter cells
                        t_cell_matrix(nCell,4) = r_p_child_1;
                        t_cell_matrix(nCell,5) = r_d_child_1;
                        t_cell_matrix(nCell,6) = t_cycle;
                        
                        t_cell_matrix(nCell,3) = c0;
                        t_cell_matrix(nCell,2) = a0;
                        
                        %type 3 -> new effector cell that has not
                        %completed division phase
                        t_cell_matrix(rec_ind_tcell_matrix,1) = 3;
                        t_cell_matrix(rec_ind_tcell_matrix,4) = r_p_child_2;
                        t_cell_matrix(rec_ind_tcell_matrix,5) = r_d_child_2;
                        t_cell_matrix(rec_ind_tcell_matrix,6) = t_cycle;
                        
                        t_cell_matrix(rec_ind_tcell_matrix,3) = c0;
                        t_cell_matrix(rec_ind_tcell_matrix,2) = a0;
                        
                        %Increase index for the next new cell
                        rec_ind_tcell_matrix = rec_ind_tcell_matrix + 1;
                    end
                end
            end
            %Next cell in population
            nCell = nCell + 1;
            
        %Memory T cell
        elseif(t_cell_matrix(nCell,1) == 2 || t_cell_matrix(nCell,1) == 5) 
            if(t_cell_matrix(nCell,6) > 0)
                %In division phase
                t_cell_matrix(nCell,6) = max(t_cell_matrix(nCell,6)-t_next,0);
                
                %Division phase completed
                if(t_cell_matrix(nCell,6)==0 && t_cell_matrix(nCell,1) == 5) 
                    N_mem=N_mem+1;
                    t_cell_matrix(nCell,1) =2;
                end
            else
                %Initial conditions
                c0_solsys =t_cell_matrix(nCell,3);
                p0_solsys =t_cell_matrix(nCell,4);
                
                %Explicit solutions for system 4.2
                [c,p] = sys_4_2_sol(t,mu_pc_mem, p0_solsys, lambda_taup_mem, lambda_pp_mem, r_tau, c0_solsys);
                
                %Division phase
                if(c <= 0)
                    delta_P_child_1 = 0.4+(0.6-0.4)*rand();
                    delta_P_child_2 = 1 - delta_P_child_1;
                   
                    r_p_child_1 = delta_P_child_1 * p;
                    r_p_child_2 = delta_P_child_2 * p;
                    
                    t_cell_matrix(nCell,4)=r_p_child_1;
                    t_cell_matrix(nCell,6)=t_cycle;
                    
                    t_cell_matrix(nCell,3) = c0_mem; 
                    
                    t_cell_matrix(rec_ind_tcell_matrix,1)=5;
                    t_cell_matrix(rec_ind_tcell_matrix,4)=r_p_child_2;
                    t_cell_matrix(rec_ind_tcell_matrix,6)=t_cycle;
                    
                    t_cell_matrix(rec_ind_tcell_matrix,3)=c0_mem;
                    
                    rec_ind_tcell_matrix = rec_ind_tcell_matrix + 1;
                else
                    t_cell_matrix(nCell,4)=p;
                    t_cell_matrix(nCell,3)=c;
      
                end
            end
            
            nCell=nCell+1;
           
            
        elseif(t_cell_matrix(nCell,1) == 4) %Effector T cell is dead
          
            if(t_cell_matrix(nCell,6) > 0)
                t_cell_matrix(nCell,6) = max(t_cell_matrix(nCell,6)-t_next,0);
                if(t_cell_matrix(nCell,6)==0) 
                    N_eff=N_eff-1;
                end
            end
            nCell=nCell+1;
        else
            break;
        end
        
    end 
    
    %Update the time
    t=t+t_next;
    
    %Record the time and the numbers of molecules
    time_vec(rec_ind) = t;
    N = N_eff + N_mem;
    rec_vector_N_eff(rec_ind) = N_eff;
    rec_vector_N_mem(rec_ind) = N_mem;
    rec_vector_Y(rec_ind) = Y;
    
    %Record the maximum number of cells
    if(N_eff > max_t_cell)
        max_t_cell = N_eff;
        time_max_t_cell = t;
    end
    
    if(Y > max_y_cell)
        max_y_cell = Y;
        time_max_y_cell = t;
    end
    
    if(N_eff == 0 && time_0_t_cell == 0)
        time_0_t_cell = t;
    end
    
    if(Y == 0 && time_0_y_cell == 0)
        time_0_y_cell = t;
    end
end

disp(max_y_cell)
disp(max_t_cell)

disp(time_max_y_cell)
disp(time_max_t_cell)

disp(time_0_y_cell)
disp(time_0_t_cell)

%Plot results
f1=figure;

figure(f1)
[hA1]=plot(time_vec,rec_vector_N_eff,'b','LineWidth', 1);

hold on
[hA2]=plot(time_vec,rec_vector_Y,'r','LineWidth', 1);


%Plot max number of cells 
[hM1] = plot(time_max_t_cell,max_t_cell,'o','MarkerFaceColor','b', 'MarkerEdgeColor', 'b');
[hM2] = plot(time_max_y_cell,max_y_cell,'o','MarkerFaceColor','r', 'MarkerEdgeColor', 'r');

[hM3] = plot(time_0_t_cell,0,'s','MarkerFaceColor','b', 'MarkerEdgeColor', 'b');
[hM4] = plot(time_0_y_cell,0,'s','MarkerFaceColor','r', 'MarkerEdgeColor', 'r');

legend([hA1,hA2,hM1,hM2,hM3,hM4],'Celulas T efectoras','Patogeno',...
    'Max celulas T', 'Max celulaspatogeno', 'Min celulas T', 'Min celulas patogeno');
xlabel('Tiempo');  ylabel('Número de celulas');
