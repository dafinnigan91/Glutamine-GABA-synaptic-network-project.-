import matplotlib.pyplot as plt 
import numpy as np
import random
import pandas as pd 
import sys
import json
import networkx as nx
  
P_r = float(sys.argv[1]) if len(sys.argv) > 1 else 0.8 # Default value for perameter sweep
           

R = 8.314  # Universal gas constant (J/(mol·K))
T = 310.15  # Temperature in Kelvin (37°C)
F = 96485  # Faraday's constant (C/mol)

# Permeabilities (Relative permeability of ions)
P_K = 1.0  
P_Na = 0.04 
P_Cl = 0.45 
P_Ca = 0.001  

class Glutamatergic_Synaps: # Glutamatergic synapse class for simulating AMPA and NMDA receptor dynamics
    def __init__(self):
        #AMPAR variables
        self.AMPAR_initial = 500 # level of AMPA receptors in post synaptic membrain
        self.AMPAR = self.AMPAR_initial #initial level of AMPA receptors in post synaptic membrain
        self.AMPAR_avalible = self.AMPAR #level of AMPA receptors in post synaptic membrain that are avalible for GlU to bind 
        self.AMPAR_bound = [] #list of bound AMPARS
        self.bound_A = 0 #number of AMPARs bound 
        self.AMPAR_insert = random.randrange(200, 300)
        #NMDAR variables 
        self.NMDAR = 50 # number of NMDA receptors in post synaptic membrain
        self.NMDAR_Mg = True #NMDA receptor blocker by Mg2+
        self.NMDAR_bound = [] # list of bound NMDARs
        self.bound_N = 0 # number of NMDARs bound
        self.GLu_mol_Vescl = 3000 #Glutamate per vesical 
        self.Vescl_doc = 2 #maximum docking cap of pre synaptic membrain
        #self.Vescl_R_prob = 0.9 #0.80 #MAIN LEVERAGE POINT FOR DYNAMICS''' # commented out for perameter sweep
        self.Vescl_R_prob = P_r # for peram sweep
        self.cleft_GLu_vol = 0 # initial level of glutamate in the cleft 
        self.V_replen_time = None #time it takes for vesical to replenish
        self.Ca_th_coef = 0.15 # coefficient for calcium threshold
        self.Ca_threshold = self.NMDAR * self.Ca_th_coef #MAIN LEVERAGE POINT FOR DYNAMICS
        self.Syn_Deprs_count = 0 # count of synaptic depression events 
        
        self.K_out = 5.0   # Extracellular K+
        self.K_in = 127.0  # Intracellular K+
        self.Na_out = 145.0  # Extracellular Na+
        self.Na_in = 15.0    # Intracellular Na+
        self.Cl_out = 110.0  # Extracellular Cl-
        self.Cl_in = 4.0     # Intracellular Cl-
        self.Ca_out = 0.1  # Extracellular Ca+2
        self.Ca_in = 0.0001     # Intracellular Ca+2
       
        self.Na_in_eq = 15.0   #equilibrium values
        self.K_in_eq = 140.0  
        self.K_out_eq = 5.0
        self.Cl_in_eq = 4.0  
        self.Cl_out_eq = 110.0  
        self.Ca_in_eq = 0.0001  
        self.Ca_out_eq = 1.8 

        self.Na_Chanel = False # Na+ channel is open state
        self.K_Chanel = False # K+ channel is open state
        
        self.AMPAR_history = [] # list to store AMPAR levels over time
        self.NMDAR_history = [] # list to store NMDAR levels over time
        self.Mem_Vlotage_History = [] # list to store membrane voltage over time
        self.Ca_history = [] # list to store calcium levels over time
        self.Na_in_History = [] # list to store Na+ levels over time
        self.Na_out_History = [] # list to store Na+ levels over time
        self.K_in_History = []  # list to store K+ levels over time
        self.K_out_History = [] # list to store K+ levels over time
        self.Cl_in_History = [] # list to store Cl- levels over time
        self.Cl_out_History = [] # list to store Cl- levels over time
        self.Ca_in_History = [] # list to store Ca2+ levels over time
        self.Ca_out_History = [] # list to store Ca2+ levels over time
        
        self.depression_threshold = self.Ca_threshold * 0.2 #MAIN LEVERAGE POINT FOR DYNAMICS and threshold for activity depression
        self.depression_duration = random.randint(500,1000) # ms of low calcium required to trigger depression

    def handle_spike(self, t): # Handle glutamate spike for glutamatergic synapse
        if self.Vescl_doc > 0:
            for _ in range(self.Vescl_doc):
                if random.random() <= self.Vescl_R_prob:
                    self.cleft_GLu_vol += self.GLu_mol_Vescl
            self.Vescl_doc -= 1
            self.V_replen_time = t

    def replenish_vesicles(self, t): # replenish vesicles for glutamatergic synapse
        if self.V_replen_time is not None and t - self.V_replen_time >= random.randint(1, 1):
            self.Vescl_doc = 1
            self.V_replen_time = None

    def handle_AMPAR_binding(self, t): # Handle AMPAR and NMDAR binding
        if self.cleft_GLu_vol > 0 and self.AMPAR_avalible > 0: # if glutamate is present in the cleft and AMPARs are available then bind glutamate
            self.bind_P = min(0.3, self.AMPAR_avalible / self.AMPAR) # binding probability is proportional to the level of available AMPARs
            self.scale_fac = (self.AMPAR / 30) # scaling factor for binding of ions
            self.max_binding_events = 100  # limit number of binding attempts per timestep
            self.binding_attempts = 0
            while self.cleft_GLu_vol > 0 and self.AMPAR_avalible > 0 and self.binding_attempts < self.max_binding_events:
                if random.random() < self.bind_P:
                    self.AMPAR_bound.append((t, self.AMPAR_avalible))  # record binding event
                    self.bound_A += 1
                    self.cleft_GLu_vol -= 1 # decrease glutamate in the cleft
                    self.AMPAR_avalible -= 1
                    self.Na_in += 2.0 * self.scale_fac # ions in and out after binding
                    self.Na_out -= 2.0 * self.scale_fac
                    self.K_out += 0.05 * self.scale_fac
                    self.K_in -= 0.05 * self.scale_fac
                    self.Cl_out -= 0.05 * self.scale_fac
                    self.Cl_in  += 0.05 * self.scale_fac
                    self.binding_attempts += 1
                if self.cleft_GLu_vol <= 0 or self.AMPAR_avalible <= 0:
                    break 
            self.Vol_AMPAR_bound = len(self.AMPAR_bound) # volum of AMPARs bound
            self.NMDAR_avalible = self.NMDAR # number of NMDARs available for binding
            if self.Vol_AMPAR_bound >= self.Ca_threshold:
                self.NMDAR_Mg = False
            if self.NMDAR_Mg == False and self.bound_N > 0:
                self.effective_P_Ca = P_Ca  # Allow Ca²⁺ through NMDARs
            else:
                self.effective_P_Ca = 0
            if self.NMDAR_Mg == False:
                self.attempts = 0
                self.max_attempts = 50
                while self.bound_N < self.NMDAR and self.cleft_GLu_vol > 0:
                    self.bind_P_NMDAR = self.NMDAR_avalible / self.NMDAR # binding probability is proportional to the level of available NMDARs
                    for Gl, NMDAR_RE in zip(range(1, int(self.cleft_GLu_vol) + 1), range(1, self.NMDAR_avalible + 1)):
                        #if random.random() < self.bind_P_NMDAR:
                            self.NMDAR_bound.append((Gl, NMDAR_RE))
                            self.bound_N += 1
                            self.cleft_GLu_vol -= 1 # decrease glutamate in the cleft
                            self.NMDAR_avalible -= 1
                            self.Na_in +=  1.0 # ions in and out after binding
                            self.Na_out -=  1.0
                            self.bind_P_NMDAR = self.NMDAR_avalible / self.NMDAR
                    if self.cleft_GLu_vol <= 0 or self.NMDAR_avalible <= 0:
                        break
                    self.attempts += 1
            self.Vol_NMDAR_bound = len(self.NMDAR_bound)
            self.Ca_in += self.Vol_NMDAR_bound
            self.Ca_out -= self.Vol_NMDAR_bound# Ca influx from NMDARs
                
        # Track how long Ca²⁺ has been low
        if self.Ca_in < self.depression_threshold:    
            self.Syn_Deprs_count += 1
        else:
            self.Syn_Deprs_count = 0  # reset if Ca spikes up again
                
    def remove_AMPAR_Ca(self):
    # If calcium stays low long enough, remove AMPARs
        if self.Syn_Deprs_count >= self.depression_duration and self.AMPAR > self.AMPAR_initial:
            self.AMPAR_ran_rem = random.randrange(50,100)
            self.AMPAR -= self.AMPAR_ran_rem 
            self.AMPAR_avalible = max(self.AMPAR_avalible - self.AMPAR_ran_rem , self.AMPAR_initial)
            self.Syn_Deprs_count = 0  
        return self.AMPAR
    
    def insert_AMPAR_Ca(self):
    # Plasticity rule: AMPAR insertion if Ca is high
        if self.Ca_in > self.Ca_threshold:
            self.delta_AMPAR = 5#int((self.Ca_in - self.Ca_threshold) * 4)
            self.AMPAR += self.delta_AMPAR
            self.AMPAR_avalible += self.delta_AMPAR
        return self.AMPAR
        
    def remove_AMPAR_vol(self):
        #Remove AMPAR if volum exceeds
        while self.AMPAR >= self.AMPAR_initial * 50: #MAIN LEVERAGE POINT FOR DYNAMICS 
            self.AMPAR -= random.randrange(5,50)
        return self.AMPAR
    
    def clamp_values(self): # function to clamp values to prevent runaway
        self.Na_in = min(max(self.Na_in, 0.01), 10000)
        self.Na_out = min(max(self.Na_out, 0.001), 50)
        self.K_in = min(max(self.K_in, 0.01), 15)
        self.K_out = min(max(self.K_out, 0.01), 2)
        self.Cl_in = min(max(self.Cl_in, 5), 50)
        self.Cl_out = min(max(self.Cl_out, 0.27), 95)
        self.Ca_in = min(max(self.Ca_in, 0.0001), 0.1)
        self.Ca_out = min(max(self.Ca_out, 0.01), 25)
        return self.Na_in, self.Na_out, self.K_in, self.K_out, self.Cl_in, self.Cl_out, self.Ca_in, self.Ca_out       
    
    
    def Ca_count(self):
    # Track how long Ca²⁺ has been low # if calcium is low
        if self.Ca_in < self.depression_threshold:    
            self.Syn_Deprs_count += 1
        else:
            self.Syn_Deprs_count = 0  # reset if Ca spikes up again      
    

    def Glutamate_appening(self): #Track ion volumes inside and out of the synaps at every step
        self.Na_in_History.append(self.Na_in)
        self.Na_out_History.append(self.Na_out)
        self.K_in_History.append(self.K_in)
        self.K_out_History.append(self.K_out)
        self.Cl_in_History.append(self.Cl_in)
        self.Cl_out_History.append(self.Cl_out)
        self.Ca_in_History.append(self.Ca_in)
        self.Ca_out_History.append(self.Ca_out)
        self.AMPAR_history.append(self.AMPAR) 
        self.NMDAR_history.append(self.NMDAR)
        
    def Ion_Pump(self): # Ion Pumps to stabelize membrain potential
        self.Na_pump_strength = (self.Na_in - self.Na_in_eq) * 0.02  
        self.K_pump_strength = (self.K_out - self.K_in_eq) * 0.02
        self.Na_out += self.Na_pump_strength * 10
        self.Na_in  -= self.Na_pump_strength * 10
        self.K_out  += self.K_pump_strength * 3
        self.K_in   -= self.K_pump_strength * 3
        return self.Na_in, self.Na_out, self.K_in, self.K_out
    
    def Ion_decay(self): # Passive ion deay with decay rates
        self.K_decay = 0.8
        self.Cl_decay = 0.9
        self.Ca_decay = 0.9
        self.K_in -= self.K_in * (1 - self.K_decay)
        self.K_out += self.K_in * (1 - self.K_decay)
        self.Cl_in -= self.Cl_in * (1 - self.Cl_decay)
        self.Cl_out += self.Cl_in * (1 - self.Cl_decay)
        self.Ca_in -= self.Ca_in * (1 - self.Ca_decay)
        self.Ca_out += self.Ca_in * (1 - self.Ca_decay)
        return self.K_in, self.K_out, self.Cl_in,self.Cl_out, self.Ca_in, self.Ca_out 

    def K_return(self, t): # Boost K+ outflow to return toward resting potential
        if t > 0 and Spike_train[t - 1] == 1 and Spike_train[t] == 0:
            self.K_in -= 5.0
            self.K_out += 5.0
        return self.K_in, self.K_out
    
    def close_channels(self, t): # Close channels if no stim
        if Spike_train[t] == 0:
            self.Na_Chanel = False  
            self.K_Chanel = False
            self.NMDAR_Mg = True

    def K_leak(self): # Stabilizing K+ leak toward resting potential
        self.K_in -= 0.005
        self.K_out += 0.005
        return self.K_in, self.K_out

    def Na_leak(self): # Stabilizing Na+ leak toward resting potential
        self.Na_in -= 0.002
        self.Na_out += 0.002
        return self.Na_in, self.Na_out

    def Cl_leak(self): # Stabilizing Cl- leak toward resting potential
        self.Cl_in -= 0.002
        self.Cl_out += 0.002
        return self.Cl_in, self.Cl_out 
        
    def Ca_leak(self): # Stabilizing Ca2+ leak toward resting potential
        self.Ca_in -= 0.0001
        self.Ca_out += 0.0001
        return self.Ca_in, self.Ca_out      

    #Glutamate cleared by transporters (at every time step, this is unrealistic)
    def Glu_leak(self):
        self.cleft_GLu_vol *= 0.9
        return self.cleft_GLu_vol
    
    def Membrain_potential(self): # Calculate membrane potential using a modified Goldman equation
        numerator = ((P_K * self.K_out) + (P_Na * self.Na_out) + (P_Cl * self.Cl_in) + (P_Ca * self.Ca_out))
        denominator = ((P_K * self.K_in) + (P_Na * self.Na_in) + (P_Cl * self.Cl_out) + (P_Ca * self.Ca_in))

        
        numerator = max(numerator, 1e-8) # Prevent divide-by-zero or log of zero
        denominator = max(denominator, 1e-8) # Prevent divide-by-zero or log of zero
        Membrain_potential = Vm = (R * T / F) * np.log(numerator / denominator) # Membrane potential in volts
        return Membrain_potential
    
    def GLU_step(self, t, per_spike): # Main function to update the synapse state at each time step
        if per_spike:
            self.handle_spike(t)
        self.replenish_vesicles(t)
        self.handle_AMPAR_binding(t)
        self.remove_AMPAR_Ca()
        self.insert_AMPAR_Ca()
        self.remove_AMPAR_vol()
        self.clamp_values()
        self.Ca_count()
        self.clamp_values()
        self.Glutamate_appening()
        self.Ion_Pump()
        self.Ion_decay()
        self.K_return(t)
        self.close_channels(t)
        self.K_leak()
        self.Na_leak()
        self.Cl_leak()
        self.Ca_leak()
        self.Glu_leak()
        Vm = self.Membrain_potential()
        self.Mem_Vlotage_History.append(Vm)
        return Vm # return the membrane potential
 
##################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################
class GABAergic_Synaps:#GABAergic synapse class for simulating GABA receptor dynamics
    def __init__(self):
        #GABAergic variables
        self.GABA_K_out = 5.0   # Extracellular K+
        self.GABA_K_in = 127.0  # Intracellular K+
        self.GABA_Na_out = 145.0  # Extracellular Na+
        self.GABA_Na_in = 15.0    # Intracellular Na+
        self.GABA_Cl_out = 110.0  # Extracellular Cl-
        self.GABA_Cl_in = 4.0     # Intracellular Cl-
        self.GABA_Ca_out = 0.1  # Extracellular Ca+2
        self.GABA_Ca_in = 0.0001     # Intracellular Ca+2
        
        self.GABA_Na_in_eq = 15.0 #equilibrium values
        self.GABA_Na_out_eq = 145.0 
        self.GABA_K_in_eq = 140.0 
        self.GABA_K_out_eq = 5.0 
        self.GABA_Cl_in_eq = 4.0 
        self.GABA_Cl_out_eq = 110.0
        self.GABA_Ca_in_eq = 0.001
        self.GABA_Ca_out_eq = 100.0 

        self.GABA_mol_Vescl = 3000 #Glutamate per vesical 
        self.GABA_Vescl_doc = 2 #maximum docking cap of pre synaptic membrain
        #self.GABA_Vescl_R_prob = 0.80 #MAIN LEVERAGE POINT FOR DYNAMICS''' # commented out for perameter sweep
        self.GABA_Vescl_R_prob = P_r # for peram sweep
        self.GABA_cleft_vol = 0 # initial level of glutamate in the cleft 
        self.GABA_V_replen_time = None #time it takes for vesical to replenish 

        self.GABA_A_initial = 500 # level of AMPA receptors in post synaptic membrain
        self.GABA_A = self.GABA_A_initial #initial level of AMPA receptors in post synaptic membrain
        self.GABA_A_avalible = self.GABA_A #level of AMPA receptors in post synaptic membrain that are avalible for GlU to bind 
        self.GABA_A_bound = [] #list of bound AMPARS
        self.bound_GABA_A = 0 #number of AMPARs bound 
        self.GABA_A_insert = random.randrange(200, 300) 

        self.GABA_B_initial = 500 # level of AMPA receptors in post synaptic membrain
        self.GABA_B = self.GABA_B_initial #initial level of AMPA receptors in post synaptic membrain
        self.GABA_B_avalible = self.GABA_B #level of AMPA receptors in post synaptic membrain that are avalible for GlU to bind 
        self.GABA_B_bound = [] #list of bound AMPARS
        self.bound_GABA_B = 0 #number of AMPARs bound 
        self.GABA_B_insert = random.randrange(200, 300) # GaBA_B receptor insertion rate
        self.GABA_B_threshold = 200 # threshold for GABA_B activation
        self.GABA_B_decay_counter = 0 # decay counter for GABA_B
        self.GABA_B_active = False # GABA_B activation state
        self.cleft_GABA_vol = 0 # initial level of glutamate in the cleft
#
        self.GABA_Ca_th_coef = 0.8 # coefficient for calcium threshold
        self.GABA_Ca_threshold = self.GABA_B_initial * self.GABA_Ca_th_coef #MAIN LEVERAGE POINT FOR DYNAMICS
        self.GABA_LTD_counter = 0 # counter for long-term depression
        #NMDAR * self.Ca_th_coef #MAIN LEVERAGE POINT FOR DYNAMICS
        
        self.depression_threshold = self.GABA_Ca_threshold * 0.2 #MAIN LEVERAGE POINT FOR DYNAMICS and threshold for activity depression
        self.depression_duration = random.randint(500,1000) # ms of low calcium required to trigger depression
        self.Syn_Deprs_count = 0 # count of synaptic depression events
        
        self.GABA_Ar_history = [] # list to store AMPAR levels over time
        self.GABA_Br_history = [] # list to store NMDAR levels over time
        self.GABA_Mem_Vlotage_History = [] # list to store membrane voltage over time
        self.GABA_Ca_history = [] # list to store calcium levels over time
        self.GABA_Na_in_History = [] # list to store Na+ levels over time
        self.GABA_Na_out_History = [] # list to store Na+ levels over time
        self.GABA_K_in_History = [] # list to store K+ levels over time
        self.GABA_K_out_History = [] # list to store K+ levels over time
        self.GABA_Cl_in_History = [] # list to store Cl- levels over time
        self.GABA_Cl_out_History = [] # list to store Cl- levels over time
        self.GABA_Ca_in_History = [] # list to store Ca2+ levels over time
        self.GABA_Ca_out_History = [] # list to store Ca2+ levels over time

    def handle_spike_GABA(self, t): # Hhandle GABA spike for GABAergic synapse
        if self.GABA_Vescl_doc > 0: 
            for _ in range(self.GABA_Vescl_doc):
                if random.random() <= self.GABA_Vescl_R_prob:
                    self.GABA_cleft_vol += self.GABA_mol_Vescl
            self.GABA_Vescl_doc -= 1
            self.GABA_V_replen_time = t

    def replenish_vesicles_GABA(self, t): # replenish vesicles for GABAergic synapse
        if self.GABA_V_replen_time is not None and t - self.GABA_V_replen_time >= random.randint(1, 1):
            self.GABA_Vescl_doc = 1
            self.GABA_V_replen_time = None

    def GABA_A_GABA_B_binding_sequence(self, t):  # GABA_A and GABA_B binding sequence
        if self.GABA_cleft_vol > 0 and self.GABA_A_avalible > 0:
            bind_P = min(0.3, self.GABA_A_avalible / self.GABA_A)
            scale_fac = self.GABA_A / 30
            self.max_binding_events = 50
            self.binding_attempts = 0
            while self.GABA_cleft_vol > 0 and self.GABA_A_avalible > 0 and self.binding_attempts < self.max_binding_events:
                if random.random() < bind_P: # binding probability
                    self.GABA_A_bound.append((t, self.GABA_A_avalible)) # record binding event
                    self.bound_GABA_A += 1
                    self.GABA_cleft_vol -= 1
                    self.GABA_A_avalible -= 1
                    self.GABA_Cl_out -= 0.05 * scale_fac
                    self.GABA_Cl_in += 0.05 * scale_fac
                    self.GABA_Ca_in += 10 * scale_fac
                    self.GABA_Ca_out -= 10 * scale_fac
                if self.GABA_cleft_vol <= 0 or self.GABA_A_avalible <= 0:
                    break

        if self.GABA_cleft_vol > self.GABA_B_threshold and self.GABA_B_decay_counter == 0: # if GABA_B threshold is reached then activate GABA_B
            self.GABA_B_active = True
            self.GABA_B_scaling = 1.0
            self.GABA_B_decay = 0.99

        if self.GABA_B_active:
            self.GABA_K_out += 30 * self.GABA_B_scaling  # Increase K+ outflow to amplify inhibition
            self.GABA_K_in -= 30 * self.GABA_B_scaling
            self.GABA_B_decay_counter += 1
            self.GABA_B_scaling *= self.GABA_B_decay
            if self.GABA_B_scaling < 0.01:
                self.GABA_B_active = False
                self.GABA_B_decay_counter = 0

    def GABA_B_enable(self):    # Track GABA_B activation 
        if self.GABA_Ca_in > self.GABA_Ca_threshold: # If calcium is high, activate GABA_B
            self.GABA_LTP_enabled = True # Enable GABA LTP
        elif self.GABA_Ca_in < (self.GABA_Ca_threshold * 0.5): # If calcium is low, deactivate GABA_B
            self.GABA_LTP_enabled = False # Deactivate GABA LTP
        return self.GABA_LTP_enabled # return the state of GABA LTP
    
    def GABA_insert(self):
        if self.GABA_LTP_enabled and self.GABA_Ca_in > self.GABA_Ca_threshold: # Track high calcium for GABA LTP
            self.delta_GABA_A = int((self.GABA_Ca_in - self.GABA_Ca_threshold) * 3)
            self.GABA_A += 10 #self.delta_GABA_A # Increase GABA_A receptors
            self.GABA_A_avalible += self.delta_GABA_A # Increase available GABA_A receptors
            self.GABA_A = min(self.GABA_A, self.GABA_A_initial * 50) # Clamp to maximum
            self.GABA_A_avalible = min(self.GABA_A_avalible, self.GABA_A_initial * 50) # Clamp to maximum available GABA_A receptors   
            self.GABA_Vescl_R_prob = min(self.GABA_Vescl_R_prob + 0.01, 1.0) # Increase vesicle release probability
        
    def remove_GABA_A(self):
        if self.GABA_LTD_counter > 30: # If calcium is low for a long time, remove GABA_A receptors
            GABA_rmv = 5 #random.randint(5, 10)  # Randomly remove GABA_A receptors
            self.GABA_A -= GABA_rmv # Decrease GABA_A receptors
            self.GABA_A_avalible = max(self.GABA_A_avalible - GABA_rmv, self.GABA_A_initial)   # Decrease available GABA_A receptors
            self.GABA_Vescl_R_prob = max(self.GABA_Vescl_R_prob - 0.01, 0.1)  # Decrease vesicle release probability
            self.GABA_LTD_counter = 0
            
    def remove_GABA_A_volume(self): #Remove GABA if volum exceeds their initial value
        while self.GABA_A >= self.GABA_A_initial * 50: #MAIN LEVERAGE POINT FOR DYNAMICS 
            self.GABA_A -= random.randrange(5,50)
        return self.GABA_A
            
    def GABA_count(self): # Track how long Ca²⁺ has been low and trigger GABA long-term depression
        if self.GABA_Ca_in < self.GABA_Ca_threshold: #
            self.GABA_LTD_counter += 1
        else:
            self.GABA_LTD_counter = 0  # reset if Ca spikes up again
        return self.GABA_LTD_counter,
            
    def GABA_clamp_values(self): # Clamp values to prevent run
        self.GABA_Na_in = min(max(self.GABA_Na_in, 0.01), 10000)
        self.GABA_Na_out = min(max(self.GABA_Na_out, 0.001), 50)
        self.GABA_K_in = min(max(self.GABA_K_in, 0.01), 15)
        self.GABA_K_out = min(max(self.GABA_K_out, 0.01), 2)
        self.GABA_Cl_in = min(max(self.GABA_Cl_in, 5), 50)
        self.GABA_Cl_out = min(max(self.GABA_Cl_out, 0.27), 95)
        self.GABA_Ca_in = min(max(self.GABA_Ca_in, 0.0001), 0.1)
        self.GABA_Ca_out = min(max(self.GABA_Ca_out, 0.01), 25)
        self.GABA_A = max(self.GABA_A, self.GABA_A_initial)
        return self.GABA_Na_in, self.GABA_Na_out, self.GABA_K_in, self.GABA_K_out, self.GABA_Cl_in, self.GABA_Cl_out, self.GABA_Ca_in, self.GABA_Ca_out, self.GABA_A        
            
    def GABA_Ca_count(self): # Track how long Ca²⁺ has been low
        if self.GABA_Ca_in < self.depression_threshold:    
            self.Syn_Deprs_count += 1
        else:
            self.Syn_Deprs_count = 0  # reset if Ca spikes up again
        return self.Syn_Deprs_count  # return the count of depression
    
    def GABA_appending(self): #Track ion volumes inside and out of the synaps at every step 
        self.GABA_Na_in_History.append(self.GABA_Na_in)
        self.GABA_Na_out_History.append(self.GABA_Na_out)
        self.GABA_K_in_History.append(self.GABA_K_in)
        self.GABA_K_out_History.append(self.GABA_K_out)
        self.GABA_Cl_in_History.append(self.GABA_Cl_in)
        self.GABA_Cl_out_History.append(self.GABA_Cl_out)
        self.GABA_Ca_in_History.append(self.GABA_Ca_in)
        self.GABA_Ca_out_History.append(self.GABA_Ca_out)
        self.GABA_Ar_history.append(self.GABA_A)
        self.GABA_Br_history.append(self.GABA_B) 

    def GABA_Ion_Pump(self): #Ion Pumps to stabelize membrain potential 
        self.Na_pump_strength = (self.GABA_Na_in - self.GABA_Na_in_eq) * 0.05  
        self.K_pump_strength = (self.GABA_K_out - self.GABA_K_in_eq) * 0.05  
        self.GABA_Na_out += self.Na_pump_strength * 10
        self.GABA_Na_in  -= self.Na_pump_strength * 10
        self.GABA_K_out  += self.K_pump_strength * 3
        self.GABA_K_in   -= self.K_pump_strength * 3
        return self.GABA_Na_in, self.GABA_Na_out, self.GABA_K_in, self.GABA_K_out
    
    
    def GABA_Ion_decay(self): # Passive ion deay with decay rates
        self.GABA_K_decay = 0.8
        self.GABA_Cl_decay = 0.9
        self.GABA_Ca_decay = 0.9
        self.GABA_K_in -= self.GABA_K_in * (1 - self.GABA_K_decay)
        self.GABA_K_out += self.GABA_K_in * (1 - self.GABA_K_decay)
        self.GABA_Cl_in -= self.GABA_Cl_in * (1 - self.GABA_Cl_decay)
        self.GABA_Cl_out += self.GABA_Cl_in * (1 - self.GABA_Cl_decay)
        self.GABA_Ca_in -= self.GABA_Ca_in * (1 - self.GABA_Ca_decay)
        self.GABA_Ca_out += self.GABA_Ca_in * (1 - self.GABA_Ca_decay)
        return self.GABA_K_in, self.GABA_K_out, self.GABA_Cl_in,self.GABA_Cl_out, self.GABA_Ca_in, self.GABA_Ca_out

    def GABA_K_return(self, t): # Boost K+ outflow to return toward resting potential
        if t > 0 and Spike_train[t - 1] == 1 and Spike_train[t] == 0:
            self.GABA_K_in -= 5.0
            self.GABA_K_out += 5.0
        return self.GABA_K_in, self.GABA_K_out

    
    def GABA_close_channels(self, t): #Close channels if no stim
        if Spike_train[t] == 0:
            self.GABA_Na_Chanel = False
            self.GABA_K_Chanel = False

    
    def GABA_K_leak(self): #Stabilizing K+ leak toward resting potential
        self.GABA_K_in -= 0.01
        self.GABA_K_out += 0.01
        return self.GABA_K_in, self.GABA_K_out

    def GABA_Na_leak(self): #Stabilizing Na+ leak toward resting potential
        self.GABA_Na_in -= 0.01
        self.GABA_Na_out += 0.01
        return self.GABA_Na_in, self.GABA_Na_out

    def GABA_Cl_leak(self): #Stabilizing Cl- leak toward resting potential
        self.GABA_Cl_in -= 0.01
        self.GABA_Cl_out += 0.01
        return self.GABA_Cl_in, self.GABA_Cl_out 
        
    def GABA_Ca_leak(self): #Stabilizing Ca2+ leak toward resting potential
        self.GABA_Ca_in -= 0.01
        self.GABA_Ca_out += 0.01
        return self.GABA_Ca_in, self.GABA_Ca_out
       
    def GABA_leak(self): #GABAergic leak toward zero
        self.cleft_GABA_vol *= 0.9
        return self.cleft_GABA_vol
    
    def Membrain_potential_GABA(self): # Calculate membrane potential using a modified Goldman equation
        numerator = ((P_K * self.GABA_K_out) + (P_Na * self.GABA_Na_out) + (P_Cl * self.GABA_Cl_in) + (P_Ca * self.GABA_Ca_out))
        denominator = ((P_K * self.GABA_K_in) + (P_Na * self.GABA_Na_in) + (P_Cl * self.GABA_Cl_out) + (P_Ca * self.GABA_Ca_in))

        
        numerator = max(numerator, 1e-8) # Prevent divide-by-zero or log of zero
        denominator = max(denominator, 1e-8)
        Membrain_potential = Vm = ((R * T / F) * np.log(numerator / denominator)) 
        return Membrain_potential
    
    def GABA_step(self, t, per_spike): # Main function to update the synapse state at each time step
        if per_spike:
            self.handle_spike_GABA(t)
        self.replenish_vesicles_GABA(t)
        self.GABA_A_GABA_B_binding_sequence(t)
        self.GABA_B_enable()
        self.GABA_insert()
        self.remove_GABA_A()
        self.remove_GABA_A_volume()
        self.GABA_count()
        self.GABA_clamp_values()
        self.GABA_Ca_count()
        self.GABA_Ion_Pump()
        self.GABA_Ion_decay()
        self.GABA_K_return(t)
        self.GABA_close_channels(t)
        self.GABA_K_leak()
        self.GABA_Na_leak()
        self.GABA_Cl_leak()
        self.GABA_Ca_leak()
        self.GABA_leak()
        self.GABA_appending()
        Vm = self.Membrain_potential_GABA()
        self.GABA_Mem_Vlotage_History.append(Vm)
        return Vm
        
Spike_train = np.zeros(9999, dtype=int) #MAIN LEVERAGE POINT FOR DYNAMICS, pluss initialization of spike train
Spike_train[500:1000:25] = 1 
Spike_train[2000:2500:25] = 1 
Spike_train[3500:4000:25] = 1 
Spike_train[random.randint(3500,6000)] = 1
spike_times = [i for i, x in enumerate(Spike_train) if x == 1]
Spike_History = []      
        
class Neuron: # Neuron class to represent a neuron in the network
    def __init__(self, name):
        self.name = name
        self.spike_train = Spike_train.copy()  # 100 timesteps
        self.Vm = -70  # Resting membrane potential
        self.Vm_history = [] # list to store Vm history
        self.synapses = [] # list to store synapses
        self.AMPAR_history = [] # list to store AMPAR history
        self.NMDAR_history = [] # list to store NMDAR history
        self.GABA_A_history = [] # list to store GABA_A history
        self.GABA_B_history = [] # list to store GABA_B history
    
    def add_synapse(self, synapse):
        self.synapses.append(synapse)

    def update_spike_train(self, t, spike):
        # Update spike train at timestep `t`
        self.spike_train[t] = spike
        
##############################################################################################################################################################################################

Net = nx.DiGraph() # Create a directed graph to represent the network (note: update to undirected graph for retro signaling)

for name in ["N1", "N2", "N3", "N4"]: # create some neurons and add them to the network
#for name in ["N1", "N2", "N3", "N4", "N5"]: # create some neurons and add them to the network
    Net.add_node(name, neuron=Neuron(name)) # 

def add_synapse_edge(net, source, target, synapse_class): # Add a synapse edge between two neurons
    synapse = synapse_class() # Create an instance of the synapse class
    net.add_edge(source, target, synapse=synapse) # Add the synapse to the source neuron
    net.nodes[source]["neuron"].add_synapse(synapse) # Add the synapse to the target neuron
    
add_synapse_edge(Net, "N1", "N2", Glutamatergic_Synaps) # Add synapse edges between neurons
add_synapse_edge(Net, "N1", "N3", Glutamatergic_Synaps)
add_synapse_edge(Net, "N2", "N4", GABAergic_Synaps)
add_synapse_edge(Net, "N3", "N4", Glutamatergic_Synaps)
#add_synapse_edge(Net, "N4", "N5", Glutamatergic_Synaps)

for t in range(len(Spike_train)): # Iterate over time steps from 0 to 9999
    for u, v, data in Net.edges(data=True): # Iterate over all edges in the network
        source_neuron = Net.nodes[u]["neuron"] # representation of the source neuron
        synapse = data["synapse"] # representations of the synapse

        if source_neuron.name == "N1": # if the source neuron is N1 
            per_spike = source_neuron.spike_train[t] == 1 # and the spike train is 1
        else: # else if the source neuron is not N1
            per_spike = source_neuron.Vm >= -0.60 # and the Vm is greater than -0.60

        if isinstance(synapse, Glutamatergic_Synaps): # if the synapse is Glutamatergic
            synapse.GLU_step(t, per_spike) # update the synapse state so that it can be used in the next time step
        elif isinstance(synapse, GABAergic_Synaps): # if the synapse is GABAergic
            synapse.GABA_step(t, per_spike) # update the synapse state so that it can be used in the next time step
            

    for node in Net.nodes: # Iterate over all nodes in the network
        neuron = Net.nodes[node]["neuron"] # Get the neuron object for the current node

        if neuron.name == "N4": # Custom logic: N4 Vm = output from N3 (GLU) - output from N2 (GABA)
            excitatory_V = 0 
            inhibitory_V = 0

            for u, v, data in Net.in_edges(node, data=True): # Iterate over incoming edges to N4 in order to calculate the Vm
                synapse = data["synapse"]
                if u == "N3" and isinstance(synapse, Glutamatergic_Synaps):
                    excitatory_V = synapse.Membrain_potential()
                elif u == "N2" and isinstance(synapse, GABAergic_Synaps):
                    current_inhibitory_V = synapse.Membrain_potential_GABA()
                    if len(synapse.GABA_Mem_Vlotage_History) > 1: # if the GABA_Mem_Vlotage_History has more than one value
                        previous_inhibitory_V = synapse.GABA_Mem_Vlotage_History[-1] # the previous value is the last value in the list 
                        if current_inhibitory_V < previous_inhibitory_V:  # Apply only if inhibitory potential is decreasing
                            inhibitory_V = current_inhibitory_V
                            neuron.Vm = excitatory_V - (5 * inhibitory_V) # apply a scaling factor to the inhibitory potential MAIN LEVERAGE POINT FOR DYNAMICS
                    else:
                        inhibitory_V = current_inhibitory_V  # No history, use current value

                        neuron.Vm = excitatory_V - inhibitory_V

        else:
            # Default logic for N1, N2, N3
            incoming_voltages = []
            for u, v, data in Net.in_edges(node, data=True):
                synapse = data["synapse"]
                if isinstance(synapse, Glutamatergic_Synaps):
                    incoming_voltages.append(synapse.Membrain_potential())
                elif isinstance(synapse, GABAergic_Synaps):
                    incoming_voltages.append(synapse.Membrain_potential_GABA())

            if incoming_voltages:
                neuron.Vm = np.mean(incoming_voltages)
            else:
                neuron.Vm = neuron.Vm_history[-1] if neuron.Vm_history else -0.70

        # Log membrane voltage
        neuron.Vm_history.append(neuron.Vm) # the 
        
        
        # Track receptor histories
        glut_found = False
        gaba_found = False
        for syn in neuron.synapses:
            if isinstance(syn, Glutamatergic_Synaps) and not glut_found:
                neuron.AMPAR_history.append(syn.AMPAR)
                neuron.NMDAR_history.append(syn.NMDAR)
                glut_found = True
            elif isinstance(syn, GABAergic_Synaps) and not gaba_found:
                neuron.GABA_A_history.append(syn.GABA_A)
                neuron.GABA_B_history.append(syn.GABA_B)
                gaba_found = True
                    
    glut_found = False
    gaba_found = False
    if isinstance(syn, Glutamatergic_Synaps) and not glut_found:
        neuron.AMPAR_history.append(syn.AMPAR)
        neuron.NMDAR_history.append(syn.NMDAR)
        glut_found = True
    elif isinstance(syn, GABAergic_Synaps) and not gaba_found:
        neuron.GABA_A_history.append(syn.GABA_A)
        neuron.GABA_B_history.append(syn.GABA_B)
        gaba_found = True
        
for node in Net.nodes:
    neuron = Net.nodes[node]["neuron"]
    filename = f"Vm_history_{neuron.name}_{P_r:.2f}.npy"
    
    filename = f"C:\\Users\\david\\OneDrive\\Desktop\\PYTHONS\\New folder\\Vm_history_{neuron.name}_{P_r:.2f}.npy"
    np.save(filename, np.array(neuron.Vm_history))

 # REMOVE EDITOUT OPERATIRS AS AND WEHN CODE IS NEEDED TO BE PLOTTED AND ANALYSED   
'''plt.figure(figsize=(12, 6))
for u, v, data in Net.edges(data=True):
    synapse = data["synapse"]
    label = f"{u}→{v}"
    if hasattr(synapse, 'AMPAR_history') and synapse.AMPAR_history:
        plt.plot(synapse.AMPAR_history, label=f"{label} AMPAR")
    #if hasattr(synapse, 'NMDAR_history') and synapse.NMDAR_history:
     #   plt.plot(synapse.NMDAR_history, label=f"{label} NMDAR")
    if hasattr(synapse, 'GABA_Ar_history') and synapse.GABA_Ar_history:
        plt.plot(synapse.GABA_Ar_history, label=f"{label} GABA_A")
    #if hasattr(synapse, 'GABA_Br_history') and synapse.GABA_Br_history:
     #   plt.plot(synapse.GABA_Br_history, label=f"{label} GABA_B")

plt.xlabel("Timestep")
plt.ylabel("Receptor Volume")
plt.title("Receptor Dynamics at Synapse Edges")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

for node in ["N2", "N3", "N4"]:
#for node in ["N2", "N3", "N4", "N5"]:
    neuron = Net.nodes[node]["neuron"]
    plt.figure(figsize=(8, 4))
    plt.plot(neuron.Vm_history, label=f"{node} Vm", color='tab:blue')
    plt.xlabel("Timestep")
    plt.ylabel("Membrane Potential (V)")
    plt.title(f"Membrane Potential of {node} Over Time")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

##############################################################################################################################################################################################
def compute_lyapunov_exponent_for_binding_prob(initial_binding_prob, delta=1e-5, timesteps=9999):
    def simulate_system(binding_prob):
        # Update the binding probability for all synapses
        for u, v, data in Net.edges(data=True):
            synapse = data["synapse"]
            if isinstance(synapse, Glutamatergic_Synaps):
                synapse.Vescl_R_prob = binding_prob
            elif isinstance(synapse, GABAergic_Synaps):
                synapse.GABA_Vescl_R_prob = binding_prob

        # Run the simulation and collect a trajectory (e.g., average calcium levels or spike rates)
        trajectory = []
        for t in range(timesteps):
            for u, v, data in Net.edges(data=True):
                synapse = data["synapse"]
                if isinstance(synapse, Glutamatergic_Synaps):
                    synapse.GLU_step(t, per_spike=False)  # Simulate without external spikes
                elif isinstance(synapse, GABAergic_Synaps):
                    synapse.GABA_step(t, per_spike=False)
            # Collect a metric (e.g., average membrane voltage across synapses)
            avg_membrainV = np.mean([synapse.Mem_Vlotage_History for _, _, data in Net.edges(data=True)
                                 if isinstance(data["synapse"], Glutamatergic_Synaps)])
            trajectory.append(avg_membrainV)
        return np.array(trajectory)

    # Simulate the system for the initial and perturbed values of the binding probability
    trajectory1 = simulate_system(initial_binding_prob)
    trajectory2 = simulate_system(initial_binding_prob + delta)

    # Calculate the divergence between the trajectories
    divergence = np.abs(trajectory1 - trajectory2)
    divergence[divergence == 0] = 1e-8  # Avoid log(0)

    # Compute the Lyapunov exponent
    lyapunov_exponent = np.mean(np.log(divergence / delta))
    return lyapunov_exponent

initial_binding_prob = 0.01  # Example initial value for Vescl_R_prob
lyapunov_exponent = compute_lyapunov_exponent_for_binding_prob(initial_binding_prob)
print(f"Lyapunov Exponent for Vescl_R_prob = {initial_binding_prob}: {lyapunov_exponent}")

binding_prob_values = np.linspace(0.01, 1.0, 5)  # Range of binding probabilities
lyapunov_exponents = []

for binding_prob in binding_prob_values:
    exponent = compute_lyapunov_exponent_for_binding_prob(binding_prob)
    lyapunov_exponents.append(exponent)

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(binding_prob_values, lyapunov_exponents, marker='o', color='blue')
plt.xlabel("Binding Probability (Vescl_R_prob)")
plt.ylabel("Lyapunov Exponent")
plt.title("Lyapunov Exponent vs. Binding Probability")
plt.grid(True)
plt.tight_layout()
plt.show()
'''