import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd 
import sys
import json

P_r = float(sys.argv[1]) if len(sys.argv) > 1 else 0.8

# membrain potentail is the Goldman Hodgkin Katz constant feild equasion
def Membrain_potential(P_K, P_Na, P_Cl, P_Ca, K_in, K_out, Na_in, Na_out, Cl_in, Cl_out, Ca_in, Ca_out):
    Membrain_potential = (R * T / F) * np.log(
                    ((P_K * K_out) + (P_Na * Na_out) + (P_Cl * Cl_in) + (P_Ca * Ca_out)) /
                    ((P_K * K_in) + (P_Na * Na_in) + (P_Cl * Cl_out) + (P_Ca * Ca_in)))
    Mem_Vlotage_History.append(Membrain_potential)
    return Membrain_potential
# Membrane potential            
#Physical perameters
R = 8.314  # Universal gas constant (J/(mol·K))
T = 310.15  # Temperature in Kelvin (37°C)
F = 96485  # Faraday's constant (C/mol)

# Permeabilities (Relative permeability of ions)
P_K = 1.0  # Relative permeability of K+
P_Na = 0.04  # Relative permeability of Na+
P_Cl = 0.45  # Relative permeability of Cl-
P_Ca = 0.001  # Relative permeability of Ca+2

AMPAR_history = []
NMDAR_history = []
Mem_Vlotage_History = []
Ca_history = []
Na_in_History = []
Na_out_History = []
K_in_History = []
K_out_History = []
Cl_in_History = []
Cl_out_History = []
Ca_in_History = []
Ca_out_History = []

#Setup
Spike_train = np.zeros(9999, dtype=int) #MAIN LEVERAGE POINT FOR DYNAMICS
Spike_train[500:1000:25] = 1 #MAIN LEVERAGE POINT FOR DYNAMICS
Spike_train[2000:2500:25] = 1 #MAIN LEVERAGE POINT FOR DYNAMICS
Spike_train[3500:4000:25] = 1 #MAIN LEVERAGE POINT FOR DYNAMICS
Spike_train[random.randint(3500,6000)] = 1
spike_times = [i for i, x in enumerate(Spike_train) if x == 1]
Spike_History = []

class Gultamatergic_Synaps:
    def __init__(self, AMPAR, NMDAR, Vescl_R_prob):
        self.AMPAR = AMPAR
        self.NMDAR = NMDAR
        self.Vescl_R_prob = Vescl_R_prob
        #AMPAR variables
        self.AMPAR_initial = 500 # level of AMPA receptors in post synaptic membrain
        self.AMPAR = self.AMPAR_initial #initial level of AMPA receptors in post synaptic membrain
        self.AMPAR_avalible = AMPAR #level of AMPA receptors in post synaptic membrain that are avalible for GlU to bind 
        self.AMPAR_bound = [] #list of bound AMPARS
        self.bound_A = 0 #number of AMPARs bound 
        self.AMPAR_insert = random.randrange(200, 300)
        #NMDAR variables 
        self.NMDAR = 50
        self.NMDAR_Mg = True
        self.NMDAR_bound = []
        self.bound_N = 0
        #Vesical veariables 
        self.GLu_mol_Vescl = 3000 #Glutamate per vesical 
        self.Vescl_doc = 2 #maximum docking cap of pre synaptic membrain
        self.Vescl_R_prob = 0.80 #MAIN LEVERAGE POINT FOR DYNAMICS''' # commented out for perameter sweep
        #self.Vescl_R_prob = P_r # for peram sweep
        self.cleft_GLu_vol = 0 # initial level of glutamate in the cleft 
        self.V_replen_time = None #time it takes for vesical to replenish
        self.Ca_th_coef = 0.8
        self.Ca_threshold = NMDAR * self.Ca_th_coef #MAIN LEVERAGE POINT FOR DYNAMICS
        self.Syn_Deprs_count = 0
        
        self.K_out = 5.0   # Extracellular K+
        self.K_in = 127.0  # Intracellular K+

        self.Na_out = 145.0  # Extracellular Na+
        self.Na_in = 15.0    # Intracellular Na+

        self.Cl_out = 110.0  # Extracellular Cl-
        self.Cl_in = 4.0     # Intracellular Cl-

        self.Ca_out = 0.1  # Extracellular Ca+2
        self.Ca_in = 0.0001     # Intracellular Ca+2
        #equilibrium values
        self.Na_in_eq = 15.0 # equilibrium values for Na
        self.Na_out_eq = 145.0 # equilibrium values for Na

        self.K_in_eq = 140.0 # equilibrium values for K
        self.K_out_eq = 5.0 # equilibrium values for K

        self.Cl_in_eq = 4.0 # equilibrium values for Cl
        self.Cl_out_eq = 110.0 # equilibrium values for cl

        self.Ca_in_eq = 0.001 # equilibrium values for Na
        self.Ca_out_eq = 100.0 # equilibrium values for Na

        self.Na_Chanel = False
        self.K_Chanel = False
        
        self.depression_threshold = self.Ca_threshold * 0.2
        self.depression_duration = random.randint(500,1000) # ms of low calcium required to trigger depression

    def __str__(self):
        return f"Synaps(AMPAR={self.AMPAR}, NMDAR={self.NMDAR}, Vescl_R_prob={self.Vescl_R_prob})"

    def AMPAR_NMDAR_binding_sequence(self, S, cleft_GLu_vol, AMPAR_avalible):
        for t, S in enumerate(Spike_train):
            Spike_History.append(S)
            # Spike triggers glutamate release
            if S == 1 and t > 0 and Spike_train[t - 1] == 0 and Vescl_doc > 0:
                for vesical in range(Vescl_doc):
                    if random.random() <= self.Vescl_R_prob:
                        cleft_GLu_vol += self.GLu_mol_Vescl
                        #K_in =+ 5
                Vescl_doc -= 1
                V_replen_time = t

            if V_replen_time is not None and t - V_replen_time >= random.randint(1,5): #MAIN LEVERAGE POINT FOR DYNAMICS
                Vescl_doc = 1
                V_replen_time = None

            
                if S == 1 and cleft_GLu_vol > 0 and AMPAR_avalible > 0:
                    bind_P = min(0.3, AMPAR_avalible / self.AMPAR)
                    scale_fac = (self.AMPAR / 30)
                    max_binding_events = 50
            if S == 1 and cleft_GLu_vol > 0 and AMPAR_avalible > 0:
                bind_P = min(0.3, AMPAR_avalible / self.AMPAR)
                scale_fac = (self.AMPAR / 30)
                max_binding_events = 50  # limit number of binding attempts per timestep
                binding_attempts = 0
                while cleft_GLu_vol > 0 and AMPAR_avalible > 0 and binding_attempts < max_binding_events:
                    if random.random() < bind_P:
                        self.AMPAR_bound.append((t, AMPAR_avalible))  # record binding event
                        bound_A += 1
                        cleft_GLu_vol -= 1
                        AMPAR_avalible -= 1
                        Na_in += 2.0 * scale_fac
                        Na_out -= 2.0 * scale_fac
                        K_out += 0.05 * scale_fac
                        K_in -= 0.05 * scale_fac
                        Cl_out -= 0.05 * scale_fac
                        Cl_in  += 0.05 * scale_fac
                    binding_attempts += 1
                    if cleft_GLu_vol <= 0 or AMPAR_avalible <= 0:
                        break 
                    
                Vol_AMPAR_bound = len(self.AMPAR_bound)
                
                NMDAR_avalible = self.NMDAR
                if Vol_AMPAR_bound >= self.Ca_threshold:
                    NMDAR_Mg = False
                    
                if NMDAR_Mg == False and bound_N > 0:
                    effective_P_Ca = P_Ca  # Allow Ca²⁺ through NMDARs
                else:
                    effective_P_Ca = 0
                    
                if NMDAR_Mg == False:
                    attempts = 0
                    max_attempts = 50
                    while bound_N < self.NMDAR and cleft_GLu_vol > 0:
                        bind_P_NMDAR = NMDAR_avalible / self.NMDAR
                        for Gl, NMDAR_RE in zip(range(1, int(cleft_GLu_vol) + 1), range(1, NMDAR_avalible + 1)):
                            if random.random() < bind_P_NMDAR:
                                self.NMDAR_bound.append((Gl, NMDAR_RE))
                                bound_N += 1
                                cleft_GLu_vol -= 1
                                NMDAR_avalible -= 1
                                Na_in +=  1.0
                                Na_out -=  1.0
                            bind_P_NMDAR = NMDAR_avalible / self.NMDAR
                        if cleft_GLu_vol <= 0 or NMDAR_avalible <= 0:
                            break
                        attempts += 1
                Vol_NMDAR_bound = len(self.NMDAR_bound)
                Ca_in += Vol_NMDAR_bound
                Ca_out -= Vol_NMDAR_bound# Ca influx from NMDARs
                
    def remove_AMPAR_Ca(self, AMPAR, AMPAR_initial):
    # If calcium stays low long enough, remove AMPARs
        if Syn_Deprs_count >= self.depression_duration and AMPAR > AMPAR_initial:
            AMPAR_ran_rem = random.randrange(200,300)
            AMPAR -= AMPAR_ran_rem 
            AMPAR_avalible = max(AMPAR_avalible - AMPAR_ran_rem , AMPAR_initial)
            Syn_Deprs_count = 0  
            return AMPAR
    
    def insert_AMPAR(self, AMPAR, AMPAR_initial):
    # Plasticity rule: AMPAR insertion if Ca is high
        if self.Ca_in > self.Ca_threshold:
            delta_AMPAR = int((self.Ca_in - self.Ca_threshold) * 4)
            AMPAR += delta_AMPAR
            AMPAR_avalible += delta_AMPAR
        return AMPAR
        
    def remove_AMPAR_volume(self, AMPAR, AMPAR_initial):
        #Remove AMPAR if volum exceeds
        while AMPAR >= AMPAR_initial * 50: #MAIN LEVERAGE POINT FOR DYNAMICS 
            AMPAR -= random.randrange(5,50)
        return AMPAR
    
    def Ca_count(self, Ca_in):
    # Track how long Ca²⁺ has been low
        if Ca_in < self.depression_threshold:    
            Syn_Deprs_count += 1
        else:
            Syn_Deprs_count = 0  # reset if Ca spikes up again

    def Glutamate_appening():
    #Track ion volumes inside and out of the synaps at every step
        Na_in_History.append(Gultamatergic_Synaps.Na_in)
        Na_out_History.append(Gultamatergic_Synaps.Na_out)
        K_in_History.append(Gultamatergic_Synaps.K_in)
        K_out_History.append(Gultamatergic_Synaps.K_out)
        Cl_in_History.append(Gultamatergic_Synaps.Cl_in)
        Cl_out_History.append(Gultamatergic_Synaps.Cl_out)
        Ca_in_History.append(Gultamatergic_Synaps.Ca_in)
        Ca_out_History.append(Gultamatergic_Synaps.Ca_out)
        AMPAR_history.append(Gultamatergic_Synaps.AMPAR) 
        NMDAR_history.append(Gultamatergic_Synaps.NMDAR)
               
##################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################

GABA_Spike_train = np.zeros(9999, dtype=int) #MAIN LEVERAGE POINT FOR DYNAMICS
GABA_Spike_train[500:1000:25] = 1 #MAIN LEVERAGE POINT FOR DYNAMICS
GABA_Spike_train[2000:2500:25] = 1 #MAIN LEVERAGE POINT FOR DYNAMICS
GABA_Spike_train[3500:4000:25] = 1 #MAIN LEVERAGE POINT FOR DYNAMICS
GABA_Spike_train[random.randint(3500,6000)] = 1
GABA_spike_times = [i for i, x in enumerate(Spike_train) if x == 1]
GABA_Spike_History = []

GABA_Ar_history = []
GABA_Br_history = []
GABA_Mem_Vlotage_History = []
GABA_Ca_history = []
GABA_Na_in_History = []
GABA_Na_out_History = []
GABA_K_in_History = []
GABA_K_out_History = []
GABA_Cl_in_History = []
GABA_Cl_out_History = []
GABA_Ca_in_History = []
GABA_Ca_out_History = []

class GABAergic_Synaps:
    def __init__(self, GABA_A, GABA_B, Vescl_R_prob, GABA_Vescl_doc ):
        self.GABA_A = GABA_A
        self.GABA_B = GABA_B    
        self.GABA_Vescl_R_prob = Vescl_R_prob
        #GABAergic variables
        
        self.GABA_K_out = 5.0   # Extracellular K+
        self.GABA_K_in = 127.0  # Intracellular K+

        self.GABA_Na_out = 145.0  # Extracellular Na+
        self.GABA_Na_in = 15.0    # Intracellular Na+

        self.GABA_Cl_out = 110.0  # Extracellular Cl-
        self.GABA_Cl_in = 4.0     # Intracellular Cl-

        self.GABA_Ca_out = 0.1  # Extracellular Ca+2
        self.GABA_Ca_in = 0.0001     # Intracellular Ca+2
        #equilibrium values
        self.GABA_Na_in_eq = 15.0 # equilibrium values for Na
        self.GABA_Na_out_eq = 145.0 # equilibrium values for Na

        self.GABA_K_in_eq = 140.0 # equilibrium values for K
        self.GABA_K_out_eq = 5.0 # equilibrium values for K

        self.GABA_Cl_in_eq = 4.0 # equilibrium values for Cl
        self.GABA_Cl_out_eq = 110.0 # equilibrium values for cl

        self.GABA_Ca_in_eq = 0.001 # equilibrium values for Na
        self.GABA_Ca_out_eq = 100.0 # equilibrium values for Na

        self.GABA_mol_Vescl = 3000 #Glutamate per vesical 
        self.GABA_Vescl_doc = 2 #maximum docking cap of pre synaptic membrain
        self.GABA_Vescl_R_prob = 0.80 #MAIN LEVERAGE POINT FOR DYNAMICS''' # commented out for perameter sweep
        self.GABA_Vescl_R_prob = P_r # for peram sweep
        self.GABA_cleft_vol = 0 # initial level of glutamate in the cleft 
        self.GABA_V_replen_time = None #time it takes for vesical to replenish 

        self.GABA_A_initial = 500 # level of AMPA receptors in post synaptic membrain
        self.GABA_A = self.GABA_A_initial #initial level of AMPA receptors in post synaptic membrain
        self.GABA_A_avalible = GABA_A #level of AMPA receptors in post synaptic membrain that are avalible for GlU to bind 
        self.GABA_A_bound = [] #list of bound AMPARS
        self.bound_GABA_A = 0 #number of AMPARs bound 
        self.GABA_A_insert = random.randrange(200, 300)

        self.GABA_B_initial = 500 # level of AMPA receptors in post synaptic membrain
        self.GABA_B = self.GABA_B_initial #initial level of AMPA receptors in post synaptic membrain
        self.GABA_B_avalible = GABA_B #level of AMPA receptors in post synaptic membrain that are avalible for GlU to bind 
        self.GABA_B_bound = [] #list of bound AMPARS
        self.bound_GABA_B = 0 #number of AMPARs bound 
        self.GABA_B_insert = random.randrange(200, 300)

        self.GABA_Ca_th_coef = 0.8
        self.GABA_Ca_threshold = self.GABA_B_initial * self.GABA_Ca_th_coef #MAIN LEVERAGE POINT FOR DYNAMICS
        #NMDAR * self.Ca_th_coef #MAIN LEVERAGE POINT FOR DYNAMICS
        
        self.depression_threshold = self.GABA_Ca_threshold * 0.2
        self.depression_duration = random.randint(500,1000) # ms of low calcium required to trigger depression

    def GABA_A_GABA_B_binding_sequence(self, S, cleft_GABA_vol, GABA_A_avalible):
        for t, S in enumerate(Mem_Vlotage_History): # GABA_Spike_train):
            GABA_Spike_History.append(S) 
            # Spike triggers glutamate release
            if S >= 55 and t > 0 and GABA_Spike_train[t - 1] == 0 and GABA_Vescl_doc > 0: #
                for vesical in range(GABA_Vescl_doc): # vesical docing
                    if random.random() <= self.GABA_Vescl_R_prob: # vesical release probability
                        self.GABA_cleft_vol += self.GABA_mol_Vescl # vesical release,
                self.GABA_Vescl_doc -= 1
                self.GABA_V_replen_time = t

            if self.GABA_V_replen_time is not None and t - self.GABA_V_replen_time >= random.randint(1,5): # if the time since vesical replenishment is greater than a random time, reset the vesical docing cap
                self.GABA_Vescl_doc = 1 # reset vesical docing cap
                self.GABA_V_replen_time = None # reset vesical replenishment time

            # GABA_A binding
            if S == 1 and self.GABA_cleft_vol > 0 and self.GABA_A_avalible > 0:
                bind_P = min(0.3, self.GABA_A_avalible / self.GABA_A)
                scale_fac = (self.GABA_A / 30)
                max_binding_events = 50  # limit number of binding attempts per timestep
                binding_attempts = 0
                while self.GABA_cleft_vol > 0 and self.GABA_A_avalible > 0 and binding_attempts < max_binding_events:
                    if random.random() < bind_P:
                        self.GABA_A_bound.append((t, self.GABA_A_avalible))  # record binding event
                        self.bound_GABA_A += 1 ## number of GABA_A bound
                        self.GABA_cleft_vol -= 1           # decrease GABA in cleft
                        self.GABA_A_avalible -= 1      # decrease available GABA_A receptors
                        self.GABA_Cl_out -= 0.05 * scale_fac  # Scale Cl- outflow
                        self.GABA_Cl_in  += 0.05 * scale_fac  # Scale Cl- inflow
                    binding_attempts += 1   # increase binding attempts
                    if self.GABA_cleft_vol <= 0 or self.GABA_A_avalible <= 0:     # if cleft GABA volume is 0 or available GABA_A receptors are 0, break
                        break 
                Vol_GABA_A_bound = len(self.GABA_A_bound)    # number of GABA_A bound
                
                if self.GABA_cleft_vol > self.GABA_B_threshold and self.GABA_B_decay_counter == 0: # Check if GABA_B is active
                    self.GABA_B_active = True 
                    self.GABA_B_scaling = 1.0
                    self.GABA_B_decay = 0.99
                    
                if self.GABA_B_active: # If GABA_B is active, apply scaling to K+ and Cl- concentrations
                    self.GABA_K_out += 0.1 * self.GABA_B_scaling # Scale K+ outflow
                    self.GABA_K_in -= 0.1 * self.GABA_B_scaling # Scale K+ inflow
                    self.GABA_B_decay_counter += 1 # Increment decay counter
                    self.GABA_B_scaling *= self.GABA_B_decay
                    if self.GABA_B_scaling < 0.01: # If scaling is very low, deactivate GABA_B
                        self.GABA_B_active = False
                        self.GABA_B_decay_counter = 0
            
    def GABA_B_enable(self, GABA_B_active, GABA_B_scaling):     
        if self.GABA_Ca_in > self.GABA_Ca_threshold: # If calcium is high, activate GABA_B
            self.GABA_LTP_enabled = True # Enable GABA LTP
        elif self.GABA_Ca_in < (self.GABA_Ca_threshold * 0.5):
            self.GABA_LTP_enabled = False
        return self.GABA_LTP_enabled
    
    def GABA_insert(self, AMPAR, AMPAR_initial):
        if self.GABA_LTP_enabled and self.GABA_Ca_in > self.GABA_Ca_threshold: # Track high calcium for GABA LTP
            delta_GABA_A = int((self.GABA_Ca_in - self.GABA_Ca_threshold) * 3)
            self.GABA_A += delta_GABA_A # Increase GABA_A receptors
            self.GABA_A_avalible += delta_GABA_A # Increase available GABA_A receptors
            self.GABA_A = min(GABA_A, self.GABA_A_initial * 50) # Clamp to maximum
            self.GABA_A_avalible = min(self.GABA_A_avalible, self.GABA_A_initial * 50) # Clamp to maximum available GABA_A receptors   
            self.GABA_Vescl_R_prob = min(self.GABA_Vescl_R_prob + 0.01, 1.0) # Increase vesicle release probability

    def GABA_count(self, GABA_Ca_in):
        # Track how long Ca²⁺ has been low
        if self.GABA_Ca_in < self.GABA_Ca_threshold:    
            GABA_LTD_counter += 1
        else:
            GABA_LTD_counter = 0  # reset if Ca spikes up again
        return GABA_LTD_counter
        
    def remove_GABA_A(self, GABA_A, GABA_A_initial):
        if self.GABA_LTD_counter > 300: # If calcium is low for a long time, remove GABA_A receptors
            GABA_rmv = random.randint(10, 50)  # Randomly remove GABA_A receptors
            self.GABA_A -= GABA_rmv # Decrease GABA_A receptors
            self.GABA_A_avalible = max(self.GABA_A_avalible - GABA_rmv, self.GABA_A_initial)   # Decrease available GABA_A receptors
            self.GABA_Vescl_R_prob = max(self.GABA_Vescl_R_prob - 0.01, 0.1)  # Decrease vesicle release probability
            self.GABA_LTD_counter = 0
    def GABA_appending(self,  ):   
        #Track ion volumes inside and out of the synaps at every step
        GABA_Na_in_History.append(self.GABA_Na_in)
        GABA_Na_out_History.append(self.GABA_Na_out)
        GABA_K_in_History.append(self.GABA_K_in)
        GABA_K_out_History.append(self.GABA_K_out)
        GABA_Cl_in_History.append(self.GABA_Cl_in)
        GABA_Cl_out_History.append(self.GABA_Cl_out)
        GABA_Ca_in_History.append(self.GABA_Ca_in)
        GABA_Ca_out_History.append(self.GABA_Ca_out)
        GABA_Ar_history.append(self.GABA_A)
        GABA_Br_history.append(self.GABA_B) 
    
        
    #Ion Pumps to stabelize membrain potential 
def Ion_Pump(Na_in, Na_out, K_in, K_out, Na_in_eq, K_in_eq):
    Na_pump_strength = (Na_in - Na_in_eq) * 0.05  
    K_pump_strength = (K_out - K_in_eq) * 0.05  
    Na_out += Na_pump_strength * 10
    Na_in  -= Na_pump_strength * 10
    K_out  += K_pump_strength * 3
    K_in   -= K_pump_strength * 3
    return Na_in, Na_out, K_in, K_out
    
# Passive ion deay
def Ion_decay(K_in, K_out, Cl_in, Cl_out, Ca_in, Ca_out):
    K_decay = 0.8
    Cl_decay = 0.9
    Ca_decay = 0.9
    K_in -= K_in * (1 - K_decay)
    K_out += K_in * (1 - K_decay)
    Cl_in -= Cl_in * (1 - Cl_decay)
    Cl_out += Cl_in * (1 - Cl_decay)
    Ca_in -= Ca_in * (1 - Ca_decay)
    Ca_out += Ca_in * (1 - Ca_decay)
    return K_in, K_out, Cl_in,Cl_out, Ca_in, Ca_out

def K_return(t, K_in, K_out):
    if t > 0 and Spike_train[t - 1] == 1 and Spike_train[t] == 0:
    #Boost K+ outflow to return toward resting potential
        K_in -= 5.0
        K_out += 5.0
    return K_in, K_out

#Close channels if no stim
def close_channels(t,Na_Chanel, K_Chanel):
    if Spike_train[t] == 0:
        Na_Chanel = False
        K_Chanel = False

#Stabilizing K+ leak toward resting potential ---
def K_leak(K_in, K_out):
    GABA_K_in -= 0.01
    GABA_K_out += 0.01
    return GABA_K_in, GABA_K_out

def Na_leak(Na_in, Na_out):
    GABA_Na_in -= 0.01
    GABA_Na_out += 0.01
    return GABA_Na_in, GABA_Na_out

def Cl_leak(Cl_in, Cl_out):
    GABA_Cl_in -= 0.01
    GABA_Cl_out += 0.01
    return GABA_Cl_in, GABA_Cl_out 
    
def Ca_leak(Ca_in, Ca_out):
    GABA_Ca_in -= 0.01
    GABA_Ca_out += 0.01
    return GABA_Ca_in, GABA_Ca_out

def GABA_leak(GABA_in, GABA_out):
    GABA_GABA_in -= 0.01
    GABA_GABA_out += 0.01
    return GABA_in, GABA_out        

#Glutamate cleared by transporters (at every time step, this is unrealistic)
def GLu_leak(cleft_GLu_vol):
    cleft_GLu_vol *= 0.9
    return cleft_GLu_vol

def GABA_leak(cleft_GABA_vol):
    cleft_GABA_vol *= 0.9
    return cleft_GABA_vol

#Clamp to minimums to avoid log erors
def clamp_values(Na_in, Na_out, K_in, K_out, Cl_in, Cl_out, Ca_in, Ca_out):
    GABA_Na_in = min(max(GABA_Na_in, 0.01), 10000)
    GABA_Na_out = min(max(GABA_Na_out, 0.001), 50)
    GABA_K_in = min(max(GABA_K_in, 0.01), 15)
    GABA_K_out = min(max(GABA_K_out, 0.01), 2)
    GABA_Cl_in = min(max(GABA_Cl_in, 5), 50)
    GABA_Cl_out = min(max(GABA_Cl_out, 0.27), 95)
    GABA_Ca_in = min(max(GABA_Ca_in, 0.0001), 0.1)
    GABA_Ca_out = min(max(GABA_Ca_out, 0.01), 25)
    GABA_A = max(GABA_A, GABAergic_Synaps.GABA_A_initial)
    return GABA_Na_in, GABA_Na_out, GABA_K_in, GABA_K_out, GABA_Cl_in, GABA_Cl_out, GABA_Ca_in, GABA_Ca_out, GABA_A




#Track Receptors every timestep
    
     
    
print(f"t={T} | Na_in={Na_in:.2f}, Na_out={Na_out:.2f}, K_in={K_in:.2f}, K_out={K_out:.2f}, Cl_in={K_in:.2f}, Cl_out={K_out:.2f}, Ca_in={K_in:.2f}, Ca_out={K_out:.2f}, Vm={Mem_Vlotage:.4f}")  

##############################################################################################################################################################################################

np.save(f"AMPAR_history_{P_r:.2f}.npy", AMPAR_history)


data = {'Time (ms)': np.arange(len(Mem_Vlotage_History)),
        'spikes 0 - 1':Spike_History,
        'Membrane Voltage (V)': Mem_Vlotage_History,
        'AMPARs': AMPAR_history,
        'Na+ In': Na_in_History,
        'Na+ Out': Na_out_History,
        'K+ In': K_in_History,
        'K+ Out': K_out_History,
        'Cl- In': Cl_in_History,
        'Cl- Out': Cl_out_History,
        'Ca2+ In': Ca_in_History,
        'Ca2+ Out': Ca_out_History,}
df = pd.DataFrame(data)
df.to_excel(r"C:\Users\david\OneDrive\Documents\synapse_simulation_output.xlsx", index=False)


# commented out for perameter sweep
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(10, 6), sharex=True)

# --- Plot spike train ---
ax1.vlines(spike_times, 1, 2)
ax1.set_ylabel('Spikes')
ax1.set_title('Neuronal Spike Train')

# --- Plot AMPAR over time ---
ax2.plot(range(len(AMPAR_history)), AMPAR_history, linewidth=1)
ax2.set_xlim([0, len(Spike_train)])
ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('AMPARs')
ax2.set_title('AMPAR Volume over Time')
ax2.grid(True)

ax3.plot(range(len(Mem_Vlotage_History)), Mem_Vlotage_History, linewidth=1)
ax3.set_xlabel('Time (ms)')
ax3.set_ylabel('Em (mV)')
ax3.set_title('Dynamic Membrane Potential (Goldman Equation)')
ax3.grid(True)

plt.tight_layout()
plt.show(block=False)

fig, (ax4, ax5, ax6, ax7) = plt.subplots(nrows=4, figsize=(10, 6), sharex=True)

ax4.plot(Na_in_History, label='Na+ in', color='red', linewidth=1)
ax4.plot(Na_out_History, label='Na+ out', color='blue', linewidth=1)
ax4.set_ylabel('[Na⁺] (mM)')
ax4.set_title('Sodium Concentrations')
ax4.legend()
ax4.grid(True)

ax5.plot(K_in_History, label='K+ in', color='red', linewidth=1)
ax5.plot(K_out_History, label='K+ out', color='blue', linewidth=1)
ax5.set_ylabel('[K+] (mM)')
ax5.set_title('Potassium Concentrations')
ax5.legend()
ax5.grid(True)

ax6.plot(Cl_in_History, label='Cl- in', color='red', linewidth=1)
ax6.plot(Cl_out_History, label='Cl- out', color='blue', linewidth=1)
ax6.set_ylabel('[Cl-] (mM)')
ax6.set_title('Chloride Concentrations')
ax6.legend()
ax6.grid(True)

ax7.plot(Ca_in_History, label='Ca2+ in', color='red', linewidth=1)
ax7.plot(Ca_out_History, label='Ca2+ out', color='blue', linewidth=1)
ax7.set_ylabel('[Ca2+] (mM)')
ax7.set_title('Calcium Concentrations')
ax7.legend()
ax7.grid(True)

plt.tight_layout()
plt.show()
