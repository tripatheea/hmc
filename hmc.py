from __future__ import division

import matplotlib.pyplot as plt
import numpy as np
import sys


N_1 = 1
N_2 = 1
N = N_1 + N_2
k = 1
delta_tau = 0.5
n_tau = 2
sigma = 2.


# Change this.
n_iterations = 1000000

debug = False

def calculate_force(x_s, gaussian_weight=0.5):
    forces = [0.0] * N

    for i in range(N):
        forces[i] = - ( forces[i] + 0.5 * np.tanh(0.5 * x_s[i]) )

    return forces

def calculate_hamiltonian(x_s, p_s, gaussian_weight=0.5):
    # The Hamiltonian is the sum of the kinetic term, the potential term and each of the individual interaction terms.
    
    reg_vander = 0.00000001
    # Compute the kinetic term first: 1/2 p^2

    kinetic_contribution = 0.
    for i in range(N):
        kinetic_contribution += p_s[i] * p_s[i]
    kinetic_contribution *= 0.5

    # Next, compute the action.
    action = 0.
    for i in range(N):
        action += np.log(np.cosh(0.5 * x_s[i]))

    # Now contribution from the interaction terms.
    potential = 0.
    for i in range(N_1):
        for j in range(i + 1, N_1):
            potential -= np.log(abs(np.tanh(gaussian_weight * (x_s[i] - x_s[j]))) + reg_vander )

    for i in range(N_1 + 1, N):
        for j in range(i + 1, N):
            potential -= np.log(abs(np.tanh(gaussian_weight * (x_s[i] - x_s[j]))) + reg_vander)

    potential = potential * 2.
    
    action = action + potential

    hamiltonian = kinetic_contribution + action

    return hamiltonian


def tests():
    x_s = [0.5, 0.6]
    forces = calculate_force(x_s, gaussian_weight=1.)

    print "The force for x = ", x_s[0], " is ", forces


def hmc(tau=1, x_s={}, p_s={}):

    x_s = {}
    p_s = {}


    for i in xrange(0, N - 1, 2):
        rand_1 = np.random.random()
        rand_2 = np.random.random()

        x_s[i] = np.sqrt( - 2. * sigma * np.log(rand_1)) * np.cos( np.pi * rand_2)
        x_s[i + 1] = np.sqrt( - 2. * sigma * np.log(rand_1)) * np.sin( np.pi * rand_2)

        # x_s[i] = 0.5  
        # x_s[i + 1] = 0.6

    rand_1 = np.random.random()
    rand_2 = np.random.random()
    x_s[N - 1] = np.sqrt( - 2. * sigma * np.log(rand_1) ) * np.cos(np.pi * rand_2)


    total_value = 0.

    all_x_s, all_p_s = [], []    
    ratio_values = []
    for iter in range(1, n_iterations + 1):

        if debug:
            print "\n"

        for i in xrange(0, N - 1, 2):
            
            rand_1 = np.random.random()
            rand_2 = np.random.random()

            p_s[i] = np.sqrt( - 2. * np.log(rand_1)) * np.cos( np.pi * rand_2)
            p_s[i + 1] = np.sqrt( - 2. * np.log(rand_1)) * np.sin( np.pi * rand_2)

            # p_s[i] = 0.1
            # p_s[i + 1] = 0.2
        
        rand_1 = np.random.random()
        rand_2 = np.random.random()
        p_s[N - 1] = np.sqrt( - 2. * np.log(rand_1) ) * np.cos(np.pi * rand_2)

        forces = calculate_force(x_s)

        # print "\nThe forces are ", forces[0], forces[1]
        
        
        p_2_s = {}
        for i in range(N):
            p_2_s[i] = p_s[i] + (delta_tau / 2) * forces[i]
        
        x_2_s = {}
        for i in range(N):
            x_2_s[i] = x_s[i] + delta_tau * p_2_s[i]

        for tau in range(1, n_tau):
            
            forces = calculate_force(x_2_s)
            
            for i in range(N):
                p_2_s[i] = p_2_s[i] + delta_tau * forces[i]
            
            for i in range(N):
                x_2_s[i] = x_2_s[i] + delta_tau * p_2_s[i]

        forces = calculate_force(x_2_s)

        for i in range(N):
            p_2_s[i] = p_2_s[i] + (delta_tau / 2) * forces[i]
        
        if debug:
            print "The new x_s are ", x_2_s[0], x_2_s[1]
            print "The new p_s are ", p_2_s[0], p_2_s[1]


        hamiltonian = calculate_hamiltonian(x_s, p_s)
        new_hamiltonian = calculate_hamiltonian(x_2_s, p_2_s)

        delta_hamiltonian = new_hamiltonian - hamiltonian

        if debug:
            print "Old Hamiltonian = ", hamiltonian
            print "New Hamiltonian = ", new_hamiltonian
        # print "delta Hamiltonian = ", delta_hamiltonian

        if float(delta_hamiltonian) < float(0.):
            transition_probability = 1.
        else:
            transition_probability = np.exp( - delta_hamiltonian)


       
        some_random_number = np.random.rand()
        # some_random_number = 0.5



        if some_random_number <= transition_probability:
            # print "Transition to new state"
            x_s = x_2_s
            # p_s = p_2_s
        else:
            # print "Stick with the current state"
            # new_x_s = x_s
            # new_p_s = p_s
            pass

        # Calculate the ratio of free energy.
        ratio_value = 1.
        for i in range(N_1):
            # print "Hey i=", i
            for j in range(N_1, N ):
                # print "Hello j=", j
                
                # print "i, j", x_s[i], x_s[j]

                ratio_value = ratio_value * np.tanh(0.5 * (x_s[i] - x_s[j]))
        
        ratio_value = ratio_value * ratio_value

        total_value += ratio_value

        ratio_values.append( ratio_value )

        # hmc(tau=(tau + 1), x_s=new_x_s, p_s=new_p_s)
        
        all_x_s.append( x_s )
        all_p_s.append( p_s )
        


    # print "ratio values", ratio_values
    for i in range(len(ratio_values)):
        # print ratio_values[i] 
        
        pass

    # plt.plot([i for i in range(len(ratio_values))], ratio_values)

    # plt.show()

    print n_iterations, total_value / n_iterations


        


def weird_distribution():
    lambda_1 = []
    lambda_2 = []
    for i in range(100000):
        rand1 = np.random.random()
        rand2 = np.random.random()

        lambda_1.append(np.sqrt( - 2 * 1 * np.log(rand1) ) * np.cos(np.pi * rand2))
        lambda_2.append(np.sqrt( - 2 * 1 * np.log(rand1) ) * np.sin(np.pi * rand2))

    # print lambda_1
    plt.hist(lambda_1, bins=50, histtype='step', label="eigen 1", lw=3)
    plt.hist(lambda_2, bins=50, histtype='step', label="eigen 2", lw=3)

    plt.legend()

    plt.show()
    

if __name__ == '__main__':
    hmc()

    # tests()