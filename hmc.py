from __future__ import division

import matplotlib.pyplot as plt
import numpy as np

N_1 = 1
N_2 = 1
N = N_1 + N_2
k = 1
delta_T = 0.5
total_steps = 750


def calculate_force(x_s, gaussian_weight=0.5):
    forces = [0.0] * N

    for i in range(N):
        forces[i] = 0. - ( forces[i] + 0.5 * np.tanh(0.5 * x_s[i]) )

    return forces

def calculate_hamiltonian(x_s, p_s, gaussian_weight=0.5):
    # The Hamiltonian is the sum of the kinetic term, the potential term and each of the individual interaction terms.
    

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
            potential -= np.log(abs(np.tanh(gaussian_weight * (x_s[i] - x_s[j]))))

    for i in range(N_1 + 1, N):
        for j in range(i + 1, N):
            potential -= np.log(abs(np.tanh(gaussian_weight * (x_s[i] - x_s[j]))))

    potential *= 2.
    
    action = action + potential

    hamiltonian = kinetic_contribution + action

    return hamiltonian
    
def hmc(n_step=1, x_s={}, p_s={}):

    if n_step >= total_steps:
        return 

    # print "Step #{}".format(n_step)

    # First, choose a position freely.
    if n_step == 1:
        x_s = {}
        p_s = {}


        # Generate initial configuration.

        # Generate positions first.

        for i in xrange(0, N - 1, 2):
            # Generate a random number in the range- [0.5, 0.5)
            x_s[i] = np.random.random() - 0.5  
            x_s[i + 1] = np.random.random() - 0.5


        # Next, generate canonical momentum.
        # Sample this from a Gaussian distribution: e^( - p^2 / 2)
        mu, sigma = 0, 1 / np.sqrt(2) # mean and standard deviation

        for i in xrange(0, N - 1, 2):
            # Scale by sqrt(pi) to make the coefficient 1.
            p_s[i] = np.sqrt(np.pi) * np.random.normal(loc=0.0, scale=sigma)
            p_s[i + 1] = np.sqrt(np.pi) * np.random.normal(loc=0.0, scale=sigma)        
    
    # Next, using the molecular dynamics, to find the state in t = n + 1.

    # The molecular dynamics comes from Hamilton-Jacobi equations.
    # dx/dt = p and dp/dt = - \frac{ \partial S(x) }{ \partial x }

    # So we need to find the force.
    forces = calculate_force(x_s)

    # Now, use the leapfrog method to propagate the system.
    p_2_s = {}              # p_2_s = p(t + delta_T / 2)
    for i in range(N):
        p_2_s[i] = p_s[i] - (delta_T / 2) * x_s[i]

    x_3_s = {}              # x_3_s = x(t + delta_T)
    for i in range(N):
        x_3_s[i] = x_s[i] + delta_T * p_2_s[i]

    p_3_s = {}
    for i in range(N):
        p_3_s[i] = p_2_s[i] - (delta_T / 2) * x_3_s[i]

    hamiltonian = calculate_hamiltonian(x_s, p_s)
    new_hamiltonian = calculate_hamiltonian(x_3_s, p_3_s)

    delta_hamiltonian = new_hamiltonian - hamiltonian

    
    if float(delta_hamiltonian) < float(0.):
        transition_probability = 1.
    else:
        transition_probability = np.exp( - delta_hamiltonian)


    # Generate a new random number between 0 and 1 and if it is less than the transition probability, accept the new state.
    some_random_number = np.random.rand()

    if some_random_number <= transition_probability:
        # print "Transition to new state"
        new_x_s = x_3_s
        new_p_s = p_3_s
    else:
        # print "Stick with the current state"
        new_x_s = x_s
        new_p_s = p_s

    # Calculate the raito of free energy.
    ratio_value = 0.
    for i in range(N_1):
        for j in range(N_1, N):
            ratio_value = np.tanh(new_x_s[i] - new_x_s[j]) ** 2
            print ratio_value

    hmc(n_step=(n_step + 1), x_s=new_x_s, p_s=new_p_s)



        


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