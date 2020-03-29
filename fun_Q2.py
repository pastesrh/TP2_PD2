import numpy as np

############### Question 1

#  pour le put up and out, S max est la barrière

################# Question 2
def CN(K,r,T,vol,n_a,n_t,type,barrier = None):

    delta_t = T/n_t
    delta_a = S_max/n_a
    w = 1.2
    v1 = delta_t/delta_s**2
    v2 = delta_t/delta_s
    r = 0.05
    a = r
    c = -0.05


    if barrier != None:
        S_max = barrier
    else:
        S_max = 200


    grid = np.zeros((n_a,n_t))

    # Boundary a optimiser 
    for i in range(n_t):

        k = T - i

        if type == "american call":
            grid[S_max][k] = 0
            grid[0][k] = S_max - K

        elif type == "american put":
            grid[S_max][k] = K
            grid[0][k] = 0

        elif type == "european call":
            grid[0][k] = 0
            grid[S_max][k] = S_max - K*exp(-r*k*delta_t)
        elif type == "european put":
            grid[0][k] = K*exp(-r*k*delta_t)
            grid[S_max][k] = 0

    grid[n_t,:] = max(range(0,S_max,delta_a)-K,0)

    vect_gk = max(range(0,S_max,delta_a)-K,0)
    for i in range(1,n_t-1,1):
        eps = 1
        gu = max(range(0,S_max,delta_a)-K,0)
        while  eps > 0.00001:
            vect_gk1 = np.zeros(n_a)
            vect_gk1[0] = grid[0][n_t-i]
            vect_gk1[S_max] = grid[S_max][n_t-i]
            gu_1 = np.zeros(n_a+1)
            gu_1[0] = grid[0][n_t-i]
            gu_1[S_max] = grid[S_max][n_t-i] 
            for j in range(1,n_a-1,1):
                b = (n_a-j)*S_max*r
                Ak = v1*a-0.5*v2*b
                Bk = -2*v1*a+delta_t*c
                Ck = v1*a+0.5*v2*b
                alphak_1 = -v1*
                betak_1 = 
                lambdak_1 = 
                h_k = -A*gu_1[j-1]-B*gu_1[j]-C*gu_1[j+1]
                gu_1[j] = gu[j] + (w/beta) 
        
            eps = 
            gu = gu_1

        gi = 
        vect_gk = vect_gk1