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

    a = vol ## a verifier
    c = -r


    if barrier != None:
        S_max = barrier
    else:
        S_max = 200


    grid = np.zeros((n_a+1,n_t+1))

    # Boundary a optimiser 
    for i in range(n_t+1):

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

    for i in range(1,n_t,1):
        vect_gk = grid[n_t-i-1,:]
        eps = 1
        gu = vect_gk 
        while  eps > 0.00001:
            gu_1 = np.zeros(n_a+1)  ## k+1
            gu_1[0] = grid[0][n_t-i]
            gu_1[S_max] = grid[S_max][n_t-i] 
            for j in range(1,n_a,1):
                bk = (n_a-j)*S_max*r     ############ A verifier
                # bk_1 = (n_a-j)*S_max*r
                Ak = v1*a-0.5*v2*b
                Bk = -2*v1*a+delta_t*c
                Ck = v1*a+0.5*v2*b
                alphak_1 = -v1*a+0.5*v2*bk
                betak_1 = 2*v1*a-delta_t*c
                lambdak_1 = -v1*a-0.5*v2*bk
                h_k = -Ak*vect_gk[j-1]-Bk*vect_gk[j]-Ck*vect_gk[j+1]
                if type == "european call" or type == "european put":
                    gu_1[j] = gu[j] + (w/betak_1)*(h_k-alphak_1*gu_1[j-1]-betak_1*gu[j]-lambdak_1*gu[j+1]) 
                else : 
                    gu_1[j] = max(gu[j] + (w/betak_1)*(h_k-alphak_1*gu_1[j-1]-betak_1*gu[j]-lambdak_1*gu[j+1]),(n_a-j)*S_max-K)  
            eps = np.sqrt(np.dot((gu-gu_1)[1:-1],(gu-gu_1)[1:-1]))          ############### A verifier
            gu = gu_1
        grid[n_t-i,:] = gu_1

    # frontier efficiente    
    frontier = 

    return grid, frontier