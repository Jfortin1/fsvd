fsvd <- function(A, k, i = 1, p = 2, method="svd"){
    l <- k + p 
    n <- ncol(A)
    m <- nrow(A)
    tall <- TRUE
    type="tall"

    # In the case a WIDE matrix is provided:
    if (m < n){
        A <- t(A)
        n <- ncol(A)
        m <- nrow(A)
        tall <- FALSE
        type="wide"
    }
    if (l > ncol(A) && method!="exact"){
        stop("(k+p) is greater than the number of columns. Please decrease the value of k.")
    }

    # Construct G to be n x l and Gaussian
    G <- matrix(rnorm(n*l,0,1), nrow=n, ncol=l)


    # SVD approach
    if (method == "svd"){
        # Power method: 
        H <- A %*% G # m x l matrix
        for (j in 1:i){
            H <- A %*% (t(A) %*% H)
        }

        # We use a SVD to find an othogonal basis Q:
        # H = F %*% Omega %*% t(S)
        svd <- svd(t(H) %*% H)
        F   <- svd$u # l x l
        omega <- diag(1/sqrt(svd$d)) # l x l
        S <- H %*% F %*% omega # m x l 

        # Define the orthogonal basis:
        Q <- S[,1:k] # m x k

        T <- t(A) %*% Q # n x k 
        T <- t(T)
    } 


    # QR approach
    if (method == "qr"){
        # Need to create a list of H matrices
        h.list <- vector("list", i+1)
        h.list[[1]] <- A %*% G 
        for (j in 2:(i+1)){
            h.list[[j]] <- A %*% (t(A) %*% h.list[[j-1]])
        }
        H <- do.call("cbind",h.list) # n x [(1+1)l] matrix

        # QR algorithm
        Q <- qr.Q(qr(H,0))
        T = t(A)%*%Q # n x [(i+1)l]
        T <- t(T) 
    }

    if (method == "svd" | method == "qr"){
        svd = svd(T)
        u <- Q %*% svd$u[,1:k]
        v <- svd$v[,1:k]
        d <- svd$d[1:k]

    # For exact SVD:
    } else {
        svd <- svd(A)
        u <- svd$u[,1:k]
        v <- svd$v[,1:k]
        d <- svd$d[1:k]
    }


    if (!tall){
        uu <- v
        v <- u
        u <- uu
    } 
    list(u=u,v=v,d=d)
}

