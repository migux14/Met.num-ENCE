#ESCOLA NACIONA DE CI�NCIAS ESTAT�STICAS - ENCE
#Autor: Miguel do Nascimento Faria Conforto

#"C�lculo Num�rico: Aspectos T�oricos e Computacionais" - Ruggiero
#ALGORITMO 1
#Resolu��o de sistemas triangulares superiores
resolucao_sistema_triangular_superior <- function(A, b){
  
  X <- rep(0, length(b)) #Criando vetor coluna X do sistema AX = b
  n <- length(X)
  
  X[n] = b[n]/A[n, n] #Valor Xn 
  
  for (k in (n - 1):1 ){
    s <- 0
    for (j in (k+1):n )( 
      s = s + (A[k, j] * X[j]) 
  )
    
    X[k] = (b[k] - s)/ A[k, k]  #Valor de Xk
  }
  
  return(X)
}

#ALGORITMO 2
#Elimina��o de Gauss em Ax = b (A: nxn, x: nx1, b: nx1)
eliminacao_de_gauss <- function(A, b) {
  
  n <- length(b) #Tanto A quanto b e X t�m o mesmo tamanho
  
  for (k in 1:(n - 1)) {
    for (i in (k + 1):n) {
      
      m = A[i, k] / A[k, k]
      A[i, k] = 0
      
      b[i] = b[i] - (m * b[k])
      
      for (j in (k + 1):n ) {
        A[i, j] = A[i, j] - (m * A[k, j])
    }
  }
}
  return(list("A" = A, "b" = b))
}

#Resolu��o de sistemas lineares atrav�s da Elimina��o de Gauss

resolucao_sistemas_gerais <- function(A, b) {
  
  n <- length(b) #Tanto A quanto b e X t�m o mesmo tamanho
  
  for (k in 1:(n - 1)) {
    for (i in (k + 1):n) {
      
      m = A[i, k] / A[k, k]
      A[i, k] = 0
      
      b[i] = b[i] - (m * b[k])
      
      for (j in (k + 1):n ) {
        A[i, j] = A[i, j] - (m * A[k, j])
      }
    }
  }
  return(resolucao_sistema_triangular_superior(A, b))
}

