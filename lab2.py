import random
import math

tasks = 'all'

if tasks == 'ask':
  tasks = input('Enter task number or "all" for all tasks:')

# 1. Pollard’s Rho Algorithm for factorization of long integers
# Function to calculate (base^exponent)%modulus
def modul_pow(base, exp, modul):
    result = 1
    while (exp > 0):
        if (exp & 1):
            result = (result * base) % modul
        exp = exp >> 1
        base = (base * base) % modul
    return result
 
# method to return prime divisor for n
def PollardRho(n):
    if (n == 1):
        return n
    if (n % 2 == 0):
        return 2

    x = (random.randint(0, 2) % (n - 2))
    y = x
    c = (random.randint(0, 1) % (n - 1))
    d = 1
 
    while (d == 1):
        x = (modul_pow(x, 2, n) + c + n)%n
        y = (modul_pow(y, 2, n) + c + n)%n
        y = (modul_pow(y, 2, n) + c + n)%n
        d = math.gcd(abs(x - y), n)
 
        # retry if the algorithm fails to find prime factor
        if (d == n):
            return PollardRho(n) 
    return d
 
if tasks == 'all' or tasks == '1':
 
    n = 189241924801294801294810243
    print("1. Pollard’s Rho Algorithm for factorization of long integers")
    print("One of the divisors for", n , "is ",PollardRho(n))
     
 
# 2. Pollard’s Rho Algorithm for logarithms
#  Extended Euclidean Algorithm
def extended_euclidean(a, b):

    if b == 0:
        return a, 1, 0
    else:
        d, xx, yy = extended_euclidean(b, a % b)
        x = yy
        y = xx - (a / b) * yy
        return d, x, y

# Inverse of a in mod n
def inverse(a, n):
    return extended_euclidean(a, n)[1]

# Pollard Step
def xab(x, a, b, G, H, P, Q):

    subset = x % 3 

    if subset == 0:
        x = x*G % P
        a = (a+1) % Q

    if subset == 1:
        x = x * H % P
        b = (b + 1) % Q

    if subset == 2:
        x = x*x % P
        a = a*2 % Q
        b = b*2 % Q

    return x, a, b


def pollard(G, H, Prime):

    Q = (Prime - 1) / 2  # sub group

    x = G*H
    a = 1
    b = 1

    X = x
    A = int(a)
    B = int(b)

    for i in range(1, Prime):
        # Hedgehog
        x, a, b = xab(x, a, b, G, H, Prime, Q)

        # Rabbit
        X, A, B = xab(X, A, B, G, H, Prime, Q)
        X, A, B = xab(X, A, B, G, H, Prime, Q)

        if x == X:
            break

    nom = int(a-A)
    denom = int(B-b)
    print (nom, denom)

    # Compute the inverse to properly compute the fraction mod q
    res = (inverse(denom, int(Q) * nom) % int(Q))
    if verify(G, H, Prime, res):
        return res
    return int(res + Q)

# Verifies a given set of g, h, p and x
def verify(g, h, p, x):
    return pow(int(g), int(x), p) == h

if tasks == 'all' or tasks == '2':
    print ("\n==============")
    print("2. Pollard’s Rho Algorithm for logarithms")
    g=int(5)
    h=int(22)
    p=int(53)

    print ("g =",g, "h =",h, "p =",p)

    print (h,"=",g,"^x (mod",p,")")
    x = int(pollard(g,h,p))
    print ("\nSolution x=",x)

    print ("Solution:",verify(g, h, p, x))
    print ("Checking h=",pow(int(g), int(x), p))


# 3. Calculation of Euler functions
def gcd(a, b):
 
    if (a == 0):
        return b
    return gcd(b % a, a)
 
# A simple method to evaluate Euler Totient Function
def phi(n):
 
    result = 1
    for i in range(2, n):
        if (gcd(i, n) == 1):   #gcd is a Greatest Common Divisor
            result+=1
    return result
 
if tasks == 'all' or tasks == '3.1':
    print ("\n==============")
    print("3.1 Calculation of Euler functions")
    for n in range(1, 11):
        print("phi(",n,") = ",
               phi(n), sep = "")

# 3.Calculation of Möbius functions 

# n is prime or not
def isPrime(n) :
 
    if (n < 2) :
        return False
    for i in range(2, n + 1) :
        if (i * i <= n and n % i == 0) :
            return False
    return True
 
def mobius(N) :

    if (N == 1) :
        return 1
 
# For a prime factor check if i^2 is also a factor
    p = 0
    for i in range(1, N + 1) :
        if (N % i == 0 and
                isPrime(i)) :
 
            # Check if N is divisible by i^2
            if (N % (i * i) == 0) :
                return 0
            else :
                p = p + 1

    if(p % 2 != 0) :
        return -1
    else :
        return 1
 
if tasks == 'all' or tasks == '3.2':
    print ("\n==============")
    print("3.2 Calculation of Möbius functions")
    N = 6
    print ("Mobius defs M(N) at N = {} is: {}" .
         format(N, mobius(N)),end = "\n")
    print ("Mobius defs M(N) at N = {} is: {}" .
        format(30, mobius(30)),end = "\n")
    print ("Mobius defs M(N) at N = {} is: {}" .
        format(210, mobius(210)),end = "\n")


#4. Computation of Legendre symbols
def Legendre(arr, p):
    e = (p - 1) // 2
    results = [pow(a, e, p) for a in arr]
    return [(r-p if r > 1 else r) for r in results]

if tasks == 'all' or tasks == '4.1':
    print ("\n==============")
    print("4.1 Calculation of Legendre symbols")
    p = 19
    arr = [64, 65, 66, 67, 68, 69, 70, 71, 72, 73]
    results = Legendre(arr, p)
    print(arr)
    print(results)


    
#4. Computation of Jacobi symbols
def Jacobi(a, n):
        assert(n > a > 0 and n%2 == 1)
        t = 1
        while a != 0:
            while a % 2 == 0:
                a /= 2
                r = n % 8
                if r == 3 or r == 5:
                    t = -t
            a, n = n, a
            if a % 4 == n % 4 == 3:
                t = -t
            a %= n
        if n == 1:
            return t
        else:
            return 0

if tasks == 'all' or tasks == '4.2':
    print ("\n==============")
    print("4.2 Calculation of Jacobi symbols")
    N = 65
    for _ in range(10):
        a = random.randrange(1, N)
        n = random.randrange(a+1, 2*N)
        if n%2 == 0:
            n += 1
            
        result = Jacobi(a, n)
        print("a = ",a,"n = ",n, "Jacobi symbols = ",result)


#5. Cipolla's algorithm for finding a discrete square root.

#Converts n to base b as a list of integers between 0 and b-1
def convertToBase(n, b):
	if(n < 2):
		return [n];
	temp = n;
	answer = [];
	while(temp != 0):
		answer = [temp % b]+ answer;
		temp /= b;
	return answer;

#Takes integer n and odd prime p and returns both square roots of n module p as a pair (a,b)
def cipolla(n,p):
	n %= p
	if(n == 0 or n == 1):
		return (n,-n%p)
	phi = p - 1
	if(pow(int(n), int(phi/2), int(p)) != 1):
		return ()
	if(p%4 == 3):
		answer = pow(int(n),int((p+1)/4),int(p))
		return (answer,-answer%p)
	aa = 0
	for i in range(1,p):
		temp = pow(int((i*i-n)%p),int(phi/2),int(p))
		if(temp == phi):
			aa = i
			break;
	exponent = convertToBase((p+1)/2,2)
	def cipollaMult(a,b,c,d,w,p):
		return ((a*c+b*d*w)%p,(a*d+b*c)%p)
	x1 = (aa,1)
	x2 = cipollaMult(x1,x1,aa*aa-n,p)
	for i in range(1,len(exponent)):
		if(exponent[i] == 0):
			x2 = cipollaMult(x2,x1,aa*aa-n,p)
			x1 = cipollaMult(x1,x1,aa*aa-n,p)
		else:
			x1 = cipollaMult(x1,x2,aa*aa-n,p)
			x2 = cipollaMult(x2,x2,aa*aa-n,p)
	return (x1[0],-x1[0]%p)

if tasks == 'all' or tasks == '5':
    print ("\n==============")
    print("5. Cipolla's algorithm")
    print("Roots of 2 mod 7: " +str(cipolla(2,7)))
    print("Roots of 8218 mod 10007: " +str(cipolla(8218,10007)))


#6 Solovey Strassen Algorithm for checking numbers for simplicity

def solovay_strassen(n, k=10):  
     if n == 2 or n == 3:  
         return True  
     if not n & 1:  
         return False  
           
     for i in range(k):  
         a = random.randrange(2, n - 1)       # choose any random number from 1 to (n-1)  
         x = Jacobi(a, n)                     # find n's jacobi number  
           
         y = pow(int(a), int((n - 1) / 2), int(n))           # calculate legendre symbol from euler criterion formula  
         if y != 1 and y != 0:                          
             y = -1  
         
         if (x == 0) or (y != x):             # if jecobi and eular criterion formula are not same (y != x) then the number is not prime  
             return False         
     return True  
   
if tasks == 'all' or tasks == '6':
    print ("\n==============")
    print("6. Strassen Algorithm")      
    n = 13
    if(n<=1):                                          
         print ('Number must be >=2')
    elif(solovay_strassen(n,10)):  
         print (n,'is Prime!' )
    else:  
         print (n,'is not Prime!')
   

#7 Elgamal Encryption using Elliptic Curve
def polynomial(LHS,RHS,n):
    for i in range(0,n):
        LHS[0].append(i)
        RHS[0].append(i)
        LHS[1].append((i*i*i + a*i + b)%n)
        RHS[1].append((i*i)%n)


def points_generate(arr_x,arr_y,n):
    count=0
    for i in range(0,n):
        for j in range(0,n):
            if(LHS[1][i]==RHS[1][j]):
                count+=1
                arr_x.append(LHS[0][i])
                arr_y.append(RHS[0][j])
    return count

if (tasks == 'all' or tasks == '7'):
  print ("\n==============")
  print("7. Elgamal Encryption using Elliptic Curve")
  LHS=[[]]
  RHS=[[]]
  LHS.append([])
  RHS.append([])
  a=15
  b=13
  n = 6
  print("n = ",n)
  print("a = ",a, " b = ",b)

  #Polynomial
  polynomial(LHS,RHS,n)

  arr_x=[]
  arr_y=[]
  #Generating base points
  count=points_generate(arr_x,arr_y,n)
      
  #Print Generated Points
  print("Generated points are:")
  for i in range(0,count):
      print(i+1," (",arr_x[i],",",arr_y[i],")")

  #Calculation of Base Point
  bx=arr_x[0]
  by=arr_y[0]
  print("Base Point taken is:(",bx,",",by,")")

  d=3
  print("d =",d)
  if(d>=n):
      print("'d' should be less than 'n'.")
  else:
      #Q i.e. sender's public key generation
      Qx=d*bx
      Qy=d*by
      print("Public key of sender is:(",Qx,",",Qy,")")

      #Encrytion
      k = 4
      print("k = ",k)
      if(k>=n):
          print("'k' should be less than 'n'")
      else:
          M=2
          print("The message to be sent:", M)

          #Cipher text 1 generation
          C1x=k*bx
          C1y=k*by
          print("Value of Cipher text 1 i.e. C1:(",C1x,",",C1y,")\n")

          #Cipher text 2 generation
          C2x=k*Qx+M
          C2y=k*Qy+M
          print("Value of Cipher text 2 i.e. C2:(",C2x,",",C2y,")\n")

      #Decryption
          Mx=C2x-d*C1x
          My=C2y-d*C1y
          print("The message recieved by reciever is:",Mx)      
