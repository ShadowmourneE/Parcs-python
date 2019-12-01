# import gmpy2
from Pyro4 import expose
import random

def det(l):
        n = len(l)
        AM = [row[:] for row in l]
 
        for fd in range(n): 
            for i in range(fd+1,n): 
                if AM[fd][fd] == 0: 
                    AM[fd][fd] == 1.0e-18
                crScaler = AM[i][fd] / (AM[fd][fd] + 1.0e-9)
                for j in range(n): 
                    
                    AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
    
        product = 1.0
        for i in range(n):
            product *= AM[i][i]
            
        return product

def getMatrixMinor(m,i,j):
        return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

class Solver:
    def __init__(self, workers=None, input_file_name=None, output_file_name=None):
        self.input_file_name = input_file_name
        self.output_file_name = output_file_name
        self.workers = workers
        print("Inited")

    def solve(self):
        print("Job Started")
        print("Workers %d" % len(self.workers))

        l = self.read_input()
        
        n = len(l)
        detA = det(l)

        parts = []
        for k in range(len(self.workers)):
            parts.append([])
        
        for i in range(len(l)):
            for j in range(len(self.workers)):
                if i%len(self.workers) == j:
                    parts[j].append(l[i])

        # map
        mapped = []
        sign = 1    
        for i in range(len(self.workers)):
            print("map %d" % i)
            mapped.append(self.workers[i].mymap(parts[i], l, sign,i).value)
            sign*=-1

        # reduce
        
        work = len(self.workers)
        
        minors = self.myreduce(mapped,work)
        self.write_output(minors)
        A1 = 1 / detA
        

        for x in range(len(minors)):
            for y in range(len(minors[x])):
                minors[x][y] = minors[x][y]*A1

        # output
        self.write_output(minors)

        print("Job Finished")


    @staticmethod
    @expose
    def mymap(part, A, sign, time):
        M = []
        k=time
        for i in range(len(part)):
            M.append([])

        for i in range(len(part)):
            
            for j in range(len(part[i])):
                
                minor = getMatrixMinor(A,k,j)
                M[i].append(det(minor)*sign)
                sign*=-1
            k+=3
            if(k>10):
                k=0    
        return M
        

    @staticmethod
    @expose
    def myreduce(partAll, workNum):
        res = []   
        row=[]
        k=0
        p=0

        i=0
        
        while i < len(partAll):
            
            try:
                row.append(partAll[i][p][k])
                
                if i == workNum-1:
                    i=-1
                    p+=1
            except IndexError:
                k+=1
                i=-1
                p=0
                res.append(row)
                print("reduce loop")
                print(res)
                row = []
            if k == 9:
                print("reduce done")
                return res
            i+=1

    def read_input(self):
        with open(self.input_file_name, 'r') as f:
            l = [[int(num) for num in line.split(",")] for line in f]
        f.close()
        return l

    def write_output(self, output):
        with open(self.output_file_name, 'w') as f:
            for item in output:
                f.write("%s\n" % item)      
        f.close()
        print("output done")

    @staticmethod
    @expose
    def is_probable_prime(n):
        """
        Miller-Rabin primality test.
        A return value of False means n is certainly not prime. A return value of
        True means n is very likely a prime.
        >>> is_probable_prime(1)
        Traceback (most recent call last):
            ...
        AssertionError
        >>> is_probable_prime(2)
        True
        >>> is_probable_prime(3)
        True
        >>> is_probable_prime(4)
        False
        >>> is_probable_prime(5)
        True
        >>> is_probable_prime(123456789)
        False
        >>> primes_under_1000 = [i for i in range(2, 1000) if is_probable_prime(i)]
        >>> len(primes_under_1000)
        168
        >>> primes_under_1000[-10:]
        [937, 941, 947, 953, 967, 971, 977, 983, 991, 997]
        >>> is_probable_prime(6438080068035544392301298549614926991513861075340134\
    3291807343952413826484237063006136971539473913409092293733259038472039\
    7133335969549256322620979036686633213903952966175107096769180017646161\
    851573147596390153)
        True
        >>> is_probable_prime(7438080068035544392301298549614926991513861075340134\
    3291807343952413826484237063006136971539473913409092293733259038472039\
    7133335969549256322620979036686633213903952966175107096769180017646161\
    851573147596390153)
        False
        """
        assert n >= 2
        # special case 2
        if n == 2:
            return True
        # ensure n is odd
        if n % 2 == 0:
            return False
        # write n-1 as 2**s * d
        # repeatedly try to divide n-1 by 2
        s = 0
        d = n - 1
        while True:
            quotient, remainder = divmod(d, 2)
            if remainder == 1:
                break
            s += 1
            d = quotient
        assert (2 ** s * d == n - 1)

        # test the base a to see whether it is a witness for the compositeness of n
        def try_composite(a):
            if pow(a, d, n) == 1:
                return False
            for i in range(s):
                if pow(a, 2 ** i * d, n) == n - 1:
                    return False
            return True  # n is definitely composite

        for i in range(1):
            a = random.randrange(2, n)
            if try_composite(a):
                return False

        return True  # no base tested showed n as composite