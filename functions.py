from math import e
from numpy import random
import matplotlib.pyplot as pl
from scipy import stats

##Chapter 3##

def mean(numbers):
    return float(sum(numbers)) / len(numbers)
    
def median(numbers):  
    sortedNumbers = sorted(numbers)
    if len(sortedNumbers) %2 == 1:
            return sortedNumbers[(len(sortedNumbers)/2)]
    elif len(sortedNumbers) %2 == 0:
            return (float(sortedNumbers[(len(sortedNumbers)/2)-1])+sortedNumbers[(len(sortedNumbers)/2)])/2

def mode(numbers):
    mode_list=[]
    counter = {}
    for number in numbers:
        if number not in counter:
            counter[number] = 1
        else:
            counter[number] += 1
    maxOccurence=max(counter.values())
    for number in counter.keys():
        if counter[number]==maxOccurence:
            mode_list.append(number)
    return mode_list
        
def percentile(numbers,p):
    if p<=0 or p>=100:
        raise Exception('Error, please enter percentile between 0 and 100')
    sortedNumbers = sorted(numbers)
    n=len(numbers)
    i=float(p)*n/100
    if (i%1==0):
        i=int(i)
        return (float(sortedNumbers[i-1])+sortedNumbers[i])/2
    else:
        i=int(i)+1
        return sortedNumbers[i-1]
      
def quartile(numbers,q):
    if q not in (1,2,3):
        raise Exception ('Error, please enter 1,2 or 3 for the percentile')
    if q==1:
        return percentile(numbers,25)
    if q==2:
        return percentile(numbers,50)
    if q==3:
        return percentile(numbers,75)
                       
def num_range(numbers):
    return float(max(numbers))-min(numbers)  
    
def iqr(numbers):
    return quartile(numbers,3)-quartile(numbers,1)  
    
def var_p(numbers):
    variance_p=0.0
    m=mean(numbers)
    for number in numbers:
        variance_p+=pow((m-number),2)
    return variance_p/len(numbers)
          
def var_s(numbers):
    variance_s=0.0
    m=mean(numbers)
    for number in numbers:
        variance_s+=pow((m-number),2)
    return variance_s/(len(numbers)-1)

def std_p(numbers):
    return pow(var_p(numbers),0.5)

def std_s(numbers):
    return pow(var_s(numbers),0.5)

def cv_p(numbers):
    return (std_p(numbers)*100/mean(numbers))

def cv_s(numbers):
    return (std_s(numbers)/mean(numbers))*100

def skewness(numbers):
    n=float(len(numbers))
    m=mean(numbers)
    s=std_s(numbers)
    sum=0.0
    for number in numbers:
        sum+=pow((number-m)/s,3)
    return sum*n/((n-1)*(n-2))

def z_score(numbers,x):
    return (x-mean(numbers))/std_s(numbers)

def outlier_z(numbers):
    outlier_z_list=[]
    m=mean(numbers)
    s=std_s(numbers)
    for number in numbers:
        z=(number-m)/s
        if z>3 or z<-3:
            outlier_z_list.append(number)
    return outlier_z_list

def outlier_iqr(numbers):
    outlier_iqr_list=[]
    left_th=quartile(numbers,1)-(1.5*float(iqr(numbers)))
    right_th=quartile(numbers,3)+(1.5*float(iqr(numbers)))
    for number in numbers:
        if number<left_th or number>right_th:
            outlier_iqr_list.append(number)
    return outlier_iqr_list

def cov_p(x,y):
    mx=mean(x)
    my=mean(y)
    sumxy=0.0
    for i in range(len(x)):
        sumxy+=(x[i]-mx)*(y[i]-my)
    return sumxy/len(x)

def cov_s(x,y):
    mx=mean(x)
    my=mean(y)
    sumxy=0.0
    for i in range(len(x)):
        sumxy+=(x[i]-mx)*(y[i]-my)
    return sumxy/(len(x)-1)
    
def r_pearson_p(x,y):
    return cov_p(x,y)/(std_p(x)*std_p(y))

def r_pearson_s(x,y):
    return cov_s(x,y)/(std_s(x)*std_s(y))
    
    
##Chapter 4 ##
    
def factorial(n):
    if n==0:
        return 1.0
    else:
        return n*factorial(n-1)
        
def combination(N,n):
    return float(factorial(N))/(factorial(n)*factorial(N-n))

def permutation(N,n):
    return factorial(N)/factorial(N-n)
    
    
    
##Chapter 5##

def binomial_dist(n,p,x,B=True):
    if B==1:
        cum=0.0
        for i in range(x+1):
            cum+=combination(n,i)*pow(p,i)*pow((1-p),(n-i))
        return cum
    if B==0:
        return combination(n,x)*pow(p,x)*pow((1-p),(n-x))
        
def poisson_dist(mu,x,B=True):
    if B==1:
        cum=0.0
        for i in range(0,x+1):
            cum+=pow(mu,i)*pow(e,-mu)/factorial(i)
        return cum
    if B==0:
        return pow(mu,x)*pow(e,-mu)/factorial(x)
            
def hypergeometric_dist(N,r,n,x,B=True):
    if B==1:
        cum=0.0
        for i in range(x+1):
            cum+=combination(r,i)*combination(N-r,n-i)/combination(N,n)
        return cum
    if B==0:
        return combination(r,x)*combination(N-r,n-x)/combination(N,n)


##Chapter 6##

def uniform_dist(a,b,x,B=True):
    if a>=b:
        raise Exception('Error, enter a value a smaller than b')
    if B==1:
        if x<a:
            return 0
        elif x>=a and x<=b:
            return float((x-a))/(b-a)
        elif x>b:
            return 1
    if B==0:
        if x<=b and x>=a:
            return 1/float((b-a))
        else:
            return 0

def exponential_dist(mu,x,B=True):
    if x<0:
        raise Exception('Error, x value must be greater than 0')
    if B==1:
        return 1-pow(e,(float(-x)/mu))
    elif B==0:
        return (1/float(mu))*pow(e,(float(-x)/mu))
        

##Chapter 7##

def test_clt_with_uniform_dist(n,t):
    meanList=[]
    sumofmeans=0.0
    sumofdifference=0.0
    standardDeviation=0.0
    for i in range(t):
        sum=0.0
        for j in range(n):
            sum+=random.uniform(0.0, 1.0) 
        meanList.append((sum/n))
        sumofmeans+=sum/n
    
    expectedValue=sumofmeans/t
    for i in range(t):
        sumofdifference+=pow((meanList[i]-expectedValue),2)
    standardDeviation=pow((sumofdifference/(t-1)),0.5)
    pl.hist(meanList, bins=100)
    return (expectedValue,standardDeviation,pl.show())
    
##Chapter 8##    

def ci_mean(L,a,B=False):
    n=len(L)
    marginofError=0.0
    if B==0:
        t_value=stats.t.ppf(1-(a/2), n-1)
        marginofError=std_s(L)* t_value/ pow(n,0.5)
        return (mean(L)-marginofError,mean(L)+marginofError)
    if B==1:
        z_value=stats.norm.ppf(1-(a/2))
        marginofError=std_p(L)* z_value/ pow(n,0.5)
        return (mean(L)-marginofError,mean(L)+marginofError)
        
def ci_proportion(pbar,n,a):
    z_value=stats.norm.ppf(1-(a/2))
    marginofError=pow((pbar*(1-pbar)/n),0.5)*z_value
    if (pbar-marginofError)<0:
        return (0,pbar+marginofError)
    else:
        return (pbar-marginofError,pbar+marginofError)


##Chapter 9##
def hypo_test_for_mean (xbar,hypoMean,n,a,sd,B,T=0):      
    if B==0:
        t=(xbar-hypoMean)/(float(sd)/pow(n,0.5))
        df=n-1
        if T==-1: 
            p_value=stats.t.cdf(t,df)
        if T==0:
            if t<0:
                p_value=2*stats.t.cdf(t,df)
            else:
                p_value=2*(1-stats.t.cdf(t,df))
        if T==1:
            p_value=(1-stats.t.cdf(t,df))
        if p_value<= a:
            return (t,a,p_value,False)
        elif p_value>a:
            return (t,a,p_value,True) 
    if B==1:
        z=(xbar-hypoMean)/(float(sd)/pow(n,0.5))
        if T==-1: 
            p_value=stats.norm.cdf(z)
        if T==0:
            if z<0:
                p_value=2*stats.norm.cdf(z)
            else:
                p_value=2*(1-stats.norm.cdf(z))
        if T==1:
            p_value=(1-stats.norm.cdf(z))

        if p_value<= a:
            return (z,a,p_value,False)
        elif p_value>a:
            return (z,a,p_value,True)    

def hypo_test_for_proportion(pbar,hypoP,n,a,T):
    sErrPbar=pow((hypoP*(1-hypoP)/n),0.5)
    z=(pbar-hypoP)/sErrPbar
    if T==-1: 
        p_value=stats.norm.cdf(z)
    if T==0:
        if z<0:
            p_value=2*stats.norm.cdf(z)
        else:
            p_value=2*(1-stats.norm.cdf(z))
    if T==1:
        p_value=(1-stats.norm.cdf(z))
    
    if p_value<= a:
        return (z,a,p_value,False)
    elif p_value>a:
        return (z,a,p_value,True)    

def power_in_hypo_test_for_mean(trueMean,hypoMean,n,a,sd,B,T=0):
    if B==0:
        df=n-1
        xbar=(stats.t.ppf(a,df)*(float(sd)/pow(n,0.5)))+hypoMean
        if T==-1: 
            z=(xbar-trueMean)/(float(sd)/pow(n,0.5))
            beta=stats.t.cdf(z,df)
        if T==0:
            z=(xbar-trueMean)/(float(sd)/pow(n,0.5))
            if z<0:
                beta=2*stats.t.cdf(z,df)
            else:
                beta=2*(1-stats.t.cdf(z,df))
        if T==1:
            z=(xbar-trueMean)/(float(sd)/pow(n,0.5))
            beta=(1-stats.t.cdf(z,df)) 
            
        return (1-beta)
        
    if B==1:
        xbar=(stats.norm.ppf(a)*(float(sd)/pow(n,0.5)))+hypoMean
        if T==-1: 
            z=(xbar-trueMean)/(float(sd)/pow(n,0.5))
            beta=stats.norm.cdf(z)
        if T==0:
            z=(xbar-trueMean)/(float(sd)/pow(n,0.5))
            if z<0:
                beta=2*stats.norm.cdf(z)
            else:
                beta=2*(1-stats.norm.cdf(z))
        if T==1:
            z=(xbar-trueMean)/(float(sd)/pow(n,0.5))
            beta=(1-stats.norm.cdf(z)) 
            
        return (1-beta)
        
##Chapter 10##

def ci_for_mean_difference(sampleMeans,sampleSizes,SD,a,B=1):
    n1=float(sampleSizes[0])
    n2=float(sampleSizes[1])
    if B==0:
        s1=float(SD[0])
        s2=float(SD[1])
        x=pow(s1,2)/n1
        y=pow(s2,2)/n2
        df = pow(x+y,2)/((1/(n1-1))*pow(x,2)+ (1/(n2-1)) * pow(y,2))
        t_value=stats.t.ppf(1-(a/2),int(df))
        marginOfError= t_value * pow(x+y, 0.5) 
        return ((sampleMeans[0]-sampleMeans[1]-marginOfError),(sampleMeans[0]-sampleMeans[1]+marginOfError))
    if B==1:
        z_value=stats.norm.ppf(1-(a/2))
        marginOfError= z_value * pow((pow(float(SD[0]),2)/n1) + (pow(float(SD[1]),2)/n2), 0.5)  
        return ((sampleMeans[0]-sampleMeans[1]-marginOfError),(sampleMeans[0]-sampleMeans[1]+marginOfError))
    
 
def hypo_test_for_mean_difference (sampleMeans,sampleSizes,SD,a,B=1,D=0,T=0):
    x=pow(float(SD[0]),2)/sampleSizes[0]
    y=pow(float(SD[1]),2)/sampleSizes[1] 
    n1=float(sampleSizes[0])
    n2=float(sampleSizes[1])    
    if B==0:
        t=(sampleMeans[0]-sampleMeans[1]-D)/ pow(x+y,0.5)
        df=pow(x+y,2)/((1/(n1-1))*pow(x,2)+ (1/(n2-1)) * pow(y,2))
        if T==-1: 
            p_value=stats.t.cdf(t,int(df))
        if T==0:
            if t<0:
                p_value=2*stats.t.cdf(t,int(df))
            else:
                p_value=2*(1-stats.t.cdf(t,int(df)))
        if T==1:
            p_value=(1-stats.t.cdf(t,int(df)))
        
        if p_value<= a:
            return (t,a,p_value,False)
        elif p_value>a:
            return (t,a,p_value,True) 
      
    elif B==1:
        z=(sampleMeans[0]-sampleMeans[1]-D)/ pow(x+y,0.5)
        if T==-1: 
            p_value=stats.norm.cdf(z)
        if T==0:
            if z<0:
                p_value=2*stats.norm.cdf(z)
            else:
                p_value=2*(1-stats.norm.cdf(z))
        if T==1:
            p_value=(1-stats.norm.cdf(z))
    
        if p_value<= a:
            return (z,a,p_value,False)
        elif p_value>a:
            return (z,a,p_value,True)   
 
 
def ci_for_proportion_difference(proportions,sampleSizes,a):
    p[0]=float(proportions[0])
    p[1]=float(proportions[1])
    z_value=stats.norm.ppf(1-(a/2))
    marginOfError=z_value*pow((p[0]*(1-p[0])/sampleSizes[0])+(p[1]*(1-p[1])/sampleSizes[1]),0.5)
    return (p[0]-p[1]-marginOfError,p[0]-p[1]+marginOfError)
    
            
def hypo_test_for_proportion_difference(proportions,sampleSizes,a,D=0,T=0):
    p1=float(proportions[0])
    p2=float(proportions[1])
    n1=float(sampleSizes[0])
    n2=float(sampleSizes[1])
    pbar=(n1*p1+n2*p2)/(n1+n2)
    z=(p1-p2)/pow(pbar*(1-pbar)*(1/n1+1/n2),0.5)
    if T==-1: 
        p_value=stats.norm.cdf(z)
    if T==0:
        if z<0:
            p_value=2*stats.norm.cdf(z)
        else:
            p_value=2*(1-stats.norm.cdf(z))
    if T==1:
        p_value=(1-stats.norm.cdf(z))
    
    if p_value<= a:
        return (z,a,p_value,False)
    elif p_value>a:
        return (z,a,p_value,True) 



##Chapter 11##    

def ci_for_population_var(svar,n,a):
    lowerBound= (n-1)*svar/stats.chi2.ppf((1-a/2),(n-1))
    upperBound= (n-1)*svar/stats.chi2.ppf(a/2,(n-1))
    return (lowerBound,upperBound)

def hypo_test_for_population_var(sVar,hVar,n,a,T): #check again!!!
    chi2=(n-1)*sVar/hVar
    if T==-1:
        p_value=stats.chi2.cdf(chi2,(n-1))
    if T==1:
        p_value=1-stats.chi2.cdf(chi2,(n-1))
    
    if p_value<= a:
        return (chi2,a,p_value,False)
    elif p_value>a:
        return (chi2,a,p_value,True) 

def hypo_test_for_two_population_var(sVar,n,a,T):
    if sVar[0]>=sVar[1]:
        var1=float(sVar[0])
        var2=float(sVar[1])
        n1=n[0]
        n2=n[1]
    else:
        var1=float(sVar[1])
        var2=float(sVar[0])
        n1=n[1]
        n2=n[0]
    F=var1/var2
    
    if T==-1:
        p_value=stats.f.cdf(F,n1-1,n2-1)
    if T==1:
        p_value=1-stats.f.cdf(F,n1-1,n2-1)
        
    if p_value<= a:
        return (F,a,p_value,False)
    elif p_value>a:
        return (F,a,p_value,True) 






