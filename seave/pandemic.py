import numpy as np
import datetime
import scipy.stats as spy
import matplotlib.pyplot as plt

class Pandemic:    
    def get_infections_function(self):
        np = len(self.__pinfections)
        ps = self.__pinfections
        return lambda x: sum([p*pow(x,np - 1 - i)
                               for i,p in list(enumerate(ps))
                              ])

    def plot(self,start=None,ndays=1000):
        if start == None:
            start = self.start_date
        x = [start+datetime.timedelta(days=i) for i in range(0,ndays)]
        y = [
            self.get_infection_rate(i)
            for i in x ]
        plt.plot(x, y,'--')
        plt.xticks(rotation=90)
        plt.show();               

    
    def days_to_x(self,d):
        return d/30.
    
    def get_infection_rate(self,date):
        days = (date-self.start_date).days
        x = self.days_to_x(days)
        retval = self.__pinfections_fn(x)
        retval = retval if retval > 0 else 0
        return 5*retval**2 * 0.0001
    
    def create_variant_fn(self,start,end):
        width = (end-start).days
        fn = spy.norm(loc=0.5*width,scale=width*0.25)
        
        def _fn(date):
            x = (date - start).days
            return fn.pdf(x)

        return lambda x : _fn(x)
        
    def get_variants(self,date):
        variants = {
            k:v(date)
            for k,v in self.variants.items()
        }
        total = sum(list(variants.values()))
        return {
            k:v/total
            for k,v in variants.items()
        }
    
    def __init__(self,start,variants):
        self.start_date = start
        self.__pinfections = np.array([-1.24000125e-09,  1.86510982e-07, -1.19526551e-05,  4.25541261e-04,
                                       -9.20077170e-03,  1.24057027e-01, -1.03106570e+00,  5.03315072e+00,
                                       -1.30288337e+01,  1.42349019e+01, -2.29713012e-01])
        self.__pinfections_fn = self.get_infections_function()
        
        self.variants = {
            k:self.create_variant_fn(start,end)
            for k,(start,end,_) in variants.items()
        }
        self.variant_severity = {
            k:s
            for k,(_,_,s) in variants.items()
        }
