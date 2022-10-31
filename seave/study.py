import scipy.stats as spy
import numpy as np
import datetime
import uuid
m_igg_scale = {
    'AZ':0.5,
    'Pfizer':1,
    'Moderna':1.2,
    'Natural':0.01,
}


def get_p(x,threshold,a=2*np.e):
    x = x/threshold
    retval = (np.power(a,x) - 1)/(a-1)
    return retval if retval < 1 else 1

symptoms = {
    'headache':lambda x: np.random.uniform() < get_p(x,threshold=200),
    'fatigue':lambda x: np.random.uniform() < get_p(x,threshold=400),
    'short breathe':lambda x: np.random.uniform() < get_p(x,threshold=600),
    'cant breathe':lambda x: np.random.uniform() < get_p(x,threshold=100000),
    'death':lambda x: np.random.uniform() < get_p(x,threshold=500000),
}

def get_symptoms(x):
    return [
        k
        for k,fn in symptoms.items()
        if fn(x)
    ]



class Study:
    
    @staticmethod 
    def age_to_dob(age):
        dob = datetime.date.today() - datetime.timedelta(days=365.25*age)
        return dob

    def get_demo(self):
        return {
            'id':self.id,
            'age':self.age,
            'sex':self.sex,
            'bmi':self.bmi,
            'simd':self.simd,
            'location':self.location,
            'comorbidities':[c.get_name() for c in self.comorbidities]
        }
    def __init__(self,
                 dob=datetime.date(1960,1,1),
                 sex='F',
                 bmi=20,
                 simd=3,
                 location=1,
                 comorbidities=[]
                ):
        self.id = uuid.uuid4()
        self.dob = dob
        self.age = None
        self.age = self.get_age()
        self.sex = sex
        
        self.date_of_vaccines = []
        self.vaccine_products = []
        
        self.date_of_infections = []
        self.infection_variants = []

        self.outcomes = []
        self.comorbidities = comorbidities
        self.bmi=bmi
        self.simd=simd
        self.location=location
        self.__immune_response_igg = lambda x: 0
        self.__virus_response = lambda x: 0
        
        self.p_infected = self.calc_p_infected()

    def get_dob(self):
        return self.dob
        
    def calc_p_infected(self):
        age = self.get_age()
        location = self.location
        #less likely with age
        # - average at 50 
        #less likely with higher location (more rural)
        p = (1.5 - (age/150)) * (1/pow(location,0.2))# *simd term
        return p
        
    def get_p_infected(self):
        return self.p_infected
        
    def set_comorbidities(self,comorbidities):
        self.comorbidities = comorbidities
        
    def add_infection_record(self,date,variant):
        self.date_of_infections.append(date)
        self.infection_variants.append(variant)
        
    def set_vaccine_record(self,vaccines):
        if len(vaccines) == 0:
            return 
        self.date_of_vaccines,self.vaccine_products = zip(*vaccines)
        self.date_of_vaccines = list(self.date_of_vaccines)
        self.vaccine_products = list(self.vaccine_products)
        
        
    def get_vaccine_record(self):
        return self.date_of_vaccines

    def create_infection_response(self):
        if len(self.date_of_infections) == 0:
            return
        #import matplotlib.pyplot as plt

        functions = []
        loc = 5
        scale = 1
        for i,d in enumerate(self.date_of_infections):
            fn = spy.norm(loc=loc,scale=scale) 

            age = self.get_age(date=d)
            age_damp = (1 - (age/150)**2)
        
            a = 1 * age_damp * np.random.uniform(0.5,2)
            functions.append(
                lambda x,d=d,a=a,fn=fn:\
                    a*(1000)*fn.pdf((x-d).days)*(1 - (self.__immune_response_igg(x)+200)/1000) 
                    #a*(1000)*fn.pdf((x-d).days) 
                    #self.__immune_response_igg(x)
                    if x > d else 0
            )
        self.__virus_response = lambda x,f=functions: \
            sum([fn(x) for fn in f])
        
        #x = [datetime.date(2020,1,1) + datetime.timedelta(days=i) for i in range(1200)]
        #y = [self.__virus_response(i) for i in x]
        #plt.plot(x,y)
        #y = [self.__immune_response_igg(i) for i in x]
        #plt.plot(x,y)
        #plt.show();

    
    def create_immune_response(self):
        
        if len(self.date_of_vaccines) == 0:
            return lambda x: 0
        
        #loc = 7
        #s = 1.2
        #scale = 170
        s = 2.7
        scale = 150000
        loc = 1
        
        a_damp = np.prod([c.get_p_immuno() for c in self.comorbidities])
        a_damp = a_damp if a_damp > 0 else 0

        events = list(zip(self.date_of_vaccines,self.vaccine_products))
                
        igg_functions = []
        for i,(d,p) in enumerate(events):
            
            fn = spy.lognorm(loc=loc,s=s,scale=scale)#s=s+i*0.1,scale=scale+i*200) 
            a = 6000000 + 5000000*i
            
            
            age = self.get_age(date=d)
            age_damp = (1 - (age/150)**2)

            a *= a_damp * age_damp * m_igg_scale[p]
            
            igg_functions.append(
                lambda x,d=d,a=a,fn=fn:\
                    a*fn.pdf((x-d).days)
                    if x > d else 0
            )
        #add infection responses
        for i,d in enumerate(self.date_of_infections):
            fn = spy.lognorm(loc=loc,s=s,scale=scale)#s=s+i*0.1,scale=scale+i*200) 
            a = 600000 + 500000*i
            age = self.get_age(date=d)
            age_damp = (1 - (age/150)**2)
            a *= a_damp * age_damp * 0.01
            igg_functions.append(
                lambda x,d=d,a=a,fn=fn:\
                    a*fn.pdf((x-d).days)
                    if x > d else 0
            )
        
            
        self.__immune_response_igg = lambda x,f=igg_functions: \
                                     sum([fn(x) for fn in f])
        #self.__immune_response_functions = igg_functions
                
    def get_immune_response(self,x):
        return self.__immune_response_igg(x)
    
    def get_infection_response(self,x):
        return self.__virus_response(x)
        
    def get_age(self,date=datetime.date.today()):
        if self.age:
            return self.age
        return int((date - self.dob).days / 365.25)

    def run_outcomes(self,days=28):
        self.date_of_death = None
        cumres = 0
        symptoms = []
        days = [i for i in range(days)]
        for date in self.date_of_infections:
            if not self.date_of_death == None:
                break
            for i in days:# random.sample(days,28):
                d = date + datetime.timedelta(days=i)
                res = self.get_infection_response(d)
                cumres += res
                cumres -= np.random.normal(loc=0.2,scale=0.1)*cumres
                ires = self.get_immune_response(d)
                cumres -= np.random.normal(loc=0.2,scale=0.1)*ires
                s = get_symptoms(cumres)
                if len(s)>0:
                    #symptoms.append({'id':self.id,'date':d,'symptoms':s})
                    symptoms.append((d,s))
                    if 'death' in s:
                        self.date_of_death = d
                        break
                if cumres<0:
                    break
        self.outcomes = symptoms
        if self.date_of_death:
            self.date_of_infections = [
                d
                for d in self.date_of_infections
                if d > self.date_of_death
            ]
            self.date_of_vaccines = [
                d
                for d in self.date_of_vaccines
                if d > self.date_of_death
            ]


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    p = Study()
    p.set_vaccine_record([(datetime.datetime(2021,1,1),'Pf'),
                          (datetime.datetime(2021,5,1),'Pf'),
                          (datetime.datetime(2022,1,1),'Pf')])
    p.create_immune_response()
    for vd in p.date_of_vaccines:
        xdata = [
            d for i in range(1000)
            if (d := datetime.datetime(2020,12,1) + datetime.timedelta(days=i)) > vd
        ]
        ydata = [p.get_immune_response(x) for x in xdata]
        plt.plot(xdata,ydata)
    plt.show();
