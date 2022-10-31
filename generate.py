import seave
import datetime
import scipy.stats as spy
import matplotlib.pyplot as plt
import numpy as np
import dill as pickle
import cProfile
import pstats
import pandas as pd

#p = seave.Study()
#print (p.get_age())

variants = {
    'Alpha':(datetime.date(2020,3,1),datetime.date(2020,7,1),1),
    'Beta':(datetime.date(2020,6,1),datetime.date(2020,11,1),1.5),
    'Delta':(datetime.date(2020,10,15),datetime.date(2021,4,1),1.1),
    'Kappa':(datetime.date(2021,3,1),datetime.date(2022,1,1),0.4)
}



covid_pandemic = seave.Pandemic(datetime.date(2020,3,1),variants=variants)

class Population:
    def __init__(self):
        self.comorbidity_map = { 
            seave.comorbidities.BloodCancer(): lambda age,bmi : 0.01*pow((age/100),0.3),
            seave.comorbidities.COPD(): lambda age,bmi: 0.02*(age/100)*pow((bmi/20),0.8),
            seave.comorbidities.Asthma(): lambda age,bmi: 0.1*pow((age/100),0.6)*pow((bmi/20),0.7),
            seave.comorbidities.RespiratoryCancer(): lambda age,bmi: 0.005*(age/100),
            seave.comorbidities.Diabetes(): lambda age,bmi: 0.3*(age/100)*(bmi/20),
            seave.comorbidities.CKD(): lambda age,bmi: 0.003*(age/100),
            seave.comorbidities.ImmunoDeficiency(): lambda age,bmi: 0.004,
            seave.comorbidities.MentalHealth(): lambda age,bmi: 0.10*pow((age/100),0.6),
            seave.comorbidities.CHD(): lambda age,bmi: 0.003*pow((age/100),1.5),
            seave.comorbidities.Epilepsy(): lambda age,bmi: 0.05*pow((age/100),0.6)*pow((bmi/20),1.1),
            seave.comorbidities.HeartFailure(): lambda age,bmi: 0.10*pow((age/100),2.2),
            seave.comorbidities.Thrombosis(): lambda age,bmi: 0.002*pow((age/100),0.6)
        }
        
    def age_distro(self,turn,end):
        totarea = turn + (end-turn)/2  
        areauptoturn = turn            
        areasloped = (end-turn)/2     
        p1 = areauptoturn/totarea
        c = np.random.choice([True,False],p=[p1,1-p1]) 
        if c:
            return np.random.uniform(low=0,high=turn)
        else:
            return np.random.triangular(left=turn,mode=turn,right=end)

    fn_age = lambda self : self.age_distro(turn=50,end=100)
    fn_sex = lambda self : np.random.choice(['M','F'],p=[0.49,0.51])
    fn_simd = lambda self: np.random.randint(1,6)
    fn_location = lambda self: np.random.choice([1,2,3],p=[0.65,0.25,0.1])

    #make these age dependent 
    def fn_comorbidites(self,age,bmi):
        return [
            c
            for c,p in self.comorbidity_map.items()
            if np.random.uniform() < p(age,bmi)
        ]
        
    fn_bmi = lambda self: spy.norm(loc=20,scale=7).rvs()
    
    def generate(self,n):
        retval = []
        for i in range(n):
            age = self.fn_age()
            sex = self.fn_sex()
            bmi = self.fn_bmi()
            comorbidities = self.fn_comorbidites(age,bmi)
            retval.append(
                {
                    'age':age,
                    'sex':sex,
                    'bmi':bmi,
                    'simd':self.fn_simd(),
                    'location':self.fn_location(),
                    'comorbidities': comorbidities
                })
        return retval

get_fdelay = lambda age: (1 - age/100) if age < 100 else 1
get_type = lambda : np.random.choice(['AZ','Pfizer','Moderna'],p=[0.6,0.3,0.1])

def generate_vaccines(p):
    vaccines = []
    vaccinated = np.random.uniform()>0.05 #95% of people
    if not vaccinated:
        return vaccines

    #do at least one vaccine
    date_first_vaccine = datetime.date(2021,2,1)

    #work out the delay from the start of the vaccine period
    fdelay = get_fdelay(p.get_age()) 
    days_delay = int(np.random.normal(loc=fdelay*200,scale=20))
    #generate a random vaccine product
    _type = get_type()
    vaccines.append((date_first_vaccine + datetime.timedelta(days=days_delay),_type))

    #up to 3 additional vaccines
    for iv in range(3):
        #randomly limit the number of vaccines someone got
        #make it less likely with more vaccines
        if np.random.uniform() < 0.10 + iv*0.2:
            break
        days_delay += int(np.random.normal(loc=80+10*iv,scale=10))
        _type = get_type()
        vaccines.append((date_first_vaccine + datetime.timedelta(days=days_delay),_type))

    return vaccines

def get_pcr_tests(p):
    return [ 
        {'id':p.id,'date':d,'result':int(p.get_infection_response(d)>5)}
        for x in p.date_of_infections
        if (d:=x+datetime.timedelta(days=np.random.randint(0,5)) )
    ]

def get_vaccines(p):
    return [ 
        {'id':p.id,'date':d,'product':prod}
        for d,prod in zip(p.date_of_vaccines,p.vaccine_products)
    ]

def get_serology_tests(p,start=datetime.date(2020,6,1),end=800,n=2):
    #def get_dose(p,d):
    #    dose = 0
    #    for i,vd in enumerate(p.date_of_vaccines):
    #        if d > vd:
    #            dose = i
    #        else:
    #            break
    #    return dose
    return [ 
        {'id':p.id,'date':d,'IgG':p.get_immune_response(d)}#,'dose':get_dose(p,d)}
        for _ in range(np.random.randint(0,n))
        if (d:=start+datetime.timedelta(days=np.random.randint(0,end)))
    ]
    

#with cProfile.Profile() as pr:
pop = Population()
demographics = pop.generate(1000)
#    stats = pstats.Stats(pr)
#    stats.sort_stats(pstats.SortKey.TIME)
#    stats.print_stats() # Print The Stats


cohort = []
for ii,demo in enumerate(demographics):
    age = demo.pop('age')
    demo['dob'] = seave.Study.age_to_dob(age)
    p = seave.Study(**demo)
    p.set_vaccine_record(generate_vaccines(p))
    #p.create_immune_response()
    
    for i in range(0,1000):
        x = datetime.date(2020,1,1) + datetime.timedelta(days=i)
        #get an infection rate probability, scale this with the person's individual
        #infection rate modifier (i.e. higher for younger people)
        pinfection = covid_pandemic.get_infection_rate(x)*p.get_p_infected()
        if pinfection > np.random.uniform():
            variants = covid_pandemic.get_variants(x) 
            variant = np.random.choice(list(variants.keys()),p=list(variants.values()))
            p.add_infection_record(x,variant)

    p.create_infection_response()
    p.create_immune_response()
    p.run_outcomes()

    
    df_demo = pd.DataFrame([p.get_demo()])
    df_demo.to_csv('raw_data/demographics.csv',
                   mode='w' if ii==0 else 'a',
                   index=False,
                   header=True if ii==0 else False)
    df_pcr = pd.DataFrame(get_pcr_tests(p))
    df_pcr.to_csv('raw_data/pcr.csv',
                   mode='w' if ii==0 else 'a',
                   index=False,
                   header=True if ii==0 else False)


    df_vaccines = pd.DataFrame(get_vaccines(p))
    df_vaccines.to_csv('raw_data/vaccines.csv',
                   mode='w' if ii==0 else 'a',
                   index=False,
                   header=True if ii==0 else False)


    
    df_serology = pd.DataFrame(get_serology_tests(p,n=10))#.sort_values('date').reset_index(drop=True).set_index('id')
    df_serology.to_csv('raw_data/serology.csv',
                       mode='w' if ii==0 else 'a',
                       index=False,
                       header=True if ii==0 else False)
    
    df_outcomes = pd.DataFrame([{'id':p.id,'date':k[0],'symptoms':k[1]} for k in p.outcomes])
    df_outcomes.to_csv('raw_data/outcomes.csv',
                  mode='w' if ii==0 else 'a',
                  index=False,
                  header=True if ii==0 else False)

    if p.date_of_death:
        df_death = pd.DataFrame([{'id':p.id,'date':p.date_of_death}])
        df_death.to_csv('raw_data/pcr.csv',
                   mode='a',
                   index=False,
                   header=False)

        
    
#with open(f"people/cohort0.pickle","wb") as f:
#    pickle.dump(cohort,f)

    
