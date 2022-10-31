import scipy.stats as spy
import numpy as np

class Comorbidity:
    def __init__(self,p_immuno=1,p_severe=1):
        self.p_immuno = spy.norm(loc=p_immuno,scale=0.2)
        self.p_severe = spy.norm(loc=p_severe,scale=0.2)

    def get_name(self):
        return type(self).__name__
    
    def get_p_immuno(self):
        return self.p_immuno.rvs()
    
    def get_p_severe(self):
        return self.p_severe.rvs()

class BloodCancer(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=0.1,p_severe=0.8)
        
    def get_p_immuno(self):
        return super().get_p_immuno()*np.random.choice([0,1], p=[0.3,0.7])

class ImmunoDeficiency(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=np.random.uniform(0,0.5),p_severe=0.8)

class COPD(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=0.6,p_severe=0.7)

class CHD(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=1,p_severe=0.75)
    
class Epilepsy(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=1,p_severe=0.98)
    
class HeartFailure(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=1,p_severe=0.87)

class Thrombosis(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=1,p_severe=0.92)
        
class Asthma(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=1,p_severe=0.75)

class RespiratoryCancer(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=0.7,p_severe=0.5)

class Diabetes(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=0.9,p_severe=0.9)
        
class CKD(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=0.95,p_severe=0.7)

class MentalHealth(Comorbidity):
    def __init__(self):
        super().__init__(p_immuno=1,p_severe=1)

        
def get_comorbidities():
    return {
        k:x
        for k,x in globals().items()
        if isinstance(x,type) and issubclass(x,Comorbidity)
    }
