from scipy import stats
import numpy as np

#obs=np.array([[35,55-35], [89,356-89]])
obs=np.array([[4,89-4], [4,17-4]])
#chi,p,df,ex=
p=stats.chi2_contingency(obs)
print(p)


if __name__ == '__main__':
    pass
