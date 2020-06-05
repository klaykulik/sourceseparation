import pymc3 as pm

with pm.Model() as model:
    parameter = pm.Exponential("poisson_param", 1.0)
    data_generator = pm.Poisson("data_generator", parameter)



with model:
    data_plus_one = data_generator + 1


with pm.Model() as model:
    theta = pm.Exponential("theta", 2.0)
    data_generator = pm.Poisson("data_generator", theta)




with pm.Model() as ab_testing:
    p_A = pm.Uniform("P(A)", 0, 1)
    p_B = pm.Uniform("P(B)", 0, 1)


print("parameter.tag.test_value =", parameter.tag.test_value)
print("data_generator.tag.test_value =", data_generator.tag.test_value)
print("data_plus_one.tag.test_value =", data_plus_one.tag.test_value)

with pm.Model() as model:
    parameter = pm.Exponential("poisson_param", 1.0, testval=0.5)

print("\nparameter.tag.test_value =", parameter.tag.test_value)



with pm.Model() as model:
    lambda_1 = pm.Exponential("lambda_1", 1.0)
    lambda_2 = pm.Exponential("lambda_2", 1.0)
    tau = pm.DiscreteUniform("tau", lower=0, upper=10)

new_deterministic_variable = lambda_1 + lambda_2


import numpy as np

n_data_points = 5  # in CH1 we had ~70 data points
idx = np.arange(n_data_points)
with model:
    lambda_ = pm.math.switch(tau >= idx, lambda_1, lambda_2)


import theano.tensor as tt

with pm.Model() as theano_test:
    p1 = pm.Uniform("p", 0, 1)
    p2 = 1 - p1
    p = tt.stack([p1, p2])

    assignment = pm.Categorical("assignment", p)



from IPython.core.pylabtools import figsize
import matplotlib.pyplot as plt
import scipy.stats as stats
# figsize(12.5, 4)


# samples = lambda_1.random(size=20000)
# plt.hist(samples, bins=70, normed=True, histtype="stepfilled")
# plt.title("Prior distribution for $\lambda_1$")
# plt.xlim(0, 8);
# plt.show()

# data = np.array([10, 5])
# with model:
#     fixed_variable = pm.Poisson("fxd", 1, observed=data)
# print("value: ", fixed_variable.tag.test_value)


# # We're using some fake data here
# data = np.array([10, 25, 15, 20, 35])
# with model:
#     obs = pm.Poisson("obs", lambda_, observed=data)
# print(obs.tag.test_value)


# tau = np.random.randint(0, 80)
# print(tau)

# alpha = 1. / 20.
# lambda_1, lambda_2 = np.random.exponential(scale=1 / alpha, size=2)
# print(lambda_1, lambda_2)

# data = np.r_[stats.poisson.rvs(mu=lambda_1, size=tau), stats.poisson.rvs(mu=lambda_2, size=80 - tau)]

# plt.bar(np.arange(80), data, color="#348ABD")
# plt.bar(tau - 1, data[tau - 1], color="r", label="user behaviour changed")
# plt.xlabel("Time (days)")
# plt.ylabel("count of text-msgs received")
# plt.title("Artificial dataset")
# plt.xlim(0, 80)
# plt.legend()
# plt.show()


def plot_artificial_sms_dataset():
    tau = stats.randint.rvs(0, 80)
    alpha = 1./20.
    lambda_1, lambda_2 = stats.expon.rvs(scale=1/alpha, size=2)
    data = np.r_[stats.poisson.rvs(mu=lambda_1, size=tau), stats.poisson.rvs(mu=lambda_2, size=80 - tau)]
    plt.bar(np.arange(80), data, color="#348ABD")
    plt.bar(tau - 1, data[tau-1], color="r", label="user behaviour changed")
    plt.xlim(0, 80);

figsize(12.5, 5)
plt.title("More example of artificial datasets")
for i in range(4):
    plt.subplot(4, 1, i+1)
    plot_artificial_sms_dataset()

plt.show()