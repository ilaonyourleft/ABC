from sklearn.datasets import make_blobs
import matplotlib.pyplot as plt
import pandas as pd


features, clusters = make_blobs(n_samples=10000, n_features=2, centers=5, cluster_std=0.2, shuffle=True)

print("feature matrix: ")
print(pd.DataFrame(features, columns=["Feature 1", "Feature 2"]).head())

df = pd.DataFrame(features)
df.to_csv('dataset10000.csv', sep=',', index=False, encoding='utf-8')

plt.scatter(features[:, 0], features[:, 1])
plt.show()

