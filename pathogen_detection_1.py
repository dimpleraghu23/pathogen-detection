import pandas as pd

from colorama import Fore, Style

# Colorful print function for better readability

def cprint(text, color):
    print(color + text + Style.RESET_ALL)

# Step 1: Load Data

cprint("Step 1: Loading Data", Fore.GREEN)
try:
    df = pd.read_csv('pathogen_detection.csv')
except FileNotFoundError:
    cprint("The CSV file was not found. Make sure it's in the correct directory.", Fore.MAGENTA) 
    exit()

# Display first 5 rows of the data

cprint("Displaying first 5 rows of the data:", Fore.MAGENTA)
print(df.head())

# Step 2: Data Cleaning

cprint("Step 2: Data Cleaning", Fore.GREEN)

# Remove rows with missing or erroneous Strains
df.dropna (subset=['Strain'], inplace=True)

# Step 3: Data Preprocessing
cprint("Step 3: Data Preprocessing", Fore.GREEN)

# Convert the string sequences to lists and then to DataFrame rows

sequence_list_series = df['Strain'].apply(list)

sequence_df = pd.DataFrame(sequence_list_series.tolist())

# One-Hot Encoding

cprint("Performing One-Hot Encoding...", Fore.MAGENTA)

stacked_sequence_df = sequence_df.stack().reset_index(level=1, drop=True)

onehot_encoded = pd.get_dummies(stacked_sequence_df, prefix='', prefix_sep='')

onehot_encoded.index.name = 'original_index'

onehot_df = onehot_encoded.groupby('original_index').sum()

# Display transformed DataFrame

cprint("Displaying the transformed DataFrame:", Fore.MAGENTA)

print(onehot_df.head())

# Step 4: Save the Processed Data

cprint("Step 4: Saving the Processed Data", Fore.GREEN)

onehot_df.to_csv('processed_strain_sequences.csv', index=False)

cprint("Data processing complete!", Fore. MAGENTA)

import matplotlib.pyplot as plt

import seaborn as sns

cprint("Step 5: Data Visualisation", Fore.GREEN)

# Get the top 10 most common locations

top_locations = df[ 'Location'].value_counts().head(10)

print(top_locations.head())

# Get the top 10 most common isolation sources

top_isolation_sources = df['Isolation source'].value_counts().head(10) 
print(top_isolation_sources.head())

# Get the top 10 most common serovars

top_serovars = df['Serovar'].value_counts().head(10) 
print(top_serovars.head())

# Plot the top 10 most common locations

plt.figure(figsize=(10, 6)) 
sns.barplot(y=top_locations.index, x=top_locations.values, palette='viridis') 
plt.title('Top 10 Most Common Locations') 
plt.xlabel('Count') 
plt.ylabel('Location') 
plt.show()
# Plot the top 10 most common isolation sources

plt.figure(figsize=(10, 6))

sns.barplot(y=top_isolation_sources.index, x=top_isolation_sources.values, palette='viridis')

plt.title('Top 10 Most Common Isolation Sources')

plt.xlabel('Count')

plt.ylabel('Isolation Source')

plt.show()

# Plot the top 10 most common serovars

plt.figure(figsize=(10, 6))

sns.barplot(y=top_serovars.index, x=top_serovars.values, palette='viridis')

plt.title('Top 10 Most Common Serovars')

plt.xlabel('Count')

plt.ylabel('Serovar')

plt.show()

#Get the top 10 most common AMR genotypes

top_amr_genotypes = df[ 'AMR genotypes'].value_counts().head(10)

print(top_amr_genotypes.head())

# Get the top 10 most common computed types

top_computed_types = df['Computed types'].value_counts().head (10)

print(top_computed_types.head())

# Plot the top 10 most common AMR genotypes

plt.figure(figsize=(10, 6))

sns.barplot(y=top_amr_genotypes.index, x=top_amr_genotypes.values, palette='viridis')

plt.title('Top 10 Most Common AMR Genotypes')

plt.xlabel('Count')

plt.ylabel('AMR Genotypes')

plt.show()

# Plot the top 10 most common computed types

plt.figure(figsize=(10, 6))

sns.barplot(y=top_computed_types.index, x=top_computed_types.values, palette='viridis')

plt.title('Top 10 Most Common Computed Types')

plt.xlabel('Count')

plt.ylabel('Computed Types')

plt.show()
# Convert 'Create date' column to datetime

df['Create date'] = pd.to_datetime(df['Create date'])

# Extract year from 'Create date'

df[ 'Year'] = df['Create date'].dt.year

# Get the number of samples collected each year

samples_per_year = df[ 'Year'].value_counts().sort_index()

print(samples_per_year)

# Plot the number of samples collected each year

plt.figure(figsize=(12, 6))

sns.lineplot(x=samples_per_year.index, y=samples_per_year.values, marker='o')

plt.title('Number of Samples Collected Each Year')

plt.xlabel('Year')

plt.ylabel('Number of Samples')

plt.grid(True)

plt.show()
# Get the top 10 most common SNP clusters

top_snp_clusters = df[ 'SNP cluster'].value_counts().head(10)

print(top_snp_clusters)

# Get the distribution of the top 10 most common SNP clusters across different locations

snp_clusters_locations = df[df[ 'SNP cluster'].isin(top_snp_clusters.index)]

snp_clusters_locations = snp_clusters_locations.groupby(['SNP cluster', 'Location']).size().unstack().fillna(0)

print(snp_clusters_locations)

# Get the distribution of the top 10 most common SNP clusters across different serovars

snp_clusters_serovars = df[df['SNP cluster'].isin(top_snp_clusters.index)]

snp_clusters_serovars = snp_clusters_serovars.groupby(['SNP cluster', 'Serovar']).size().unstack().fillna(0) 
print(snp_clusters_serovars)

# Plot the top 10 most common SNP clusters

plt.figure(figsize=(10, 6))

sns.barplot(y=top_snp_clusters.index, x=top_snp_clusters.values, palette='viridis')

plt.title('Top 10 Most Common SNP Clusters')

plt.xlabel('Count')

plt.ylabel('SNP Cluster')

plt.show()
#Get the distribution of different isolation types

isolation_types = df['Isolation type'].value_counts()

print(isolation_types)

# Plot the distribution of different isolation types

plt.figure(figsize=(8, 6))

sns.barplot(y=isolation_types.index, x=isolation_types.values, palette='viridis')

plt.title('Distribution of Different Isolation Types')

plt.xlabel('Count')

plt.ylabel('Isolation Type')

plt.show()

# Get the distribution of 'Min-same' and 'Min-diff' columns

min_same = df['Min-same'].dropna()

print(min_same)

min_diff = df[ 'Min-diff'].dropna()

print(min_diff)

# Plot the distribution of 'Min-same' and 'Min-diff' columns

plt.figure(figsize=(12, 6))

sns.histplot(min_same, kde= True, color='b', label='Min-same')

sns.histplot(min_diff, kde=True, color='r', label='Min-diff')

plt.title('Distribution of Min-same and Min-diff')

plt.xlabel('Value')

plt.ylabel('Density')

plt.legend()

plt.show()

from sklearn.cluster import MiniBatchKMeans

from sklearn.decomposition import IncrementalPCA

from sklearn.metrics import silhouette_score

from sklearn.preprocessing import StandardScaler

from sklearn.utils import shuffle

def cprint(text, color):
    print(color + text + Style.RESET_ALL)

# Feature Engineering

cprint("Feature Engineering", Fore.GREEN)

# Sampling a subset of the data for quick experiments

onehot_df_sampled = shuffle(onehot_df).sample(frac=0.1, random_state=42)

# Scaling the features

scaler = StandardScaler()

scaled_onehot_df = scaler.fit_transform(onehot_df_sampled)

# Dimensionality reduction using Incremental PCA

incremental_pca = IncrementalPCA(n_components=2, batch_size=200)

principal_components = incremental_pca.fit_transform(scaled_onehot_df)
# Model Building

cprint("Model Building", Fore.GREEN)

# K-means clustering with MiniBatchKMeans

mini_batch_kmeans = MiniBatchKMeans (n_clusters=3, batch_size=200)

kmeans_clusters = mini_batch_kmeans.fit_predict(principal_components)

# Model Evaluation

cprint("Model Evaluation", Fore.GREEN)

# Evaluate clustering by silhouette score

kmeans_score = silhouette_score (principal_components, kmeans_clusters)

print (f"K-means silhouette score: {kmeans_score}")

cprint("Model building complete!", Fore.MAGENTA)

from sklearn.model_selection import train_test_split

from sklearn.ensemble import RandomForestClassifier

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

# Step 7: Machine Learning Algorithm (MLA) Design

cprint("Step 6: Machine Learning Algorithm (MLA) Design", Fore.GREEN)

# Splitting the data into training and testing sets

X_train, X_test, y_train, y_test = train_test_split(principal_components, kmeans_clusters, test_size=0.2, random_state=42 )

# Initialize a Random ForestClassifier (you can choose a different classifier)

clf = RandomForestClassifier(random_state=42)

# Training the classifier 
clf.fit(X_train, y_train)

# Predicting labels on the test data 
y_pred = clf.predict(X_test)
# Step 8: Model Evaluation

cprint("Step 7: Model Evaluation", Fore.GREEN)

# Evaluate the model's performance

accuracy = accuracy_score(y_test, y_pred)

precision = precision_score(y_test, y_pred, average='weighted')

recall = recall_score(y_test, y_pred, average='weighted')

f1 = f1_score(y_test, y_pred, average='weighted')

print (f"Accuracy: {accuracy:.2f}")
print(f"Precision: {precision:.2f}")
print (f"Recall: {recall:.2f}")
print(f"F1 Score: {f1:.2f}")
cprint("MLA design, training, and evaluation complete!", Fore.MAGENTA)
