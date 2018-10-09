import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def manhattan_plot(in_file):
    data = pd.read_table(in_file, sep="\t")
    data = data.dropna()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    data.plot(kind='scatter', x='POS', y='WEIR_AND_COCKERHAM_FST', ax=ax)
    ax.set_xlabel('Position')
    fig.savefig(in_file + ".pdf", bbox_inches='tight')


def admixture_barplot(in_file, numinds):
    matrix = []
    with open(in_file, 'r') as f:
        for line in f:
            if 'Miss' in line:
                for i in range(numinds):
                    row = (next(f).strip().replace(':', '').split())
                    matrix.append(row)
    data = pd.DataFrame(matrix, columns=['Index', 'Label', '(%Miss)', 'Pop', 'A', 'B'])
    data = data.drop('Index',  axis=1)

    fig = plt.figure()
    r = data.index.tolist()

    # convert series A and B dataset to list
    series_a = data[data.columns[3]].apply(pd.to_numeric).tolist()
    series_b = data[data.columns[4]].apply(pd.to_numeric).tolist()
    # Create bar plot for series A and B data
    plt.bar(r, series_a,  width=1.0, color='#b5ffb9')
    plt.bar(r, series_b, bottom=series_a, width=1.0, color='#f9bc86')
    # Add label
    plt.xlabel("Label")
    plt.ylabel("Frequency")
    # Save plot as pdf file
    fig.savefig(in_file + ".pdf", bbox_inches='tight')


def pca_plot(in_file):
    df = pd.read_table(in_file, index_col=0, header=0)

    features = list(df.columns.values)
    # Separating out the features
    data = df.loc[:, features].values

    # Standardizing the features
    data = StandardScaler().fit_transform(data)
    pca = PCA(n_components=2)
    prncpl_comp = pca.fit_transform(data)
    prncpl_df = pd.DataFrame(data=prncpl_comp, columns=['PC1', 'PC2'])
    fig = plt.figure()
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('2 component PCA', fontsize=20)

    plt.scatter(prncpl_df['PC1'], prncpl_df['PC2'])
    plt.legend('Samples')
    plt.grid()
    fig.savefig(in_file + ".pdf", bbox_inches='tight')
