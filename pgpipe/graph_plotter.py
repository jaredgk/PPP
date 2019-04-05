import pandas as pd
import matplotlib
matplotlib.use('agg')
import subprocess
from io import StringIO
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Color blind color cycle hex codes
CB_color_cycle = ['#ffb000', '#dc267f', '#648fff', '#fe6100', '#785ef0']


def manhattan_plot(in_file):
    try:
        data = pd.read_table(in_file, sep="\t")
        data = data.dropna()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        data.plot(kind='scatter', x='POS', y='WEIR_AND_COCKERHAM_FST', ax=ax)
        ax.set_xlabel('Position')
        fig.savefig(in_file + ".pdf", bbox_inches='tight')
    except FileNotFoundError as e:
        print(e+' not found')


def admixture_barplot(in_file, numinds):
    matrix = []
    try:
        with open(in_file, 'r') as file:
            for line in file:
                if 'Miss' in line:
                    for i in range(numinds):
                        row = (next(file).strip().replace(':', '').split())
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
    except FileNotFoundError as e:
        print(e)


def pca_plot(in_file):
    try:
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
    except FileNotFoundError as e:
        print(e)


# Function for plotting graphs for admixture output data set
def bar_plot(in_file):
    try:
        a = in_file
        df = pd.read_table(a, header=None, delim_whitespace=True)
        df.plot.bar(stacked=True, figsize=(23, 5), color=CB_color_cycle)
        plt.savefig(in_file + ".pdf", bbox_inches='tight')
    except FileNotFoundError as e:
        print(e)


# Function for plotting graphs for ima2p output data sets
def do_plots(in_file):
    #try:
        output1 = subprocess.run("grep -B 1000 -m 1 \"SumP\" " + in_file + " | grep -v \"SumP\"", shell=True,
                                 stdout=subprocess.PIPE, universal_newlines=True)
        output2 = subprocess.run("grep -B 1000 -m 2 \"SumP\" " + in_file + " | grep -v \"SumP\" | tail -n1000",
                                 shell=True, stdout=subprocess.PIPE, universal_newlines=True)

        t0 = pd.read_table(StringIO(output1.stdout), sep="\t", header=None)
        t0.drop(t0.columns[[0]], axis=1, inplace=True)

        q0 = pd.read_table(StringIO(output2.stdout), sep="\t", header=None)
        q0.drop(q0.columns[[0]], axis=1, inplace=True)

        fig, ax = plt.subplots(nrows=3, ncols=2)

        ax[0, 0].plot(t0[2], t0[3])
        ax[0, 0].set_xlabel("Divergence times (t0)")
        ax[0, 0].set_ylabel("Frequency")

        ax[0, 1].plot(q0[2], q0[3])
        ax[0, 1].set_xlabel("Population size (q0)")
        ax[0, 1].set_ylabel("Frequency")

        ax[1, 0].plot(q0[4], q0[5])
        ax[1, 0].set_xlabel("Population size (q1)")
        ax[1, 0].set_ylabel("Frequency")

        ax[1, 1].plot(q0[6], q0[7])
        ax[1, 1].set_xlabel("Population size (q2)")
        ax[1, 1].set_ylabel("Frequency")

        ax[2, 0].plot(q0[8], q0[9])
        ax[2, 0].set_xlabel("Migration rate (m0->1)")
        ax[2, 0].set_ylabel("Frequency")

        ax[2, 1].plot(q0[10], q0[11])
        ax[2, 1].set_xlabel("Migration rate (m1->0)")
        ax[2, 1].set_ylabel("Frequency")

        plt.tight_layout()
        plt.savefig(in_file + ".pdf")
    #except:
    #    print(in_file + ": No such file or directory")
